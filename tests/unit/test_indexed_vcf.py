"""Tests for the indexed remote-VCF transport used by laptop-safe acquisition.

Everything here is offline: BGZF blocks and a tabix index are synthesized in
memory so the block decoder, the index math, and the transfer budget are all
exercised without touching the network. The one live-source test is marked
``network`` and deselected by default.
"""

from __future__ import annotations

import gzip
import struct
import zlib

import pytest

from canidae.io.indexed_vcf import (
    BGZF_MAX_BLOCK,
    SOURCE_VCF_URL,
    TransferBudget,
    chunks_for,
    decompress_bgzf,
    genotype,
    merge_chunks,
    parse_tabix,
    reg2bins,
    source_header,
)


def make_bgzf_block(payload: bytes) -> bytes:
    """A single spec-conformant BGZF block (gzip member with the BC extra field)."""
    compressor = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = compressor.compress(payload) + compressor.flush()
    total = 18 + len(body) + 8
    header = struct.pack(
        "<BBBBIBBHBBHH",
        0x1F, 0x8B, 8, 4,  # magic, deflate, FEXTRA
        0, 0, 0xFF,        # mtime, xfl, os
        6, 66, 67, 2,      # XLEN, SI1='B', SI2='C', SLEN
        total - 1,         # BSIZE at byte offset 16
    )
    trailer = struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))
    return header + body + trailer


def make_tabix(names: list[str], bins: dict[int, list[tuple[int, int]]],
               linear: list[int]) -> bytes:
    """A minimal single-reference TBI index, gzip-compressed as tabix stores it."""
    raw = bytearray(b"TBI\x01")
    raw += struct.pack("<i", len(names))
    raw += struct.pack("<6i", 2, 1, 2, 0, ord("#"), 0)
    encoded = b"".join(name.encode() + b"\x00" for name in names)
    raw += struct.pack("<i", len(encoded)) + encoded
    for _ in names:
        raw += struct.pack("<i", len(bins))
        for bin_id, chunks in bins.items():
            raw += struct.pack("<I", bin_id) + struct.pack("<i", len(chunks))
            for begin, end in chunks:
                raw += struct.pack("<QQ", begin, end)
        raw += struct.pack("<i", len(linear))
        raw += struct.pack(f"<{len(linear)}Q", *linear) if linear else b""
    return gzip.compress(bytes(raw))


# -- transfer budget -------------------------------------------------------------------


def test_transfer_budget_tracks_and_refuses_overrun() -> None:
    budget = TransferBudget(limit=1000)
    budget.add(400)
    assert budget.used == 400
    budget.reserve(600)  # exactly at the limit is allowed
    with pytest.raises(RuntimeError, match="safety limit"):
        budget.reserve(601)
    with pytest.raises(RuntimeError, match="safety limit"):
        budget.add(-1)
    assert budget.used == 400  # a refused reservation must not consume budget


# -- BGZF ------------------------------------------------------------------------------


def test_decompress_bgzf_concatenates_blocks() -> None:
    data = make_bgzf_block(b"first;") + make_bgzf_block(b"second;")
    assert decompress_bgzf(data) == b"first;second;"


def test_decompress_bgzf_honours_the_first_block_offset() -> None:
    data = make_bgzf_block(b"header\nbody\n")
    assert decompress_bgzf(data, first_uncompressed_offset=7) == b"body\n"


def test_decompress_bgzf_stops_at_a_truncated_trailing_block() -> None:
    """A byte range rarely ends on a block boundary; the tail must be dropped."""
    complete = make_bgzf_block(b"kept")
    truncated = make_bgzf_block(b"discarded")[:-4]
    assert decompress_bgzf(complete + truncated) == b"kept"


def test_decompress_bgzf_stops_at_non_gzip_bytes() -> None:
    assert decompress_bgzf(make_bgzf_block(b"kept") + b"\x00\x01\x02") == b"kept"


# -- tabix index math ------------------------------------------------------------------


def test_reg2bins_covers_every_level_for_a_point_query() -> None:
    bins = reg2bins(0, 1)
    assert bins[0] == 0
    # One bin per level: 1, 9, 73, 585, 4681 are the level start offsets.
    assert {1, 9, 73, 585, 4681} <= set(bins)


def test_parse_tabix_round_trip() -> None:
    bins = {4681: [(100 << 16, 200 << 16)]}
    names, refs = parse_tabix(make_tabix(["chr1"], bins, [0]))
    assert names == ["chr1"]
    assert len(refs) == 1
    assert refs[0].bins[4681] == [(100 << 16, 200 << 16)]
    assert refs[0].linear == [0]


def test_parse_tabix_rejects_a_foreign_index() -> None:
    with pytest.raises(RuntimeError, match="not a tabix"):
        parse_tabix(gzip.compress(b"CSI\x01padding"))


def test_chunks_for_drops_chunks_below_the_linear_offset() -> None:
    _, refs = parse_tabix(
        make_tabix(
            ["chr1"],
            {4681: [(10 << 16, 20 << 16), (200 << 16, 300 << 16)]},
            [100 << 16],
        )
    )
    # The linear index says nothing before virtual offset 100<<16 can match.
    assert chunks_for(refs[0], 1) == [(200 << 16, 300 << 16)]


def test_merge_chunks_joins_adjacent_ranges_and_splits_distant_ones() -> None:
    near_end = (BGZF_MAX_BLOCK // 2) << 16
    merged = merge_chunks({(0, near_end), (near_end, (2 * BGZF_MAX_BLOCK) << 16)})
    assert merged == [(0, (2 * BGZF_MAX_BLOCK) << 16)]

    far = (10 * BGZF_MAX_BLOCK) << 16
    assert len(merge_chunks({(0, 1 << 16), (far, far + (1 << 16))})) == 2


# -- genotype parsing ------------------------------------------------------------------


def test_genotype_accepts_biallelic_calls_only() -> None:
    assert genotype("0/1:30:99", 0) == "0/1"
    assert genotype("1|1", 0) == "1|1"
    assert genotype("./.", 0) is None       # missing
    assert genotype("0/2", 0) is None       # multiallelic
    assert genotype("0", 0) is None         # haploid
    assert genotype("30:99", 5) is None     # GT column absent


# -- live source (opt-in) --------------------------------------------------------------


@pytest.mark.network
@pytest.mark.slow
def test_live_source_header_is_range_fetched_within_budget() -> None:
    """Byte-range fetch against the real NHGRI source; run with `pytest -m network`."""
    budget = TransferBudget(limit=32 * 1024 * 1024)
    meta, columns = source_header(SOURCE_VCF_URL, budget)
    assert columns[:9] == [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    ]
    assert len(columns) - 9 == 722  # the 722-genome panel
    assert any(line.startswith("##fileformat=VCF") for line in meta)
    # The whole point of the indexed path: a header costs megabytes, not gigabytes.
    assert budget.used < 32 * 1024 * 1024
