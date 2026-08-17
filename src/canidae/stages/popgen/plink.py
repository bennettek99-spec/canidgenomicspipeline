"""Minimal PLINK 1 binary (.bed/.bim/.fam) writer.

Written so the ADMIXTURE-binary backend is self-contained (ADMIXTURE reads a PLINK
fileset) without depending on the ``plink`` executable to produce it. Implements the
SNP-major .bed layout: a 3-byte magic header then, per variant, 2 bits per sample packed
4-samples-per-byte, LSB first.

Genotype ALT-count encoding:  0 -> 00 (A1/A1), 1 -> 10 (het), 2 -> 11 (A2/A2),
missing -> 01. A2 is the ALT allele.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from canidae.stages.popgen.store import Genotypes

_MAGIC = bytes([0x6C, 0x1B, 0x01])  # PLINK .bed magic + SNP-major mode
# ALT allele count -> 2-bit code
_ENCODE = {0: 0b00, 1: 0b10, 2: 0b11, -1: 0b01}


@dataclass(slots=True)
class PlinkFileset:
    bed: Path
    bim: Path
    fam: Path


def write_plink_bed(genotypes: Genotypes, prefix: Path) -> PlinkFileset:
    """Write ``{prefix}.bed/.bim/.fam`` from a :class:`Genotypes`. Returns the paths."""
    prefix = Path(prefix)
    n_alt = np.asarray(genotypes.calls.to_n_alt(fill=-1), dtype=np.int8)  # (n_var, n_samp)
    bed = prefix.with_suffix(".bed")
    bed.write_bytes(_MAGIC + _pack_snp_major(n_alt))
    _write_bim(prefix.with_suffix(".bim"), genotypes)
    _write_fam(prefix.with_suffix(".fam"), genotypes)
    return PlinkFileset(bed=bed, bim=prefix.with_suffix(".bim"),
                        fam=prefix.with_suffix(".fam"))


def _pack_snp_major(n_alt: np.ndarray) -> bytes:
    n_var, n_samp = n_alt.shape
    codes = np.select(
        [n_alt == 0, n_alt == 1, n_alt == 2],
        [_ENCODE[0], _ENCODE[1], _ENCODE[2]],
        default=_ENCODE[-1],
    ).astype(np.uint8)
    pad = (-n_samp) % 4
    if pad:
        codes = np.pad(codes, ((0, 0), (0, pad)), constant_values=_ENCODE[-1])
    quads = codes.reshape(n_var, -1, 4)
    weights = np.array([1, 4, 16, 64], dtype=np.uint16)
    packed = (quads * weights).sum(axis=2).astype(np.uint8)
    return packed.tobytes()


def _write_bim(path: Path, g: Genotypes) -> None:
    chrom_codes = _numeric_chrom(g.chrom)
    lines = [
        f"{chrom_codes[i]}\t{g.chrom[i]}:{g.pos[i]}\t0\t{int(g.pos[i])}\tA\tG"
        for i in range(g.n_variants)
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_fam(path: Path, g: Genotypes) -> None:
    lines = [f"{s}\t{s}\t0\t0\t0\t-9" for s in g.samples]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _numeric_chrom(chrom: np.ndarray) -> list[int]:
    """Map contig names to integer codes ADMIXTURE/PLINK accept (order-stable)."""
    mapping: dict[str, int] = {}
    codes: list[int] = []
    for c in chrom:
        c = str(c)
        if c.isdigit():
            codes.append(int(c))
        else:
            codes.append(mapping.setdefault(c, 90 + len(mapping)))
    return codes


def decode_bed(fileset: PlinkFileset, n_samples: int) -> np.ndarray:
    """Decode a .bed back to an ALT-count matrix (n_variants, n_samples). For testing."""
    raw = fileset.bed.read_bytes()
    if raw[:3] != _MAGIC:
        raise ValueError("not a SNP-major PLINK .bed")
    body = np.frombuffer(raw[3:], dtype=np.uint8)
    bytes_per_var = (n_samples + 3) // 4
    body = body.reshape(-1, bytes_per_var)
    decode = {0b00: 0, 0b10: 1, 0b11: 2, 0b01: -1}
    out = np.empty((body.shape[0], bytes_per_var * 4), dtype=np.int8)
    for b in range(4):  # b-th 2-bit field within each byte -> samples b, b+4, b+8, ...
        codes = (body >> (2 * b)) & 0b11
        out[:, b::4] = np.vectorize(decode.get)(codes)
    return out[:, :n_samples]
