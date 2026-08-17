"""I/O helpers shared by stages and analysis scripts."""

from __future__ import annotations

from canidae.io.indexed_vcf import (
    BGZF_MAX_BLOCK,
    DEFAULT_LIMIT,
    SOURCE_VCF_URL,
    ReferenceIndex,
    TransferBudget,
    chunks_for,
    decompress_bgzf,
    download_small,
    genotype,
    merge_chunks,
    parse_tabix,
    reg2bins,
    request,
    source_header,
)

__all__ = [
    "BGZF_MAX_BLOCK",
    "DEFAULT_LIMIT",
    "SOURCE_VCF_URL",
    "ReferenceIndex",
    "TransferBudget",
    "chunks_for",
    "decompress_bgzf",
    "download_small",
    "genotype",
    "merge_chunks",
    "parse_tabix",
    "reg2bins",
    "request",
    "source_header",
]
