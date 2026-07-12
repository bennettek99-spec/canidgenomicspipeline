"""Minimal VCF reading/writing helpers used by the harmonization stage.

Reading uses scikit-allel; writing is a small, dependency-free VCF 4.2 emitter sufficient
for biallelic-SNP genotype matrices (which is what the downstream analyses consume). Not a
general-purpose VCF library — just enough to round-trip harmonized callsets.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import allel
import numpy as np


@dataclass(slots=True)
class VcfSites:
    """A biallelic-SNP callset loaded into arrays."""

    chrom: np.ndarray
    pos: np.ndarray
    ref: np.ndarray
    alt: np.ndarray
    samples: np.ndarray
    gt: np.ndarray  # (n_sites, n_samples, 2), missing = -1

    @property
    def keys(self) -> np.ndarray:
        """Per-site identity key 'chrom:pos:ref:alt' for cross-dataset matching."""
        return np.array([f"{c}:{p}:{r}:{a}" for c, p, r, a in
                         zip(self.chrom, self.pos, self.ref, self.alt, strict=True)])


def read_biallelic_snps(path: Path) -> VcfSites:
    """Read a VCF and keep only biallelic SNP sites."""
    cs = allel.read_vcf(
        str(path),
        fields=["samples", "calldata/GT", "variants/CHROM", "variants/POS",
                "variants/REF", "variants/ALT"],
    )
    if cs is None or "calldata/GT" not in cs:
        raise ValueError(f"no genotypes in {path}")
    ref = np.asarray(cs["variants/REF"], dtype=str)
    alt = np.asarray(cs["variants/ALT"], dtype=str)
    n_alt = (alt != "").sum(axis=1)
    is_snp = (np.char.str_len(ref) == 1) & (np.char.str_len(alt[:, 0]) == 1)
    keep = (n_alt == 1) & is_snp
    return VcfSites(
        chrom=np.asarray(cs["variants/CHROM"], dtype=str)[keep],
        pos=np.asarray(cs["variants/POS"], dtype=np.int64)[keep],
        ref=ref[keep],
        alt=alt[keep, 0],
        samples=np.asarray(cs["samples"], dtype=str),
        gt=np.asarray(cs["calldata/GT"])[keep],
    )


def write_minimal_vcf(path: Path, chrom: np.ndarray, pos: np.ndarray, ref: np.ndarray,
                      alt: np.ndarray, samples: np.ndarray, gt: np.ndarray) -> Path:
    """Write a biallelic-SNP VCF 4.2 file from arrays. ``gt`` is (n_sites, n_samples, 2)."""
    path = Path(path)
    contigs = list(dict.fromkeys(str(c) for c in chrom))
    with path.open("w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for c in contigs:
            fh.write(f"##contig=<ID={c}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(map(str, samples)) + "\n")
        for i in range(len(pos)):
            calls = "\t".join(_gt_str(gt[i, s]) for s in range(len(samples)))
            fh.write(f"{chrom[i]}\t{int(pos[i])}\t.\t{ref[i]}\t{alt[i]}\t.\t.\t.\tGT\t"
                     f"{calls}\n")
    return path


def _gt_str(call: np.ndarray) -> str:
    a, b = int(call[0]), int(call[1])
    sa = "." if a < 0 else str(a)
    sb = "." if b < 0 else str(b)
    return f"{sa}/{sb}"
