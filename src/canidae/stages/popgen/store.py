"""Materialized genotype store shared by the population-genomics stages.

Rather than re-parse the VCF in every analysis, one ``load_genotypes`` stage reads it once
into a compact on-disk matrix; PCA, F_ST, and diversity all read that matrix back. This is
the small-scale stand-in for the chunked Zarr/sgkit store used at the thousands-of-genomes
scale — same contract (load once, analyse many), different backing format.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import allel
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


@dataclass(slots=True)
class Genotypes:
    """An in-memory genotype matrix plus its coordinate and sample context."""

    calls: allel.GenotypeArray   # (n_variants, n_samples, ploidy)
    pos: np.ndarray              # (n_variants,) int
    chrom: np.ndarray            # (n_variants,) str
    samples: np.ndarray          # (n_samples,) str

    @property
    def n_variants(self) -> int:
        return self.calls.shape[0]

    @property
    def n_samples(self) -> int:
        return self.calls.shape[1]

    def allele_counts(self, subpop: list[int] | None = None) -> allel.AlleleCountsArray:
        return self.calls.count_alleles(subpop=subpop)

    def sample_index(self) -> dict[str, int]:
        return {s: i for i, s in enumerate(self.samples)}


def save_genotypes(path: Path, genotypes: Genotypes) -> Path:
    """Persist a :class:`Genotypes` to a compressed ``.npz``."""
    path = Path(path)
    np.savez_compressed(
        path,
        gt=np.asarray(genotypes.calls, dtype=np.int8),
        pos=np.asarray(genotypes.pos, dtype=np.int64),
        chrom=np.asarray(genotypes.chrom, dtype="U32"),
        samples=np.asarray(genotypes.samples, dtype="U64"),
    )
    # np.savez appends .npz if absent; normalize the returned path to the real file.
    return path if path.suffix == ".npz" else path.with_suffix(".npz")


def load_genotypes(path: Path) -> Genotypes:
    """Load a :class:`Genotypes` previously written by :func:`save_genotypes`."""
    with np.load(path, allow_pickle=False) as data:
        return Genotypes(
            calls=allel.GenotypeArray(data["gt"]),
            pos=data["pos"],
            chrom=data["chrom"].astype(str),
            samples=data["samples"].astype(str),
        )


def allele_difference_matrix(gn: np.ndarray) -> np.ndarray:
    """Pairwise allele-difference distance in [0, 1] from an alt-dosage matrix.

    ``gn`` is (n_sites, n_samples); returns an (n_samples, n_samples) matrix equal to the
    mean per-site absolute dosage difference divided by 2.
    """
    X = np.asarray(gn, dtype=float).T
    n_sites = max(X.shape[1], 1)
    return squareform(pdist(X, metric="cityblock") / (2.0 * n_sites))


def load_sample_labels(sample_sheet_csv: Path) -> pd.DataFrame:
    """Read the normalized sample sheet, indexed by ``sample_id``.

    Guaranteed columns: ``taxon``, ``population`` (plus any optional metadata columns).
    """
    df = pd.read_csv(sample_sheet_csv, dtype=str)
    return df.set_index("sample_id")


def align_labels(genotypes: Genotypes, labels: pd.DataFrame) -> pd.DataFrame:
    """Return per-sample labels aligned to the genotype-matrix sample order.

    Samples absent from the sample sheet are filled with ``"unknown"`` so downstream code
    never key-errors on an unlabeled sample.
    """
    rows = []
    for sample_id in genotypes.samples:
        if sample_id in labels.index:
            row = labels.loc[sample_id]
            rows.append({
                "sample_id": sample_id,
                "taxon": str(row.get("taxon", "unknown")),
                "population": str(row.get("population", "unknown")),
            })
        else:
            rows.append({"sample_id": sample_id, "taxon": "unknown",
                         "population": "unknown"})
    return pd.DataFrame(rows)


def population_indices(
    genotypes: Genotypes, labels: pd.DataFrame, *, column: str = "population"
) -> dict[str, list[int]]:
    """Map each population label to the genotype-matrix column indices of its samples.

    Only samples present in *both* the genotype matrix and the sample sheet are included;
    the mapping is keyed to the genotype sample order so allele counts align.
    """
    idx = genotypes.sample_index()
    groups: dict[str, list[int]] = {}
    for sample_id, row in labels.iterrows():
        col = idx.get(str(sample_id))
        if col is None:
            continue
        groups.setdefault(str(row[column]), []).append(col)
    return {k: sorted(v) for k, v in sorted(groups.items())}
