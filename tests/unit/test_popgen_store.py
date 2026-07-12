from __future__ import annotations

from pathlib import Path

import allel
import numpy as np
import pandas as pd

from canidae.stages.popgen.store import (
    Genotypes,
    align_labels,
    load_genotypes,
    population_indices,
    save_genotypes,
)


def _demo_genotypes() -> Genotypes:
    gt = np.array(
        [
            [[0, 0], [0, 1], [1, 1]],
            [[0, 1], [1, 1], [0, 0]],
            [[0, 0], [0, 0], [0, 1]],
        ],
        dtype="i1",
    )
    return Genotypes(
        calls=allel.GenotypeArray(gt),
        pos=np.array([100, 200, 300], dtype=np.int64),
        chrom=np.array(["1", "1", "1"]),
        samples=np.array(["A", "B", "C"]),
    )


def test_save_load_roundtrip(tmp_path: Path) -> None:
    geno = _demo_genotypes()
    out = save_genotypes(tmp_path / "g.npz", geno)
    loaded = load_genotypes(out)
    assert loaded.n_variants == 3
    assert loaded.n_samples == 3
    assert list(loaded.samples) == ["A", "B", "C"]
    np.testing.assert_array_equal(np.asarray(loaded.calls), np.asarray(geno.calls))


def test_population_indices_aligns_to_matrix() -> None:
    geno = _demo_genotypes()
    labels = pd.DataFrame(
        {"taxon": ["gray_wolf", "gray_wolf", "coyote"],
         "population": ["wolf", "wolf", "coyote"]},
        index=["A", "B", "C"],
    )
    labels.index.name = "sample_id"
    groups = population_indices(geno, labels)
    assert groups == {"coyote": [2], "wolf": [0, 1]}


def test_align_labels_fills_unknown() -> None:
    geno = _demo_genotypes()
    labels = pd.DataFrame(
        {"taxon": ["gray_wolf"], "population": ["wolf"]}, index=["A"]
    )
    labels.index.name = "sample_id"
    aligned = align_labels(geno, labels)
    assert list(aligned["sample_id"]) == ["A", "B", "C"]
    assert aligned.loc[aligned["sample_id"] == "B", "population"].item() == "unknown"
