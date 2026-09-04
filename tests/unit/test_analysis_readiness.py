from __future__ import annotations

import json

import allel
import numpy as np
import pandas as pd

from canidae.core.executor import NativeExecutor
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.stage import RunContext
from canidae.stages.popgen.store import Genotypes, load_genotypes, save_genotypes
from canidae.stages.qc.readiness import (
    AnalysisReadinessConfig,
    AnalysisReadinessStage,
    _ld_prune,
)


def test_readiness_filters_non_autosomes_and_emits_small_population_warning(
    tmp_context: RunContext,
) -> None:
    source = tmp_context.datastore.stage_dir("source")
    gt = np.array(
        [
            [[0, 0], [0, 1], [1, 1], [0, 0]],
            [[0, 0], [0, 1], [1, 1], [0, 0]],
            [[0, 1], [0, 1], [0, 0], [1, 1]],
        ],
        dtype=np.int8,
    )
    store = save_genotypes(
        source / "genotypes.npz",
        Genotypes(
            allel.GenotypeArray(gt),
            np.array([100, 200, 300]),
            np.array(["1", "1", "X"]),
            np.array(["A", "B", "C", "D"]),
        ),
    )
    sheet = source / "samples.csv"
    pd.DataFrame(
        {
            "sample_id": ["A", "B", "C", "D"],
            "taxon": ["wolf"] * 4,
            "population": ["wolf", "wolf", "wolf", "jackal"],
        }
    ).to_csv(sheet, index=False)
    tmp_context.datastore.register(
        Artifact(
            ArtifactKind.GENOTYPES,
            "genotypes",
            store,
            FileFormat.NPZ,
            metadata={"panel_relative": True, "panel_name": "test_panel"},
        )
    )
    tmp_context.datastore.register(
        Artifact(ArtifactKind.SAMPLE_SHEET, "sample_sheet", sheet, FileFormat.CSV)
    )

    stage = AnalysisReadinessStage(AnalysisReadinessConfig(ld_prune=False))
    assert NativeExecutor().run([stage], tmp_context).ok
    ready = load_genotypes(
        tmp_context.datastore.get(ArtifactKind.GENOTYPES, "analysis_genotypes").path
    )
    assert set(ready.chrom) == {"1"}
    audit = json.loads(
        tmp_context.datastore.get(ArtifactKind.QC_TABLE, "analysis_readiness").path.read_text(
            encoding="utf-8"
        )
    )
    assert audit["status"] == "warning"
    assert {issue["code"] for issue in audit["issues"]} >= {
        "small_population",
        "ascertainment_panel_relative",
    }


def test_ld_prune_drops_perfectly_linked_duplicates_within_window() -> None:
    calls = np.zeros((4, 6, 2), dtype=np.int8)
    calls[0, :, 0] = [0, 0, 0, 1, 1, 1]
    calls[0, :, 1] = [0, 0, 0, 0, 1, 1]
    calls[1] = calls[0]  # perfect duplicate of variant 0
    calls[3, :, 0] = [1, 1, 1, 0, 0, 0]
    calls[3, :, 1] = [1, 1, 0, 0, 0, 0]
    geno = Genotypes(
        allel.GenotypeArray(calls),
        np.array([100, 150, 200, 300]),
        np.array(["1", "1", "1", "1"]),
        np.array([f"s{i}" for i in range(6)]),
    )
    assert list(_ld_prune(geno, np.arange(4), window_bp=500_000, r2_threshold=0.8)) == [0, 2, 3]
    assert list(_ld_prune(geno, np.arange(4), window_bp=10, r2_threshold=0.8)) == [0, 1, 2, 3]
