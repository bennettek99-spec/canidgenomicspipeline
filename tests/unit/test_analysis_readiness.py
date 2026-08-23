from __future__ import annotations

import json

import allel
import numpy as np
import pandas as pd

from canidae.core.executor import NativeExecutor
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.stage import RunContext
from canidae.stages.popgen.store import Genotypes, load_genotypes, save_genotypes
from canidae.stages.qc.readiness import AnalysisReadinessConfig, AnalysisReadinessStage


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
