"""Golden regression for introgression statistics and the NJ tree."""

from __future__ import annotations

import json
import os
import sys
from io import StringIO
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.model import ArtifactKind
from canidae.pipeline import run_pipeline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_cohort import simulate_introgression_cohort

pytestmark = [pytest.mark.integration, pytest.mark.slow]

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden" / "introgression.json"
RESULT = ArtifactKind.ANALYSIS_RESULT


def _round(value: object, places: int = 6) -> object:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (int, bool, str)):
        return value
    return round(cast(float, value), places)


def _split_signature(newick: str) -> list[list[str]]:
    tree = Phylo.read(StringIO(newick), "newick")
    signatures = []
    for clade in tree.get_nonterminals():
        leaves = sorted(str(terminal.name) for terminal in clade.get_terminals())
        if 1 < len(leaves) < tree.count_terminals():
            signatures.append(leaves)
    return sorted(signatures)


@pytest.fixture(scope="module")
def observed(tmp_path_factory) -> dict:
    workspace = tmp_path_factory.mktemp("introgression-golden")
    vcf, sheet = simulate_introgression_cohort(
        workspace / "input", seed=11, admixture_proportion=0.2
    )
    cfg = GlobalConfig.load(
        overrides={
            "project_name": "introgression-golden",
            "paths.root": str(workspace),
            "pipeline": ["ingest", "load_genotypes", "nj_tree", "f3", "dstats"],
            "logging.level": "WARNING",
            "stages.ingest.sample_sheet": str(sheet),
            "stages.ingest.callset": str(vcf),
            "stages.nj_tree.n_bootstrap": 30,
            "stages.f3.outgroup": "jackal",
            "stages.dstats.outgroup": "jackal",
            "stages.dstats.block_mode": "site_count",
        }
    )
    report = run_pipeline(cfg)
    assert report.ok, report.failed
    store = DataStore(cfg.paths.data_root / "store")
    dstats = pd.read_csv(store.get(RESULT, "dstats").path)
    f3 = pd.read_csv(store.get(RESULT, "f3").path, index_col=0)
    top = dstats.iloc[dstats["Z"].abs().argmax()]
    return {
        "dstats": {
            "top_quartet": [top.P1, top.P2, top.P3, top.O],
            "top_D": _round(top.D, 5),
            "top_Z": _round(top.Z, 3),
            "top_q_value_bh": _round(top.q_value_bh, 6),
            "n_quartets": len(dstats),
        },
        "f3": {
            "outgroup": store.get(RESULT, "f3").metadata["outgroup"],
            "matrix": {
                row: {column: _round(value) for column, value in values.items()}
                for row, values in f3.iterrows()
            },
        },
        "nj": {"splits": _split_signature(store.get(ArtifactKind.TREE, "nj").path.read_text())},
    }


def test_introgression_estimators_match_golden(observed: dict) -> None:
    if os.environ.get("CANIDAE_UPDATE_GOLDEN"):
        GOLDEN.write_text(json.dumps(observed, indent=2) + "\n", encoding="utf-8")
        pytest.skip(f"golden snapshot rewritten: {GOLDEN}")
    assert GOLDEN.exists(), f"missing golden snapshot; regenerate {GOLDEN}"
    assert observed == json.loads(GOLDEN.read_text(encoding="utf-8"))
