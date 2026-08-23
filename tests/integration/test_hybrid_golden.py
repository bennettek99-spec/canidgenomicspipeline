"""Golden snapshots of the hybrid stages on the deterministic synthetic panel.

These pin the numbers the three stages currently produce so a refactor cannot
silently move published-style results. They deliberately snapshot *point*
estimates only: bootstrap intervals depend on the NumPy random stream, which is
not a stable contract across releases, and their shape is already covered by
``test_hybrid_pipeline``.

To re-bless after an intended change, inspect the diff and then run::

    CANIDAE_UPDATE_GOLDEN=1 pytest tests/integration/test_hybrid_golden.py
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import cast

import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.executor import build_context
from canidae.core.provenance import ProvenanceWriter
from canidae.core.runtime import LocalRunner
from canidae.stages.hybrid.breed_assign import BreedAssignConfig, BreedAssignStage
from canidae.stages.hybrid.multiway_admixture import (
    MultiwayAdmixtureConfig,
    MultiwayAdmixtureStage,
)
from canidae.stages.hybrid.reference_mixture import (
    ReferenceMixtureConfig,
    ReferenceMixtureStage,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_hybrid_panel import build_hybrid_panel

pytestmark = [pytest.mark.integration]

GOLDEN = Path(__file__).resolve().parents[1] / "golden" / "hybrid_synthetic.json"
BOOTSTRAP_N = 25
SEED = 20260817


def _context(tmp_path: Path, pipeline: list[str]):
    cfg = GlobalConfig.load(
        include_defaults=True,
        overrides={
            "project_name": "hybrid-golden",
            "paths.root": str(tmp_path),
            "pipeline": pipeline,
            "logging.level": "ERROR",
        },
    )
    run_dir = tmp_path / f"run-{pipeline[0]}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return build_context(
        cfg,
        DataStore(tmp_path / "store"),
        LocalRunner(),
        ProvenanceWriter(run_dir, config_digest=cfg.digest(), seed=SEED),
        run_dir=run_dir,
    )


def _round(value: object, places: int = 4) -> object:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, (int, bool, str)):
        return value
    return round(cast(float, value), places)


def _observed(tmp_path: Path) -> dict:
    panel = build_hybrid_panel(tmp_path / "panel")
    pedigree = [s for s in panel.two_source_query_ids() if s.startswith("NYFIXTURE")]

    mixture = ReferenceMixtureStage(
        ReferenceMixtureConfig(
            bridge_vcf=panel.bridge_vcf,
            reference_genotypes=panel.reference_genotypes,
            calls_dir=panel.calls_dir,
            query_samples=pedigree,
            wgs_samples=panel.reference_samples,
            coyote_samples=panel.coyote_samples,
        )
    ).run(_context(tmp_path, ["reference_mixture"]))
    mixture_table = pd.read_csv(mixture.artifacts[0].path).set_index("sample_id")

    breed = BreedAssignStage(
        BreedAssignConfig(
            bridge_vcf=panel.bridge_vcf,
            wgs_genotypes=panel.wgs_genotypes,
            calls_dir=panel.calls_dir,
            dog_fractions=panel.dog_fractions,
            locus_filter_json=panel.reference_genotypes,
            query_samples=pedigree,
            min_group=3,
        )
    ).run(_context(tmp_path, ["breed_assign"]))
    breed_manifest = json.loads(
        Path(breed.artifacts[0].metadata["manifest"]).read_text(encoding="utf-8")
    )

    multiway = MultiwayAdmixtureStage(
        MultiwayAdmixtureConfig(
            bridge_vcf=panel.bridge_vcf,
            wgs_genotypes=panel.wgs_genotypes,
            western_ids=panel.query_ids(region="western"),
            bootstrap_n=BOOTSTRAP_N,
            seed=SEED,
            min_group=3,
        )
    ).run(_context(tmp_path, ["multiway_admixture"]))
    multiway_table = pd.read_csv(multiway.artifacts[0].path).set_index("sample_id")
    multiway_manifest = json.loads(
        Path(multiway.artifacts[0].metadata["manifest"]).read_text(encoding="utf-8")
    )

    return {
        "reference_mixture": {
            "dog_fraction": {
                sample: _round(mixture_table.loc[sample, "dog_fraction"]) for sample in pedigree
            },
            "diagnostic_sites_called": {
                sample: int(mixture_table.loc[sample, "diagnostic_sites_called"])
                for sample in pedigree
            },
        },
        "breed_assign": {
            "calibration": {
                "n_tested": breed_manifest["calibration"]["n_tested"],
                "top1_accuracy": _round(breed_manifest["calibration"]["top1_accuracy"]),
            },
            "results": {
                sample: {
                    "best_breed": info["best_breed"],
                    "single_breed_gap": _round(info["single_breed_gap"], 2),
                    "single_breed_supported": info["single_breed_supported"],
                    "loci_used": info["loci_used"],
                }
                for sample, info in sorted(breed_manifest["results"].items())
            },
        },
        "multiway_admixture": {
            "loci_informative": multiway_manifest["references"]["loci_informative"],
            "per_sample": {
                sample: [_round(row["f_wolf"], 3), _round(row["f_dog"], 3)]
                for sample, row in multiway_table.sort_index().iterrows()
            },
            "cohort": {
                name: [
                    _round(multiway_manifest["cohort_estimates"][name]["f_wolf"], 3),
                    _round(multiway_manifest["cohort_estimates"][name]["f_dog"], 3),
                ]
                for name in ("eastern", "western_control")
            },
        },
    }


def test_hybrid_stage_outputs_match_the_golden_snapshot(tmp_path: Path) -> None:
    observed = _observed(tmp_path)
    if os.environ.get("CANIDAE_UPDATE_GOLDEN"):
        GOLDEN.parent.mkdir(parents=True, exist_ok=True)
        GOLDEN.write_text(json.dumps(observed, indent=2) + "\n", encoding="utf-8")
        pytest.skip(f"golden snapshot rewritten: {GOLDEN}")

    assert GOLDEN.exists(), f"missing golden snapshot; regenerate {GOLDEN}"
    expected = json.loads(GOLDEN.read_text(encoding="utf-8"))
    assert observed == expected
