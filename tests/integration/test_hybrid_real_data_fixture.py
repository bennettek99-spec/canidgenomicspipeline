"""CI-runnable regression tests over the vendored public bridge-panel fixture.

The larger local ``data/`` panels remain the promotion gate. This small derived
fixture keeps the real-data regression surface available on every CI runner
without downloading public sources.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import cast

import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.executor import build_context
from canidae.core.provenance import ProvenanceWriter
from canidae.core.registry import STAGES
from canidae.core.runtime import LocalRunner
from canidae.stages import load_builtin_stages

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "fixtures" / "real_bridge"
GOLDEN = ROOT / "tests" / "golden"
pytestmark = [pytest.mark.integration]


def _manifest(result) -> dict:
    return json.loads(Path(result.artifacts[0].metadata["manifest"]).read_text("utf-8"))


def _round(value: object, places: int = 6) -> object:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, (bool, int, str)):
        return value
    return round(cast(float, value), places)


def _assert_golden(name: str, observed: dict) -> None:
    path = GOLDEN / name
    if os.environ.get("CANIDAE_UPDATE_GOLDEN"):
        path.write_text(json.dumps(observed, indent=2) + "\n", encoding="utf-8")
        pytest.skip(f"golden snapshot rewritten: {path}")
    assert path.exists(), f"missing golden snapshot; regenerate {path}"
    assert observed == json.loads(path.read_text(encoding="utf-8"))


def _run_stage(preset: str, stage_name: str, tmp_path: Path):
    load_builtin_stages()
    overrides: dict[str, object] = {
        "paths.root": str(tmp_path),
        "logging.level": "ERROR",
        f"stages.{stage_name}.bridge_vcf": str(FIXTURE / "bridge.vcf.gz"),
        f"stages.{stage_name}.wgs_genotypes": str(FIXTURE / "wgs_genotypes.json"),
    }
    if stage_name == "multiway_admixture":
        overrides["stages.multiway_admixture.bootstrap_n"] = 25
    cfg = GlobalConfig.load(
        ROOT / "configs/examples" / preset,
        overrides=overrides,
    )
    if stage_name in {"reference_mixture", "breed_assign"}:
        values = {
            "bridge_vcf": FIXTURE / "bridge.vcf.gz",
            "calls_dir": FIXTURE,
        }
        if stage_name == "reference_mixture":
            values.update(
                {
                    "reference_genotypes": FIXTURE / "reference_genotypes.json",
                }
            )
        else:
            values.update(
                {
                    "wgs_genotypes": FIXTURE / "wgs_genotypes.json",
                    "dog_fractions": FIXTURE / "validation.json",
                    "locus_filter_json": FIXTURE / "reference_genotypes.json",
                }
            )
        cfg = GlobalConfig.load(
            ROOT / "configs/examples" / preset,
            overrides={
                "paths.root": str(tmp_path),
                "logging.level": "ERROR",
                **{f"stages.{stage_name}.{key}": str(value) for key, value in values.items()},
            },
        )
    stage_cls = STAGES.get(stage_name)
    stage = stage_cls(cfg.parse_stage_config(stage_name, stage_cls.config_model))
    run_dir = tmp_path / f"run-{stage_name}"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = build_context(
        cfg,
        DataStore(tmp_path / "store"),
        LocalRunner(),
        ProvenanceWriter(run_dir, config_digest=cfg.digest(), seed=cfg.seed),
        run_dir=run_dir,
    )
    return stage.run(ctx)


@pytest.fixture(scope="module")
def real_fixture_runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("real-bridge-fixture")
    return {
        "mixture": _run_stage("nyc_coydog_validation.yaml", "reference_mixture", root / "mixture"),
        "breed": _run_stage("nyc_coydog_validation.yaml", "breed_assign", root / "breed"),
        "multiway": _run_stage(
            "eastern_coyote_ancestry.yaml", "multiway_admixture", root / "multiway"
        ),
    }


def test_fixture_reference_mixture_matches_real_golden(real_fixture_runs) -> None:
    result = real_fixture_runs["mixture"]
    table = pd.read_csv(result.artifacts[0].path).set_index("sample_id")
    manifest = _manifest(result)
    observed = {
        "dog_fraction": {
            sample: _round(table.loc[sample, "dog_fraction"]) for sample in table.index
        },
        "diagnostic_sites_called": {
            sample: int(table.loc[sample, "diagnostic_sites_called"]) for sample in table.index
        },
        "n_reference_loci": manifest["method"]["n_reference_loci"],
        "validation_passed": manifest["validation"]["passed"],
    }
    _assert_golden("nyc_reference_mixture_fixture.json", observed)


def test_fixture_breed_assignment_matches_real_golden(real_fixture_runs) -> None:
    manifest = _manifest(real_fixture_runs["breed"])
    observed = {
        "calibration": {
            "n_tested": manifest["calibration"]["n_tested"],
            "top1_correct": manifest["calibration"]["top1_correct"],
            "top1_accuracy": _round(manifest["calibration"]["top1_accuracy"]),
        },
        "sources": manifest["sources"],
        "results": {
            sample: {
                "best_breed": info["best_breed"],
                "single_breed_gap": _round(info["single_breed_gap"], 2),
                "single_breed_supported": info["single_breed_supported"],
                "loci_used": info["loci_used"],
            }
            for sample, info in sorted(manifest["results"].items())
        },
    }
    _assert_golden("nyc_breed_assign_fixture.json", observed)


def test_fixture_multiway_matches_real_golden(real_fixture_runs) -> None:
    result = real_fixture_runs["multiway"]
    table = pd.read_csv(result.artifacts[0].path).set_index("sample_id")
    manifest = _manifest(result)
    observed = {
        "queries": manifest["queries"],
        "references": manifest["references"],
        "group_summary": manifest["group_summary"],
        "per_sample": {
            sample: [_round(row["f_wolf"], 3), _round(row["f_dog"], 3)]
            for sample, row in table.sort_index().iterrows()
        },
    }
    _assert_golden("eastern_multiway_admixture_fixture.json", observed)
