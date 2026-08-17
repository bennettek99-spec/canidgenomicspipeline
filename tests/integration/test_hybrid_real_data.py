"""Regression guard on the real NYC coydog / eastern coyote results.

The study inputs live under ``data/``, which is gitignored, so these tests skip
everywhere the prepared panels are absent (CI included) and act as a local
promotion gate: they run the shipped presets over the real bridge panel and
compare against a committed snapshot, so a library change cannot quietly move
the published numbers.

To re-bless after an intended change, inspect the diff and then run::

    CANIDAE_UPDATE_GOLDEN=1 pytest tests/integration/test_hybrid_real_data.py
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.executor import build_context
from canidae.core.provenance import ProvenanceWriter
from canidae.core.registry import STAGES
from canidae.core.runtime import LocalRunner
from canidae.stages import load_builtin_stages

pytestmark = [pytest.mark.integration, pytest.mark.slow]

ROOT = Path(__file__).resolve().parents[2]
GOLDEN_DIR = ROOT / "tests" / "golden"
BRIDGE = ROOT / "data/eastern_coyote_wgs_bridge/eastern_coyote_wgs_bridge_1000.vcf.gz"
NYC_CALLS = ROOT / "data/nyc_coydog_validation/reference_genotypes.json"
BREED_PANEL = ROOT / "data/nyc_coydog_breeds/all_sample_genotypes.json"

_MISSING = [p for p in (BRIDGE, NYC_CALLS, BREED_PANEL) if not p.exists()]
pytestmark.append(
    pytest.mark.skipif(
        bool(_MISSING),
        reason=f"prepared study panels absent: {[p.name for p in _MISSING]}",
    )
)


def _run_preset_stage(preset: str, stage_name: str, tmp_path: Path):
    """Run one stage exactly as its shipped preset configures it."""
    load_builtin_stages()
    cfg = GlobalConfig.load(
        ROOT / "configs/examples" / preset,
        overrides={"paths.root": str(ROOT), "logging.level": "ERROR"},
    )
    stage_cls = STAGES.get(stage_name)
    stage = stage_cls(cfg.parse_stage_config(stage_name, stage_cls.config_model))
    run_dir = tmp_path / f"run-{stage_name}"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = build_context(
        cfg, DataStore(tmp_path / "store"), LocalRunner(),
        ProvenanceWriter(run_dir, config_digest=cfg.digest(), seed=cfg.seed),
        run_dir=run_dir,
    )
    return stage.run(ctx)


def _manifest(result) -> dict:
    return json.loads(Path(result.artifacts[0].metadata["manifest"]).read_text("utf-8"))


# The real panel takes minutes per stage, so each stage runs once per module and
# its result is shared by the tests that assert on different parts of it.


@pytest.fixture(scope="module")
def nyc_mixture(tmp_path_factory):
    return _run_preset_stage(
        "nyc_coydog_validation.yaml", "reference_mixture",
        tmp_path_factory.mktemp("nyc_mixture"),
    )


@pytest.fixture(scope="module")
def nyc_breeds(tmp_path_factory):
    return _run_preset_stage(
        "nyc_coydog_validation.yaml", "breed_assign",
        tmp_path_factory.mktemp("nyc_breeds"),
    )


@pytest.fixture(scope="module")
def eastern_multiway(tmp_path_factory):
    return _run_preset_stage(
        "eastern_coyote_ancestry.yaml", "multiway_admixture",
        tmp_path_factory.mktemp("eastern"),
    )


def _check(name: str, observed: dict) -> None:
    path = GOLDEN_DIR / name
    if os.environ.get("CANIDAE_UPDATE_GOLDEN"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(observed, indent=2) + "\n", encoding="utf-8")
        pytest.skip(f"golden snapshot rewritten: {path}")
    assert path.exists(), f"missing golden snapshot; regenerate {path}"
    assert observed == json.loads(path.read_text(encoding="utf-8"))


def _round(value: object, places: int = 6) -> object:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, (bool, int, str)):
        return value
    return round(float(value), places)


def test_nyc_dog_fractions_match_the_published_snapshot(nyc_mixture) -> None:
    result = nyc_mixture
    table = pd.read_csv(result.artifacts[0].path).set_index("sample_id")
    manifest = _manifest(result)
    _check("nyc_reference_mixture.json", {
        "dog_fraction": {
            sample: _round(table.loc[sample, "dog_fraction"])
            for sample in table.index
        },
        "diagnostic_sites_called": {
            sample: int(table.loc[sample, "diagnostic_sites_called"])
            for sample in table.index
        },
        "n_reference_loci": manifest["method"]["n_reference_loci"],
        "validation_passed": manifest["validation"]["passed"],
    })


def test_nyc_breed_assignment_matches_the_published_snapshot(nyc_breeds) -> None:
    manifest = _manifest(nyc_breeds)
    _check("nyc_breed_assign.json", {
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
    })


def test_no_nyc_sample_is_assigned_a_single_breed(nyc_breeds) -> None:
    """The honest-labelling invariant for a 252-locus bridge panel."""
    manifest = _manifest(nyc_breeds)
    assert manifest["results"], "expected at least one scored query"
    assert not any(
        info["single_breed_supported"] for info in manifest["results"].values()
    )


def test_eastern_coyote_summary_matches_the_published_snapshot(eastern_multiway) -> None:
    result = eastern_multiway
    table = pd.read_csv(result.artifacts[0].path).set_index("sample_id")
    manifest = _manifest(result)
    _check("eastern_multiway_admixture.json", {
        "queries": manifest["queries"],
        "references": manifest["references"],
        "group_summary": manifest["group_summary"],
        "per_sample": {
            sample: [_round(row["f_wolf"], 3), _round(row["f_dog"], 3)]
            for sample, row in table.sort_index().iterrows()
        },
    })


def test_western_controls_carry_no_wolf_or_dog_ancestry(eastern_multiway) -> None:
    """The negative control that makes the eastern signal interpretable."""
    table = pd.read_csv(eastern_multiway.artifacts[0].path)
    western = table[table["region"] == "western"]
    assert not western.empty
    assert western["f_wolf"].max() < 0.05
    assert western["f_dog"].max() < 0.05
