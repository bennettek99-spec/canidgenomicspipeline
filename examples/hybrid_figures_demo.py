"""Demonstrate the auto-generated hybrid-diagnostic figures.

Builds the deterministic synthetic bridge panel and runs all three hybrid
stages (reference_mixture -> breed_assign, plus multiway_admixture) together
with the report stage. The report embeds the figures as base64; this script
also writes the three figures as standalone PNGs so they can be inspected
directly.

Run it:  python examples/hybrid_figures_demo.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "tests" / "simulated"))
from make_hybrid_panel import build_hybrid_panel  # noqa: E402

from canidae.core.config import GlobalConfig  # noqa: E402
from canidae.core.datastore import DataStore  # noqa: E402
from canidae.core.model import ArtifactKind  # noqa: E402
from canidae.pipeline import run_pipeline  # noqa: E402
from canidae.stages import load_builtin_stages  # noqa: E402
from canidae.stages.reporting import figures  # noqa: E402

R = ArtifactKind.ANALYSIS_RESULT


def main() -> None:
    load_builtin_stages()
    workspace = Path(tempfile.mkdtemp(prefix="canis-hybrid-demo-"))
    panel = build_hybrid_panel(workspace / "panel")
    pedigree = [s for s in panel.two_source_query_ids() if s.startswith("NYFIXTURE")]

    cfg = GlobalConfig.load(overrides={
        "project_name": "hybrid-figures-demo",
        "paths.root": str(workspace),
        "pipeline": ["reference_mixture", "breed_assign", "multiway_admixture", "report"],
        "logging.level": "WARNING",
        "stages.reference_mixture.bridge_vcf": str(panel.bridge_vcf),
        "stages.reference_mixture.reference_genotypes": str(panel.reference_genotypes),
        "stages.reference_mixture.calls_dir": str(panel.calls_dir),
        "stages.reference_mixture.query_samples": pedigree,
        "stages.reference_mixture.wgs_samples": panel.reference_samples,
        "stages.reference_mixture.coyote_samples": panel.coyote_samples,
        "stages.breed_assign.bridge_vcf": str(panel.bridge_vcf),
        "stages.breed_assign.wgs_genotypes": str(panel.wgs_genotypes),
        "stages.breed_assign.calls_dir": str(panel.calls_dir),
        "stages.breed_assign.dog_fractions": str(panel.dog_fractions),
        "stages.breed_assign.locus_filter_json": str(panel.reference_genotypes),
        "stages.breed_assign.query_samples": pedigree,
        "stages.breed_assign.min_group": 3,
        "stages.multiway_admixture.bridge_vcf": str(panel.bridge_vcf),
        "stages.multiway_admixture.wgs_genotypes": str(panel.wgs_genotypes),
        "stages.multiway_admixture.western_ids": panel.query_ids(region="western"),
        "stages.multiway_admixture.bootstrap_n": 25,
        "stages.multiway_admixture.min_group": 3,
        "stages.report.title": "Hybrid figures demo (synthetic bridge panel)",
    })
    report = run_pipeline(cfg)
    print("executed:", report.executed)

    store = DataStore(cfg.paths.data_root / "store")
    report_html = store.get(ArtifactKind.REPORT, "html").path

    out_dir = _REPO / "runs" / "hybrid-figures-demo"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    if store.has(R, "reference_mixture"):
        figures.dog_fraction_barplot(
            store.get(R, "reference_mixture").path, fig_dir / "dog_fraction.png")
    if store.has(R, "breed_assign"):
        breed_art = store.get(R, "breed_assign")
        figures.breed_gap_barplot(
            breed_art.path, fig_dir / "breed_gap.png")
        scores = breed_art.metadata.get("scores_csv")
        if scores and Path(scores).exists():
            figures.breed_ranking_barplot(
                Path(scores), breed_art.path, fig_dir / "breed_ranking.png")
    if store.has(R, "multiway_admixture"):
        figures.admixture_scatter(
            store.get(R, "multiway_admixture").path, fig_dir / "admixture_scatter.png")

    print("report:", report_html)
    print("figures:")
    for png in sorted(fig_dir.glob("*.png")):
        print("  ", png)


if __name__ == "__main__":
    main()
