"""End-to-end demonstration of the CANIS foundation with two tiny, pure-Python stages.

This uses *no* external bioinformatics tools — it exists to show the full vertical slice
working: config -> stage registry -> DAG executor -> datastore artifacts -> provenance
manifest. Real stages (acquisition, qc, popgen, ...) follow exactly this pattern but shell
out to samtools/bcftools/PLINK/etc. through the tool runner.

Run it:  python examples/demo_pipeline.py
"""

from __future__ import annotations

import csv
import json
import tempfile
from collections import Counter
from pathlib import Path

from canidae.core.config import GlobalConfig
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.pipeline import run_pipeline

_DEMO_SAMPLES = [
    ("W001", "gray_wolf", "eurasia"),
    ("W002", "gray_wolf", "north_america"),
    ("C001", "coyote", "north_america"),
    ("C002", "coyote", "north_america"),
    ("D001", "domestic_dog", "europe"),
    ("V001", "village_dog", "africa"),
    ("DG1", "dingo", "australia"),
]


@STAGES.register("demo_ingest")
class DemoIngestStage(Stage):
    """Pretend to acquire a cohort: writes a sample sheet artifact."""

    name = "demo_ingest"
    config_model = StageConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.QC_TABLE, "sample_sheet")]

    def run(self, ctx: RunContext) -> StageResult:
        out = ctx.datastore.path_for(self.name, "samples.csv")
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_id", "taxon", "region"])
            writer.writerows(_DEMO_SAMPLES)
        art = Artifact(ArtifactKind.QC_TABLE, "sample_sheet", out, fmt=FileFormat.CSV)
        return StageResult(artifacts=[art], metrics={"n_samples": len(_DEMO_SAMPLES)})


@STAGES.register("demo_summarize")
class DemoSummarizeStage(Stage):
    """Pretend comparative analysis: counts individuals per taxon."""

    name = "demo_summarize"
    config_model = StageConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.QC_TABLE, "sample_sheet")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "taxon_counts")]

    def run(self, ctx: RunContext) -> StageResult:
        sheet = self.gather_inputs(ctx)[0]
        with sheet.path.open(encoding="utf-8") as fh:
            counts = Counter(row["taxon"] for row in csv.DictReader(fh))
        out = ctx.datastore.path_for(self.name, "taxon_counts.json")
        out.write_text(json.dumps(dict(counts), indent=2), encoding="utf-8")
        art = Artifact(ArtifactKind.ANALYSIS_RESULT, "taxon_counts", out,
                       fmt=FileFormat.JSON)
        return StageResult(artifacts=[art], metrics={"n_taxa": len(counts)})


def main() -> None:
    workspace = Path(tempfile.mkdtemp(prefix="canis-demo-"))
    cfg = GlobalConfig.load(
        overrides={
            "project_name": "canis-demo",
            "paths.root": str(workspace),
            "pipeline": ["demo_ingest", "demo_summarize"],
        }
    )

    report = run_pipeline(cfg)

    print("\n=== RESULT ===")
    print("executed:", report.executed)
    store_dir = cfg.paths.data_root / "store"
    counts_file = store_dir / "demo_summarize" / "taxon_counts.json"
    print("taxon counts:", counts_file.read_text(encoding="utf-8"))
    manifest = next((cfg.paths.run_root).glob("*/manifest.json"))
    print("manifest at:", manifest)


if __name__ == "__main__":
    main()
