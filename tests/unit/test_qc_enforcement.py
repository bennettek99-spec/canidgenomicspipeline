"""QC is an enforced input gate, not merely a flagging report."""

from __future__ import annotations

import numpy as np
import pandas as pd

from canidae.core.executor import NativeExecutor
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.stage import RunContext
from canidae.stages.popgen.load import LoadGenotypesStage
from canidae.stages.popgen.store import load_genotypes
from canidae.stages.processing.vcf_io import write_minimal_vcf
from canidae.stages.qc.callrate import SampleQCConfig, SampleQCStage


def test_qc_removes_failed_samples_and_sites_before_genotype_loading(
    tmp_context: RunContext,
) -> None:
    source_dir = tmp_context.datastore.stage_dir("source")
    vcf = source_dir / "input.vcf"
    samples = np.array(["A", "B", "C"])
    gt = np.array([
        [[0, 0], [-1, -1], [0, 1]],
        [[0, 1], [-1, -1], [-1, -1]],
        [[1, 1], [1, 1], [1, 1]],
    ], dtype=np.int8)
    write_minimal_vcf(
        vcf, np.array(["1", "1", "1"]), np.array([100, 200, 300]),
        np.array(["A", "A", "C"]), np.array(["G", "G", "T"]), samples, gt,
    )
    sheet = source_dir / "samples.csv"
    pd.DataFrame({
        "sample_id": samples,
        "taxon": ["gray_wolf", "gray_wolf", "coyote"],
        "population": ["wolf", "wolf", "coyote"],
    }).to_csv(sheet, index=False)
    tmp_context.datastore.register(
        Artifact(ArtifactKind.CALLSET, "callset", vcf, FileFormat.VCF)
    )
    tmp_context.datastore.register(
        Artifact(ArtifactKind.SAMPLE_SHEET, "sample_sheet", sheet, FileFormat.CSV)
    )

    qc = SampleQCStage(SampleQCConfig(min_sample_call_rate=0.6, min_site_call_rate=0.75))
    report = NativeExecutor().run([qc, LoadGenotypesStage()], tmp_context)

    assert report.ok, report.failed
    exclusions = pd.read_csv(tmp_context.datastore.get(ArtifactKind.QC_TABLE, "qc_exclusions").path)
    assert set(exclusions["entity"]) == {"B", "1:200:A:G"}
    assert "sample_call_rate_below_threshold" in set(exclusions["exclusion_reason"])
    assert "site_call_rate_below_threshold" in set(exclusions["exclusion_reason"])

    filtered = load_genotypes(tmp_context.datastore.get(ArtifactKind.GENOTYPES, "genotypes").path)
    assert filtered.samples.tolist() == ["A", "C"]
    assert filtered.pos.tolist() == [100, 300]
