from __future__ import annotations

from pathlib import Path

from canidae.core.model import ArtifactKind
from canidae.stages.sv.calling import (
    StructuralVariationStage,
    SvConfig,
    _resolve,
    sv_command,
)


def test_sv_command_chains_delly_and_bcftools() -> None:
    cmd = sv_command(["A.cram", "B.cram"], "ref.fasta", "out/sv.vcf")
    assert cmd.startswith("delly call -g ref.fasta -o sv.bcf A.cram B.cram")
    assert "bcftools view sv.bcf -Ov -o out/sv.vcf" in cmd
    assert "&&" in cmd  # bcftools conversion must not run if delly fails


def test_sv_command_single_cram() -> None:
    cmd = sv_command(["only.cram"], "r.fa", "sv.vcf")
    assert "-o sv.bcf only.cram" in cmd


def test_resolve_absolute_path_untouched(tmp_path: Path) -> None:
    absolute = tmp_path / "reference.fasta"
    assert _resolve(absolute, tmp_path) == absolute


def test_resolve_relative_path_anchored_at_root(tmp_path: Path) -> None:
    assert _resolve(Path("reference.fasta"), tmp_path) == tmp_path / "reference.fasta"


def test_stage_registration_and_io_contract() -> None:
    stage = StructuralVariationStage(SvConfig())
    assert stage.name == "sv"
    assert [(spec.kind, spec.role) for spec in stage.required_inputs()] == [
        (ArtifactKind.ALIGNMENT, "crams")
    ]
    assert [(spec.kind, spec.role) for spec in stage.produced_outputs()] == [
        (ArtifactKind.ANALYSIS_RESULT, "sv")
    ]


def test_config_defaults_match_documented_caller() -> None:
    cfg = SvConfig()
    assert cfg.sv_types == ["DEL", "DUP", "INV"]
    assert cfg.threads == 4
    assert cfg.reference == Path("reference.fasta")
