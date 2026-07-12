from __future__ import annotations

import gzip
import hashlib
import struct
from pathlib import Path

import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.model import ArtifactKind
from canidae.pipeline import instantiate_stages
from canidae.stages.acquisition import reduced_panel as panel

ROOT = Path(__file__).resolve().parents[2]


def _tiny_panel() -> bytes:
    return gzip.compress(
        b"##fileformat=VCFv4.2\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        b"1\t100\t.\tA\tG\t.\t.\t.\n"
        b"1\t200\t.\tC\tT\t.\t.\t.\n"
    )


def _tiny_tbi() -> bytes:
    """A valid one-contig tabix skeleton; no chunks are needed for confirmation testing."""
    raw = bytearray(b"TBI\x01")
    raw.extend(struct.pack("<i", 1))
    raw.extend(struct.pack("<6i", 2, 1, 2, 0, ord("#"), 0))
    names = b"1\x00"
    raw.extend(struct.pack("<i", len(names)))
    raw.extend(names)
    raw.extend(struct.pack("<i", 0))  # n_bin
    raw.extend(struct.pack("<i", 0))  # n_linear
    return gzip.compress(bytes(raw))


def test_presets_and_checksum_validation() -> None:
    assert panel.PANEL_PRESETS == {"2k": 2_000, "10k": 10_000, "25k": 25_000}
    assert panel.preset_site_count("10k") == 10_000
    with pytest.raises(panel.ReducedPanelError, match="unknown reduced-panel preset"):
        panel.preset_site_count("1m")

    data = b"canis"
    digest = hashlib.sha256(data).hexdigest()
    assert panel.verify_checksum(data, f"sha256:{digest}", label="fixture") == f"sha256:{digest}"
    with pytest.raises(Exception, match="checksum mismatch"):
        panel.verify_checksum(data, "md5:00000000000000000000000000000000", label="fixture")


def test_transfer_policy_requires_confirmation_and_enforces_ceiling() -> None:
    estimate = panel.TransferEstimate(100, 500, 100, 3, 2_000)
    with pytest.raises(panel.TransferConfirmationRequiredError, match="confirm_large_transfer"):
        panel.enforce_transfer_policy(
            estimate, ceiling_bytes=1_000, confirmation_threshold_bytes=600, confirmed=False
        )
    panel.enforce_transfer_policy(
        estimate, ceiling_bytes=1_000, confirmation_threshold_bytes=600, confirmed=True
    )
    with pytest.raises(panel.TransferLimitExceededError, match="safety ceiling"):
        panel.enforce_transfer_policy(
            estimate, ceiling_bytes=699, confirmation_threshold_bytes=600, confirmed=True
        )


def test_choose_samples_preserves_explicit_order() -> None:
    sheet = pd.DataFrame(
        {
            "sample_id": ["wolf", "jackal", "coyote"],
            "taxon": ["red_wolf", "golden_jackal", "coyote"],
            "population": ["red_wolf", "jackal", "coyote"],
        }
    )
    chosen = panel.choose_samples(sheet, ["jackal", "wolf"])
    assert chosen["sample_id"].tolist() == ["jackal", "wolf"]
    with pytest.raises(panel.ReducedPanelError, match="missing"):
        panel.choose_samples(sheet, ["not-in-sheet"])


def test_confirmation_happens_before_source_vcf_ranges_and_temporary_files_are_cleaned(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[str, tuple[int, int] | None]] = []
    panel_data, index_data = _tiny_panel(), _tiny_tbi()

    def fake_request(url, budget, byte_range=None, **kwargs):
        calls.append((url, byte_range))
        data = panel_data if url == "https://example.test/panel.vcf.gz" else index_data
        budget.add(len(data))
        return data

    monkeypatch.setattr(panel, "_request", fake_request)
    with pytest.raises(panel.TransferConfirmationRequiredError):
        panel.extract_indexed_panel(
            source_vcf_url="https://example.test/source.vcf.gz",
            panel_url="https://example.test/panel.vcf.gz",
            panel_checksum="",
            index_checksum="",
            sample_ids=["A"],
            target_sites=2_000,
            preset="2k",
            output_path=tmp_path / "out.vcf.gz",
            manifest_path=tmp_path / "out.manifest.json",
            max_download_bytes=9_000_000_000,
            confirmation_threshold_bytes=0,
            confirm_large_transfer=False,
            temporary_dir=tmp_path,
        )

    # Only the small panel and tabix index were fetched; source-header/range traffic was blocked.
    assert calls == [
        ("https://example.test/panel.vcf.gz", None),
        ("https://example.test/source.vcf.gz.tbi", None),
    ]
    assert not list(tmp_path.glob(".reduced-panel-*"))
    assert not (tmp_path / "out.vcf.gz").exists()


def test_atomic_compact_vcf_contains_only_selected_samples(tmp_path: Path) -> None:
    output = tmp_path / "compact.vcf.gz"
    panel._atomic_gzip_vcf(
        output,
        ["##fileformat=VCFv4.2"],
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "A", "B"],
        ["B"],
        {
            ("2", 20): "2\t20\t.\tC\tT\t.\t.\t.\tGT\t0/1",
            ("1", 10): "1\t10\t.\tA\tG\t.\t.\t.\tGT\t1/1",
        },
        preset="2k",
    )
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert lines[-3].endswith("\tB")
    assert lines[-2].startswith("1\t10")
    assert lines[-1].startswith("2\t20")


def test_stage_wires_compact_callset_and_normalized_sheet(
    tmp_context, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "sample_id,taxon,population\nWolf25,red_wolf,red_wolf\nGoldenJackal01,golden_jackal,jackal\n",
        encoding="utf-8",
    )

    def fake_extract(**kwargs):
        kwargs["output_path"].write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n"))
        kwargs["manifest_path"].write_text("{}\n", encoding="utf-8")
        return panel.ReducedPanelResult(
            vcf_path=kwargs["output_path"],
            manifest_path=kwargs["manifest_path"],
            retained_sites=1_234,
            transfer_estimate=panel.TransferEstimate(1, 2, 3, 4, 10_000),
            downloaded_bytes=6,
            skipped_records={},
            output_sha256="abc",
        )

    monkeypatch.setattr(panel, "extract_indexed_panel", fake_extract)
    stage = panel.ReducedPanelStage(
        panel.ReducedPanelConfig(
            sample_sheet=sheet, selected_samples=["GoldenJackal01"], preset="10k"
        )
    )
    result = stage.run(tmp_context)
    assert {artifact.key for artifact in result.artifacts} == {
        (ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        (ArtifactKind.CALLSET, "callset"),
    }
    assert pd.read_csv(result.artifacts[0].path)["sample_id"].tolist() == ["GoldenJackal01"]
    callset = result.artifacts[1]
    assert callset.metadata["preset"] == "10k"
    assert callset.metadata["retained_sites"] == 1_234


def test_integrated_example_registers_reduced_panel_stage() -> None:
    cfg = GlobalConfig.load(ROOT / "configs/examples/redwolf_jackal_reduced_panel.yaml")
    stages = instantiate_stages(cfg)
    assert stages[0].name == "reduced_panel"
    assert [stage.name for stage in stages[1:4]] == ["qc", "load_genotypes", "distance"]
