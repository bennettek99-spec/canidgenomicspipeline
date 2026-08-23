"""Tests for preset citation bundles and how the report consumes them."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from canidae.core.citations import (
    CitationBundle,
    SourceCitation,
    citation_rows,
    load_citation_bundle,
    load_citation_bundles,
)
from canidae.core.config import GlobalConfig
from canidae.core.errors import ConfigError

ROOT = Path(__file__).resolve().parents[2]
BUNDLE_DIR = ROOT / "configs" / "citations"
SHIPPED = sorted(BUNDLE_DIR.glob("*.yaml"))


def test_repository_ships_citation_bundles() -> None:
    assert SHIPPED, "expected citation bundles under configs/citations/"


@pytest.mark.parametrize("path", SHIPPED, ids=lambda p: p.stem)
def test_every_shipped_bundle_validates(path: Path) -> None:
    bundle = load_citation_bundle(path)
    assert bundle.id == path.stem, "bundle id must match its filename"
    assert bundle.label
    assert bundle.sources, "a bundle with no sources cites nothing"
    keys = [source.key for source in bundle.sources]
    assert len(keys) == len(set(keys)), "source keys must be unique within a bundle"
    for source in bundle.sources:
        assert source.label
        assert source.used_for, f"{source.key} does not say what it is used for"
        # Every source must be resolvable by a reader: a URL, a DOI, or an
        # accession we know how to link.
        assert source.resolved_url, f"{source.key} has no resolvable link"


def test_bundles_resolve_by_bare_id_and_by_path() -> None:
    by_id = load_citation_bundles([Path("nhgri_722g_wgs")], ROOT)
    by_path = load_citation_bundles([Path("configs/citations/nhgri_722g_wgs.yaml")], ROOT)
    assert by_id[0].id == by_path[0].id == "nhgri_722g_wgs"


def test_duplicate_bundles_are_listed_once() -> None:
    bundles = load_citation_bundles([Path("nhgri_722g_wgs"), Path("nhgri_722g_wgs")], ROOT)
    assert len(bundles) == 1


def test_a_missing_bundle_is_a_config_error() -> None:
    with pytest.raises(ConfigError, match="not found"):
        load_citation_bundles([Path("no_such_bundle")], ROOT)


def test_an_invalid_bundle_is_a_config_error(tmp_path: Path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("id: bad\n", encoding="utf-8")  # missing required label
    with pytest.raises(ConfigError, match="invalid citation bundle"):
        load_citation_bundle(bad)

    not_a_mapping = tmp_path / "list.yaml"
    not_a_mapping.write_text("- one\n- two\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="must be a mapping"):
        load_citation_bundle(not_a_mapping)


def test_doi_prefixes_are_normalized() -> None:
    for raw in ("https://doi.org/10.1000/x", "doi:10.1000/x", "10.1000/x"):
        source = SourceCitation(key="k", label="l", doi=raw)
        assert source.doi == "10.1000/x"
        assert source.resolved_url == "https://doi.org/10.1000/x"


def test_resolved_url_falls_back_to_the_archive_link() -> None:
    source = SourceCitation(
        key="k", label="l", accession="PRJNA448733", accession_kind="bioproject"
    )
    assert source.resolved_url.endswith("/bioproject/PRJNA448733")
    assert SourceCitation(key="k", label="l").resolved_url == ""


def test_citation_rows_report_missing_identifiers_honestly() -> None:
    bundle = CitationBundle(
        id="b",
        label="Bundle",
        sources=[SourceCitation(key="k", label="Source", url="https://example.org")],
    )
    row = citation_rows([bundle])[0]
    assert row["dataset"] == "Bundle"
    assert row["accession"] == "not recorded"
    assert row["doi"] == "not recorded"


# -- preset wiring ---------------------------------------------------------------------


PRESETS_WITH_CITATIONS = [
    "nyc_coydog_validation.yaml",
    "eastern_coyote_ancestry.yaml",
    "redwolf_jackal_aadr.yaml",
    "redwolf_jackal_reduced_panel.yaml",
    "redwolf_laptop.yaml",
    "redwolf_real.yaml",
]


@pytest.mark.parametrize("preset", PRESETS_WITH_CITATIONS)
def test_preset_declares_resolvable_citation_bundles(preset: str) -> None:
    raw = yaml.safe_load((ROOT / "configs/examples" / preset).read_text(encoding="utf-8"))
    declared = raw["stages"]["report"]["citations"]
    assert declared, f"{preset} declares no citation bundles"
    bundles = load_citation_bundles([Path(entry) for entry in declared], ROOT)
    assert len(bundles) == len(set(declared))


def test_report_stage_accepts_the_citations_field(tmp_path: Path) -> None:
    from canidae.stages import load_builtin_stages
    from canidae.stages.reporting.report import ReportConfig

    load_builtin_stages()
    cfg = GlobalConfig.load(
        ROOT / "configs/examples/nyc_coydog_validation.yaml",
        overrides={"paths.root": str(tmp_path)},
    )
    report_cfg = cfg.parse_stage_config("report", ReportConfig)
    assert report_cfg.citations == ["nhgri_722g_wgs", "eastern_coyote_radseq"]
