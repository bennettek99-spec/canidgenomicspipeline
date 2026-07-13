from pathlib import Path

from canidae.core.config import GlobalConfig
from canidae.pipeline import instantiate_stages
from canidae.stages.acquisition.sample_sheet import read_sample_sheet

ROOT = Path(__file__).resolve().parents[2]


def test_redwolf_jackal_aadr_config_and_samples_are_wired() -> None:
    cfg = GlobalConfig.load(ROOT / "configs/examples/redwolf_jackal_aadr.yaml")
    stages = instantiate_stages(cfg)
    sheet = read_sample_sheet(ROOT / "configs/examples/redwolf_jackal_aadr_samples.csv")

    assert [stage.name for stage in stages] == [
        "ingest",
        "qc",
        "load_genotypes",
        "analysis_readiness",
        "distance",
        "pca",
        "fst",
        "diversity",
        "nj_tree",
        "report",
    ]
    assert sheet["sample_id"].tolist() == ["Wolf25", "Wolf26", "GoldenJackal01"]
    assert sheet["taxon"].tolist() == ["red_wolf", "red_wolf", "golden_jackal"]
    assert cfg.stage_config("ingest")["callset"].endswith("caninehd25000.vcf.gz")


def test_redwolf_jackal_preparation_script_uses_indexed_panel_sources() -> None:
    script = (ROOT / "scripts/prepare_redwolf_jackal_aadr.sh").read_text(encoding="utf-8")

    assert "722g.990.SNP.INDEL.chrAll.vcf.gz" in script
    assert "CFA31_IlluminaHD.vcf.gz" in script
    assert 'curl --fail --silent --show-error --head "${SOURCE_VCF_URL}.tbi"' in script
    assert "Wolf25 Wolf26 GoldenJackal01" in script
    assert "bcftools view" in script
    assert "canidae run" not in script.split('echo ">> next command', maxsplit=1)[0]


def test_windows_preparation_script_has_a_sub_10gb_hard_limit() -> None:
    script = (ROOT / "scripts/prepare_redwolf_jackal_aadr.py").read_text(encoding="utf-8")

    assert "DEFAULT_LIMIT = 9_000_000_000" in script
    assert "status != 206" in script
    assert "server ignored byte-range request" in script
