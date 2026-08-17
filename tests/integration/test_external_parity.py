"""Parity checks against PLINK2 and Dsuite, when those binaries are installed.

CANIS reimplements the parts of these tools it needs so the pipeline stays
laptop-installable with no bioinformatics stack. That trade is only safe if the
reimplementations agree with the reference tools, so these tests run the real
binaries side by side whenever they are on PATH and skip otherwise. Nothing in
the default suite depends on them.

Run them explicitly with::

    pytest tests/integration/test_external_parity.py -v
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import allel
import numpy as np
import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.model import ArtifactKind
from canidae.pipeline import run_pipeline
from canidae.stages.popgen.plink import write_plink_bed
from canidae.stages.popgen.store import Genotypes, chunked_alt_frequency
from canidae.stages.processing.vcf_io import write_minimal_vcf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_cohort import simulate_introgression_cohort

pytestmark = [pytest.mark.integration, pytest.mark.slow]

PLINK2 = shutil.which("plink2")
DSUITE = shutil.which("Dsuite") or shutil.which("dsuite")
BCFTOOLS = shutil.which("bcftools")

requires_plink2 = pytest.mark.skipif(PLINK2 is None, reason="plink2 not on PATH")
requires_dsuite = pytest.mark.skipif(DSUITE is None, reason="Dsuite not on PATH")
requires_bcftools = pytest.mark.skipif(BCFTOOLS is None, reason="bcftools not on PATH")


def _demo_cohort(n_variants: int = 400, n_samples: int = 12, seed: int = 5) -> Genotypes:
    rng = np.random.default_rng(seed)
    freqs = rng.uniform(0.05, 0.95, n_variants)
    gt = np.stack(
        [rng.binomial(1, freqs, size=(n_samples, n_variants)).T for _ in range(2)],
        axis=-1,
    ).astype("i1")
    return Genotypes(
        calls=allel.GenotypeArray(gt),
        pos=np.arange(1, n_variants + 1, dtype=np.int64) * 1000,
        chrom=np.array(["1"] * n_variants),
        samples=np.array([f"S{i:02d}" for i in range(n_samples)]),
    )


# -- PLINK2 --------------------------------------------------------------------------


@requires_plink2
def test_plink2_reads_our_bed_and_agrees_on_allele_frequencies(tmp_path: Path) -> None:
    """Our .bed packing and our frequency math must both match PLINK2's reading."""
    geno = _demo_cohort()
    write_plink_bed(geno, tmp_path / "cohort")
    subprocess.run(
        [PLINK2, "--bfile", str(tmp_path / "cohort"), "--freq",
         "--allow-extra-chr", "--out", str(tmp_path / "freq")],
        check=True, capture_output=True, text=True,
    )
    report = pd.read_csv(tmp_path / "freq.afreq", sep="\t")
    assert len(report) == geno.n_variants

    ours = chunked_alt_frequency(geno)
    # Our .bim labels every site A/G with G as the ALT allele. Orient by the
    # letters PLINK2 echoes back rather than by column position, so the check
    # survives a change in how plink2 assigns REF/ALT to a PLINK1 fileset.
    theirs = report["ALT_FREQS"].to_numpy(dtype=float)
    alt_is_ours = report["ALT"].astype(str).to_numpy() == "G"
    expected = np.where(alt_is_ours, ours, 1.0 - ours)
    observed_sites = report["OBS_CT"].to_numpy(dtype=float) > 0
    assert np.allclose(theirs[observed_sites], expected[observed_sites], atol=1e-6)


@requires_plink2
def test_plink2_agrees_on_per_sample_missingness(tmp_path: Path) -> None:
    """The 0b01 missing code must survive the round trip into PLINK2."""
    geno = _demo_cohort(n_variants=200, n_samples=8)
    calls = np.asarray(geno.calls).copy()
    calls[:20, 0] = -1  # 20 missing sites for the first sample only
    geno = Genotypes(
        calls=allel.GenotypeArray(calls), pos=geno.pos,
        chrom=geno.chrom, samples=geno.samples,
    )
    write_plink_bed(geno, tmp_path / "cohort")
    subprocess.run(
        [PLINK2, "--bfile", str(tmp_path / "cohort"), "--missing", "sample-only",
         "--allow-extra-chr", "--out", str(tmp_path / "miss")],
        check=True, capture_output=True, text=True,
    )
    report = pd.read_csv(tmp_path / "miss.smiss", sep="\t").set_index("IID")
    assert int(report.loc[geno.samples[0], "MISSING_CT"]) == 20
    assert int(report.loc[geno.samples[1], "MISSING_CT"]) == 0


@requires_plink2
def test_plink2_pca_agrees_on_first_component(tmp_path: Path) -> None:
    """PLINK2 and scikit-allel must recover the same leading genotype axis."""
    geno = _demo_cohort(n_variants=500, n_samples=12)
    write_plink_bed(geno, tmp_path / "cohort")
    subprocess.run(
        [PLINK2, "--bfile", str(tmp_path / "cohort"), "--pca", "2",
         "--allow-extra-chr", "--out", str(tmp_path / "pca")],
        check=True, capture_output=True, text=True,
    )
    theirs = pd.read_csv(tmp_path / "pca.eigenvec", sep=r"\s+")
    ours, _model = allel.pca(
        np.asarray(geno.calls.to_n_alt(), dtype=np.float32),
        n_components=2,
        scaler="patterson",
    )
    order = {str(sample): index for index, sample in enumerate(geno.samples)}
    theirs_pc1 = np.array([float(value) for value in theirs["PC1"]])
    ours_pc1 = np.array([ours[order[str(sample)], 0] for sample in theirs["IID"]])
    assert abs(float(np.corrcoef(theirs_pc1, ours_pc1)[0, 1])) > 0.9


@requires_bcftools
def test_bcftools_stats_agrees_on_snp_count(tmp_path: Path) -> None:
    """The external callset summary must see the same biallelic SNP count as CANIS."""
    geno = _demo_cohort(n_variants=200, n_samples=8)
    vcf = write_minimal_vcf(
        tmp_path / "cohort.vcf", geno.chrom, geno.pos,
        np.full(geno.n_variants, "A"), np.full(geno.n_variants, "G"),
        geno.samples, np.asarray(geno.calls),
    )
    result = subprocess.run(
        [BCFTOOLS, "stats", str(vcf)], check=True, capture_output=True, text=True,
    )
    match = re.search(r"^SN\t[^\n]*number of SNPs:\t(\d+)", result.stdout, re.MULTILINE)
    assert match, result.stdout
    assert int(match.group(1)) == geno.n_variants


# -- Dsuite --------------------------------------------------------------------------


@requires_dsuite
def test_dsuite_agrees_on_the_simulated_introgression_signal(tmp_path: Path) -> None:
    """Both implementations must find the same COYOTE->DOG pulse in the same cohort."""
    vcf, sheet = simulate_introgression_cohort(
        tmp_path / "sim", seed=11, admixture_proportion=0.25
    )
    labels = pd.read_csv(sheet)

    # Dsuite wants sample<TAB>population, with the outgroup named "Outgroup".
    sets = tmp_path / "sets.txt"
    sets.write_text(
        "".join(
            f"{row.sample_id}\t"
            f"{'Outgroup' if row.population == 'JACKAL' else row.population}\n"
            for row in labels.itertuples(index=False)
        ),
        encoding="utf-8",
    )
    subprocess.run(
        [DSUITE, "Dtrios", "-o", str(tmp_path / "dsuite"), str(vcf), str(sets)],
        check=True, capture_output=True, text=True, cwd=tmp_path,
    )
    theirs = pd.read_csv(tmp_path / "dsuite_BBAA.txt", sep="\t")

    cfg = GlobalConfig.load(overrides={
        "project_name": "parity", "paths.root": str(tmp_path),
        "pipeline": ["ingest", "load_genotypes", "dstats"],
        "logging.level": "WARNING",
        "stages.ingest.sample_sheet": str(sheet), "stages.ingest.callset": str(vcf),
        "stages.dstats.outgroup": "JACKAL",
    })
    run_pipeline(cfg)
    store = DataStore(cfg.paths.data_root / "store")
    ours = pd.read_csv(store.get(ArtifactKind.ANALYSIS_RESULT, "dstats").path)

    def _key(p1: str, p2: str, p3: str) -> frozenset:
        return frozenset({p1, p2}), p3

    our_by_key = {_key(r.P1, r.P2, r.P3): abs(float(r.D)) for r in ours.itertuples()}
    compared = 0
    for row in theirs.itertuples(index=False):
        key = _key(row.P1, row.P2, row.P3)
        if key not in our_by_key:
            continue
        compared += 1
        # Sign depends on each tool's P1/P2 ordering, so compare magnitudes.
        assert our_by_key[key] == pytest.approx(abs(float(row.Dstatistic)), abs=0.03)
    assert compared >= 3, "no overlapping trios were compared"

    # And both must localize the signal to the same trio.
    our_top = ours.loc[ours["D"].abs().idxmax()]
    their_top = theirs.loc[theirs["Dstatistic"].abs().idxmax()]
    assert {our_top.P1, our_top.P2, our_top.P3} == {
        their_top.P1, their_top.P2, their_top.P3
    }
