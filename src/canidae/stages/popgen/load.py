"""Load a callset into a materialized genotype matrix for downstream analyses.

Reads the VCF once, restricts to analysis-appropriate sites (biallelic SNPs by default,
with optional MAF and per-site missingness filters), and writes a compact ``.npz`` the PCA,
F_ST, and diversity stages share. Keeping this a separate stage means the (potentially
expensive) parse happens once and the site-filtering policy is captured in one place and in
provenance.
"""

from __future__ import annotations

import allel
import numpy as np

from canidae.core.errors import ExternalToolError
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import Genotypes, save_genotypes


class LoadGenotypesConfig(StageConfig):
    biallelic_snps_only: bool = True
    min_maf: float = 0.0        # minor-allele-frequency floor (0 disables)
    max_missing: float = 1.0    # max per-site missing fraction (1 disables)


@STAGES.register("load_genotypes")
class LoadGenotypesStage(Stage):
    name = "load_genotypes"
    config_model = LoadGenotypesConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        # ``qc_callset`` is optional so old recipes that intentionally omit the QC stage
        # still work; when QC is present the DAG wires it ahead of this loader and the
        # enforced, filtered callset below is selected.
        return [
            ArtifactSpec(ArtifactKind.CALLSET, "callset"),
            ArtifactSpec(ArtifactKind.CALLSET, "qc_callset", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: LoadGenotypesConfig = self.config  # type: ignore[assignment]
        source = (
            ctx.datastore.get(ArtifactKind.CALLSET, "qc_callset")
            if ctx.datastore.has(ArtifactKind.CALLSET, "qc_callset")
            else ctx.datastore.get(ArtifactKind.CALLSET, "callset")
        )
        vcf = source.path

        callset = allel.read_vcf(
            str(vcf),
            fields=["samples", "calldata/GT", "variants/CHROM", "variants/POS",
                    "variants/REF", "variants/ALT"],
        )
        if callset is None or "calldata/GT" not in callset:
            raise ExternalToolError(f"no genotypes found in callset {vcf}")

        samples = np.asarray(callset["samples"], dtype=str)
        ga = allel.GenotypeArray(callset["calldata/GT"])
        pos = np.asarray(callset["variants/POS"], dtype=np.int64)
        chrom = np.asarray(callset["variants/CHROM"], dtype=str)
        n_input = ga.shape[0]

        keep = np.ones(n_input, dtype=bool)
        if cfg.biallelic_snps_only:
            keep &= _biallelic_snp_mask(callset)
        keep &= _site_filters(ga, samples.size, cfg)

        ga, pos, chrom = ga[keep], pos[keep], chrom[keep]
        if ga.shape[0] == 0:
            raise ExternalToolError(
                f"no sites survived filtering (from {n_input} input variants); "
                "loosen min_maf/max_missing or check the callset"
            )

        genotypes = Genotypes(calls=ga, pos=pos, chrom=chrom, samples=samples)
        out = save_genotypes(ctx.datastore.path_for(self.name, "genotypes.npz"), genotypes)

        art = Artifact(
            ArtifactKind.GENOTYPES, "genotypes", out, fmt=FileFormat.NPZ,
            produced_by=self.name,
            metadata={
                "n_samples": int(samples.size),
                "n_variants": int(ga.shape[0]),
                "source_callset_role": source.role,
                "qc_enforced": bool(source.metadata.get("qc_enforced", False)),
                "panel_relative": bool(source.metadata.get("panel_relative", False)),
                "panel_name": source.metadata.get("panel_name"),
                "callable_sites": source.metadata.get("callable_sites"),
                "reference_build": source.metadata.get("reference_build"),
            },
        )
        metrics = {
            "n_variants_input": int(n_input),
            "n_variants_kept": int(ga.shape[0]),
            "n_samples": int(samples.size),
        }
        return StageResult(artifacts=[art], metrics=metrics)


def _biallelic_snp_mask(callset: dict) -> np.ndarray:
    ref = np.asarray(callset["variants/REF"], dtype=str)
    alt = np.asarray(callset["variants/ALT"], dtype=str)
    n_alt = (alt != "").sum(axis=1)
    ref_is_base = np.char.str_len(ref) == 1
    alt0_is_base = np.char.str_len(alt[:, 0]) == 1
    return (n_alt == 1) & ref_is_base & alt0_is_base


def _site_filters(ga: allel.GenotypeArray, n_samples: int,
                  cfg: LoadGenotypesConfig) -> np.ndarray:
    mask = np.ones(ga.shape[0], dtype=bool)
    if cfg.max_missing < 1.0:
        call_rate = ga.is_called().sum(axis=1) / max(n_samples, 1)
        mask &= call_rate >= (1.0 - cfg.max_missing)
    if cfg.min_maf > 0.0:
        ac = ga.count_alleles()
        # Guard against monomorphic sites (allele total could be 0 after masking).
        totals = ac.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            alt_freq = np.where(totals > 0, ac[:, 1] / totals, 0.0)
        maf = np.minimum(alt_freq, 1.0 - alt_freq)
        mask &= maf >= cfg.min_maf
    return mask
