"""Admixture / ancestry-proportion estimation.

Estimates per-sample ancestry proportions (a Q matrix) over a range of K and selects the
best K by cross-validation. Two pluggable backends behind one stage — mirroring the
variant-calling design:

* ``nmf``    — in-Python sNMF-style weighted-NMF estimator (default; runs anywhere).
* ``binary`` — the ADMIXTURE executable via the tool runner (used where installed).
* ``auto``   — ``binary`` if the ADMIXTURE tool is on PATH, else ``nmf``.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from canidae.core.errors import ExternalToolError, StageInputError
from canidae.core.logging import get_logger
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.runtime import ResourceSpec, ToolSpec
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen import ancestry_nmf
from canidae.stages.popgen.plink import write_plink_bed
from canidae.stages.popgen.store import (
    Genotypes,
    align_labels,
    load_genotypes,
    load_sample_labels,
)

_log = get_logger("stage.admixture")
_ADMIXTURE = ToolSpec(name="admixture", version_args=("--version",))


class AdmixtureConfig(StageConfig):
    backend: str = "auto"           # auto | nmf | binary
    k_min: int = 2
    k_max: int = 6
    max_sites: int = 20000          # subsample sites above this (0 = use all)
    cv_holdout: float = 0.1
    n_iter: int = 250
    n_replicates: int = 3


@STAGES.register("admixture")
class AdmixtureStage(Stage):
    name = "admixture"
    config_model = AdmixtureConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "admixture")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: AdmixtureConfig = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        labels = align_labels(
            geno, load_sample_labels(
                ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path))

        k_values = self._k_range(cfg, geno.n_samples)
        backend = self._resolve_backend(cfg, ctx)
        _log.info("admixture backend=%s, K in %s", backend, k_values)

        if backend == "binary":
            best_k, Q, cv_errors, diagnostics = self._run_binary(cfg, ctx, geno, k_values)
        else:
            best_k, Q, cv_errors, diagnostics = self._run_nmf(cfg, ctx, geno, k_values)

        q_cols = [f"Q{i + 1}" for i in range(best_k)]
        table = pd.DataFrame(Q, columns=q_cols)
        table.insert(0, "population", labels["population"].to_numpy())
        table.insert(0, "sample_id", geno.samples)

        out = ctx.datastore.path_for(self.name, "admixture_Q.csv")
        table.to_csv(out, index=False)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "admixture", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "admixture",
                "backend": backend,
                "best_k": int(best_k),
                "k_values": k_values,
                "cv_errors": {int(k): round(float(v), 6) for k, v in cv_errors.items()},
                "q_columns": q_cols,
                "replicate_diagnostics": diagnostics,
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={"best_k": int(best_k), "backend": backend},
        )

    # -- backends ----------------------------------------------------------------------

    def _run_nmf(self, cfg, ctx, geno: Genotypes, k_values):
        X = self._allele_count_matrix(geno, cfg, seed=ctx.config.seed)
        replicate_cv = [
            ancestry_nmf.cross_validate_k(
                X, k_values, holdout=cfg.cv_holdout, n_iter=cfg.n_iter,
                seed=ctx.config.seed + replicate,
            )
            for replicate in range(cfg.n_replicates)
        ]
        cv_errors = {
            K: float(np.mean([values[K] for values in replicate_cv])) for K in k_values
        }
        best_k = ancestry_nmf.select_k(cv_errors)
        fits = [
            ancestry_nmf.fit_admixture(
                X, best_k, seed=ctx.config.seed + replicate, n_iter=cfg.n_iter
            )
            for replicate in range(cfg.n_replicates)
        ]
        fit = min(fits, key=lambda value: value.reconstruction_error)
        diagnostics = {
            "n_replicates": cfg.n_replicates,
            "reconstruction_errors": [round(value.reconstruction_error, 6) for value in fits],
            "cv_errors_by_replicate": [
                {int(k): round(float(v), 6) for k, v in values.items()}
                for values in replicate_cv
            ],
        }
        return best_k, fit.Q, cv_errors, diagnostics

    def _run_binary(self, cfg, ctx, geno: Genotypes, k_values):
        runner = ctx.runner
        runner.ensure(_ADMIXTURE)
        stage_dir = ctx.datastore.stage_dir(self.name)
        fileset = write_plink_bed(geno, stage_dir / "cohort")
        record = ctx.scratch.get("_record")
        cv_errors: dict[int, float] = {}
        q_by_k: dict[int, np.ndarray] = {}
        for K in k_values:
            result = runner.run(
                _ADMIXTURE, ["--cv", fileset.bed.name, str(K)],
                resources=ResourceSpec(cpus=ctx.config.resources.cpus),
                record=record, cwd=stage_dir,
                expect_outputs=[stage_dir / f"cohort.{K}.Q"],
            )
            cv_errors[K] = _parse_cv_error(result.stdout, K)
            q_by_k[K] = np.loadtxt(stage_dir / f"cohort.{K}.Q")
        best_k = min(cv_errors, key=lambda k: cv_errors[k])
        return best_k, np.atleast_2d(q_by_k[best_k]), cv_errors, {
            "n_replicates": 1,
            "backend": "admixture_binary",
        }

    # -- helpers -----------------------------------------------------------------------

    def _resolve_backend(self, cfg: AdmixtureConfig, ctx: RunContext) -> str:
        if cfg.backend == "nmf":
            return "nmf"
        if cfg.backend == "binary":
            return "binary"
        return "binary" if ctx.runner.which(_ADMIXTURE) else "nmf"

    def _k_range(self, cfg: AdmixtureConfig, n_samples: int) -> list[int]:
        hi = min(cfg.k_max, n_samples - 1)
        lo = max(2, cfg.k_min)
        if hi < lo:
            raise StageInputError(
                f"invalid K range for {n_samples} samples (k_min={cfg.k_min}, k_max={cfg.k_max})")
        return list(range(lo, hi + 1))

    def _allele_count_matrix(self, geno: Genotypes, cfg: AdmixtureConfig,
                             *, seed: int) -> np.ndarray:
        ac = geno.calls.count_alleles()
        seg = ac.is_segregating()
        X = geno.calls.to_n_alt().T[:, seg].astype(float)  # ALT counts, (n_samples, n_seg_sites)
        if X.shape[1] == 0:
            raise StageInputError("no segregating sites available for admixture")
        if cfg.max_sites and X.shape[1] > cfg.max_sites:
            rng = np.random.default_rng(seed)
            cols = rng.choice(X.shape[1], size=cfg.max_sites, replace=False)
            X = X[:, np.sort(cols)]
        return X


def _parse_cv_error(stdout: str, K: int) -> float:
    match = re.search(rf"CV error \(K={K}\):\s*([0-9.]+)", stdout)
    if not match:
        raise ExternalToolError(f"could not parse ADMIXTURE CV error for K={K}")
    return float(match.group(1))
