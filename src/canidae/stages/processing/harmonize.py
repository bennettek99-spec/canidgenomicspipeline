"""Laptop-safe harmonization of pre-called VCFs on one reference build.

This stage deliberately does not launch a whole-genome liftover: cross-assembly conversion
is a preparation task that deserves independent review.  Instead it verifies the declared
build, normalizes biallelic SNP representations, reconciles safe REF/ALT swaps and strand
complements, and writes transparent batch/concordance diagnostics alongside the compact
merged VCF.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field

from canidae.core.errors import StageInputError
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.acquisition.sample_sheet import read_sample_sheet
from canidae.stages.processing.vcf_io import VcfSites, read_biallelic_snps, write_minimal_vcf


class DatasetInput(BaseModel):
    vcf: Path
    sample_sheet: Path
    dataset_id: str = "dataset"
    reference_build: str = ""  # takes precedence over a VCF ##reference/##assembly header


class HarmonizeConfig(StageConfig):
    datasets: list[DatasetInput] = Field(default_factory=list)
    reference_id: str = "unknown"
    cohort_id: str = "harmonized"
    require_reference_build: bool = False
    require_sample_concordance: bool = True
    normalize_variants: bool = True
    allow_strand_complement: bool = True
    liftover_requested: bool = False


@STAGES.register("harmonize")
class HarmonizeStage(Stage):
    name = "harmonize"
    config_model = HarmonizeConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []  # entry point: reads pre-called VCFs named in configuration

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.CALLSET, "callset"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
            ArtifactSpec(ArtifactKind.QC_TABLE, "harmonization_diagnostics"),
            ArtifactSpec(ArtifactKind.QC_TABLE, "harmonization_variant_orientation"),
        ]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: HarmonizeConfig = self.config  # type: ignore[assignment]
        if len(cfg.datasets) < 2:
            raise StageInputError("harmonize needs >= 2 datasets")
        if cfg.liftover_requested:
            raise StageInputError(
                "automatic liftover is deliberately disabled for laptop runs; prepare all "
                "VCFs on one reference build first, then rerun harmonize"
            )
        root = ctx.config.paths.root
        paths = [_resolve(ds.vcf, root) for ds in cfg.datasets]
        sheets = [read_sample_sheet(_resolve(ds.sample_sheet, root)) for ds in cfg.datasets]
        builds = [
            _reference_build(path, ds.reference_build)
            for path, ds in zip(paths, cfg.datasets, strict=True)
        ]
        reference_status = _verify_reference_builds(builds, cfg)

        sites = [read_biallelic_snps(path) for path in paths]
        if cfg.normalize_variants:
            sites = [_normalize_sites(s) for s in sites]
        sample_diagnostics = _sample_concordance(cfg.datasets, sites, sheets, builds)
        if cfg.require_sample_concordance:
            failures = sample_diagnostics.query("missing_metadata > 0 or sheet_only_samples > 0")
            if not failures.empty:
                bad = ", ".join(failures["dataset_id"].astype(str))
                raise StageInputError(
                    "VCF/sample-sheet IDs are not concordant for "
                    + bad
                    + "; see input metadata before merging"
                )
        _check_unique_samples(sites)

        shared, orientations, excluded = _shared_oriented_sites(
            sites, allow_complement=cfg.allow_strand_complement
        )
        if not shared:
            raise StageInputError(
                "datasets share no orientation-compatible biallelic SNPs on the same build"
            )
        merged_gt, samples = _merge_oriented(sites, shared, orientations)
        template = sites[0]
        idx0: np.ndarray = np.asarray(shared, dtype=int)
        stage_dir = ctx.datastore.stage_dir(self.name)
        out_vcf = stage_dir / "harmonized.vcf"
        _write_vcf_atomic(
            out_vcf,
            template.chrom[idx0],
            template.pos[idx0],
            template.ref[idx0],
            template.alt[idx0],
            samples,
            merged_gt,
        )

        sheet_df = _merge_sample_sheets(sheets)
        merged_sheet = _atomic_csv(sheet_df, stage_dir / "sample_sheet.csv")
        orientation_table = _orientation_report(cfg.datasets, orientations, excluded)
        orientation_out = _atomic_csv(orientation_table, stage_dir / "variant_orientation.csv")
        diagnostics = _batch_diagnostics(sample_diagnostics, sites, cfg.datasets, reference_status)
        diagnostics_out = _atomic_csv(diagnostics, stage_dir / "harmonization_diagnostics.csv")

        metadata = {
            "reference_id": cfg.reference_id,
            "reference_build_status": reference_status,
            "n_datasets": len(cfg.datasets),
            "n_shared_sites": len(shared),
            "n_samples": len(samples),
            "normalization": "biallelic SNP upper-case / canonical contigs"
            if cfg.normalize_variants
            else "disabled",
            "allow_strand_complement": cfg.allow_strand_complement,
        }
        callset_art = Artifact(
            ArtifactKind.CALLSET,
            "callset",
            out_vcf,
            FileFormat.VCF,
            produced_by=self.name,
            metadata=metadata,
        )
        sheet_art = Artifact(
            ArtifactKind.SAMPLE_SHEET,
            "sample_sheet",
            merged_sheet,
            FileFormat.CSV,
            produced_by=self.name,
            metadata={"n_samples": len(sheet_df)},
        )
        diagnostics_art = Artifact(
            ArtifactKind.QC_TABLE,
            "harmonization_diagnostics",
            diagnostics_out,
            FileFormat.CSV,
            produced_by=self.name,
            metadata={"reference_build_status": reference_status},
        )
        orientation_art = Artifact(
            ArtifactKind.QC_TABLE,
            "harmonization_variant_orientation",
            orientation_out,
            FileFormat.CSV,
            produced_by=self.name,
        )
        return StageResult(
            artifacts=[callset_art, sheet_art, diagnostics_art, orientation_art],
            metrics={
                "n_datasets": len(cfg.datasets),
                "n_shared_sites": len(shared),
                "n_samples": len(samples),
                "reference_build_status": reference_status,
                "excluded_incompatible_sites": int(excluded),
            },
        )


def _resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else (root / path)


def _reference_build(path: Path, explicit: str) -> str:
    if explicit.strip():
        return explicit.strip()
    opener = gzip.open if path.name.endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    break
                match = re.match(r"##(?:reference|assembly)=([^\s]+)", line.strip(), flags=re.I)
                if match:
                    return match.group(1)
    except OSError as exc:
        raise StageInputError(f"could not inspect VCF header {path}: {exc}") from exc
    return "unknown"


def _verify_reference_builds(builds: list[str], cfg: HarmonizeConfig) -> str:
    known = {build for build in builds if build and build.lower() != "unknown"}
    if cfg.require_reference_build and len(known) != len(builds):
        raise StageInputError(
            "reference build is missing for one or more VCFs; set dataset reference_build or "
            "add a ##reference header before harmonization"
        )
    if len(known) > 1:
        raise StageInputError(f"VCFs declare conflicting reference builds: {sorted(known)}")
    if cfg.reference_id.lower() != "unknown" and known and known != {cfg.reference_id}:
        raise StageInputError(
            f"VCF reference build {sorted(known)} does not match requested {cfg.reference_id}"
        )
    return "verified" if known and len(known) == 1 else "unverified"


def _normalize_sites(site: VcfSites) -> VcfSites:
    """Canonicalize the safe representation of biallelic SNPs without liftover."""
    return VcfSites(
        chrom=np.asarray([_canonical_chrom(c) for c in site.chrom], dtype=str),
        pos=site.pos.copy(),
        ref=np.char.upper(site.ref.astype(str)),
        alt=np.char.upper(site.alt.astype(str)),
        samples=site.samples.copy(),
        gt=site.gt.copy(),
    )


def _canonical_chrom(value: str) -> str:
    text = str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def _sample_concordance(
    datasets: list[DatasetInput],
    sites: list[VcfSites],
    sheets: list[pd.DataFrame],
    builds: list[str],
) -> pd.DataFrame:
    rows = []
    for ds, site, sheet, build in zip(datasets, sites, sheets, builds, strict=True):
        vcf_samples = set(map(str, site.samples))
        sheet_samples = set(sheet["sample_id"].astype(str))
        rows.append(
            {
                "dataset_id": ds.dataset_id,
                "reference_build": build,
                "n_vcf_samples": len(vcf_samples),
                "n_sheet_samples": len(sheet_samples),
                "matched_samples": len(vcf_samples & sheet_samples),
                "missing_metadata": len(vcf_samples - sheet_samples),
                "sheet_only_samples": len(sheet_samples - vcf_samples),
            }
        )
    return pd.DataFrame(rows)


def _check_unique_samples(sites: list[VcfSites]) -> None:
    seen: set[str] = set()
    for s in sites:
        dup = seen & set(map(str, s.samples))
        if dup:
            raise StageInputError(f"sample id(s) appear in multiple datasets: {sorted(dup)}")
        seen |= set(map(str, s.samples))


def _coord_key(site: VcfSites, index: int) -> str:
    return f"{site.chrom[index]}:{int(site.pos[index])}"


def _site_maps(site: VcfSites) -> dict[str, int]:
    mapping: dict[str, int] = {}
    for i in range(len(site.pos)):
        key = _coord_key(site, i)
        if key in mapping:
            raise StageInputError(f"duplicate biallelic SNP coordinate within a VCF: {key}")
        mapping[key] = i
    return mapping


def _orientation(
    template_ref: str, template_alt: str, ref: str, alt: str, allow_complement: bool
) -> tuple[str, bool] | None:
    if (ref, alt) == (template_ref, template_alt):
        return "same", False
    if (ref, alt) == (template_alt, template_ref):
        return "swap", True
    if allow_complement:
        # A/T and C/G pairs are strand-ambiguous without an independent reference/allele
        # frequency check. Keep direct matches above, but never guess a complemented match.
        if _is_palindromic(template_ref, template_alt):
            return None
        comp_ref, comp_alt = _complement(ref), _complement(alt)
        if (comp_ref, comp_alt) == (template_ref, template_alt):
            return "complement", False
        if (comp_ref, comp_alt) == (template_alt, template_ref):
            return "complement_swap", True
    return None


def _complement(base: str) -> str:
    return base.translate(str.maketrans("ACGT", "TGCA"))


def _is_palindromic(ref: str, alt: str) -> bool:
    return {ref, alt} in ({"A", "T"}, {"C", "G"})


def _shared_oriented_sites(
    sites: list[VcfSites], *, allow_complement: bool
) -> tuple[list[int], list[dict[int, tuple[int, str, bool]]], int]:
    """Return template indices and per-dataset orientation maps for compatible sites."""
    maps = [_site_maps(site) for site in sites]
    selected: list[int] = []
    orientations: list[dict[int, tuple[int, str, bool]]] = [dict() for _ in sites]
    incompatible = 0
    for i in range(len(sites[0].pos)):
        key = _coord_key(sites[0], i)
        matches: list[tuple[int, str, bool]] = [(i, "same", False)]
        valid = True
        for dataset_index, site in enumerate(sites[1:], start=1):
            other_i = maps[dataset_index].get(key)
            if other_i is None:
                valid = False
                break
            orientation = _orientation(
                str(sites[0].ref[i]),
                str(sites[0].alt[i]),
                str(site.ref[other_i]),
                str(site.alt[other_i]),
                allow_complement,
            )
            if orientation is None:
                valid = False
                incompatible += 1
                break
            mode, flip = orientation
            matches.append((other_i, mode, flip))
        if not valid:
            continue
        selected.append(i)
        for dataset_index, match in enumerate(matches):
            orientations[dataset_index][i] = match
    return selected, orientations, incompatible


def _merge_oriented(
    sites: list[VcfSites],
    template_indices: list[int],
    orientations: list[dict[int, tuple[int, str, bool]]],
) -> tuple[np.ndarray, np.ndarray]:
    matrices, samples = [], []
    for site, mapping in zip(sites, orientations, strict=True):
        selected = [mapping[index] for index in template_indices]
        indices: np.ndarray = np.asarray([entry[0] for entry in selected], dtype=int)
        flip: np.ndarray = np.asarray([entry[2] for entry in selected], dtype=bool)
        gt = site.gt[indices].copy()
        if flip.any():
            values = gt[flip]
            gt[flip] = np.where(values >= 0, 1 - values, -1)
        matrices.append(gt)
        samples.append(site.samples)
    return np.concatenate(matrices, axis=1), np.concatenate(samples)


def _orientation_report(
    datasets: list[DatasetInput],
    orientations: list[dict[int, tuple[int, str, bool]]],
    incompatible: int,
) -> pd.DataFrame:
    rows = []
    for dataset_index, (dataset, mapping) in enumerate(zip(datasets, orientations, strict=True)):
        counts = pd.Series([entry[1] for entry in mapping.values()]).value_counts().to_dict()
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "same_orientation_sites": int(counts.get("same", 0)),
                "ref_alt_swapped_sites": int(counts.get("swap", 0)),
                "strand_complemented_sites": int(counts.get("complement", 0)),
                "complemented_and_swapped_sites": int(counts.get("complement_swap", 0)),
                "incompatible_coordinate_alleles": incompatible if dataset_index else 0,
            }
        )
    return pd.DataFrame(rows)


def _batch_diagnostics(
    sample_diagnostics: pd.DataFrame,
    sites: list[VcfSites],
    datasets: list[DatasetInput],
    reference_status: str,
) -> pd.DataFrame:
    rows = []
    for row, site, _dataset in zip(
        sample_diagnostics.to_dict("records"), sites, datasets, strict=True
    ):
        called = np.all(site.gt >= 0, axis=2)
        sample_rates = called.mean(axis=0) if called.shape[1] else np.empty(0)
        site_rates = called.mean(axis=1) if called.shape[0] else np.empty(0)
        rows.append(
            {
                **row,
                "n_biallelic_snps": len(site.pos),
                "mean_sample_call_rate": round(float(sample_rates.mean()), 6)
                if sample_rates.size
                else None,
                "min_sample_call_rate": round(float(sample_rates.min()), 6)
                if sample_rates.size
                else None,
                "mean_site_call_rate": round(float(site_rates.mean()), 6)
                if site_rates.size
                else None,
                "reference_build_status": reference_status,
                "normalization_scope": "biallelic_snp_only",
                "liftover": "not_run",
            }
        )
    return pd.DataFrame(rows)


def _merge_sample_sheets(sheets: list[pd.DataFrame]) -> pd.DataFrame:
    merged = pd.concat(sheets, ignore_index=True)
    if merged["sample_id"].duplicated().any():
        raise StageInputError("duplicate sample_id across dataset sample sheets")
    return merged


def _atomic_csv(table: pd.DataFrame, destination: Path) -> Path:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    table.to_csv(temporary, index=False)
    temporary.replace(destination)
    return destination


def _write_vcf_atomic(
    path: Path,
    chrom: np.ndarray,
    pos: np.ndarray,
    ref: np.ndarray,
    alt: np.ndarray,
    samples: np.ndarray,
    gt: np.ndarray,
) -> Path:
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_minimal_vcf(temporary, chrom, pos, ref, alt, samples, gt)
    temporary.replace(path)
    return path
