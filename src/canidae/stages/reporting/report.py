# ruff: noqa: E501
"""Assemble a self-contained HTML report from the population-genomics results.

Reads the available result artifacts (PCA, admixture, clustering, F_ST, diversity, ROH,
geography, QC), renders the shared figures, and embeds everything — images as base64 data
URIs, tables as HTML — into one portable ``report.html`` via a Jinja2 template. Phase-2
sections render only when their result artifacts are present, so the same report stage
serves both a minimal and a full pipeline. Reporting only *reads* results.
"""

from __future__ import annotations

import base64
import html
import os
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
from jinja2 import Template

from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.reporting import figures

_TEMPLATE = Template(
    """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>{{ title }}</title>
<style>
  body { font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
         max-width: 1080px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }
  h1 { border-bottom: 3px solid #0072B2; padding-bottom: .3rem; }
  h2 { margin-top: 2.2rem; color: #0072B2; }
  .meta { color: #666; font-size: .85rem; }
  table { border-collapse: collapse; font-size: .85rem; margin: .5rem 0; }
  th, td { border: 1px solid #ddd; padding: 4px 8px; text-align: right; }
  th { background: #f4f4f4; }
  img { max-width: 100%; height: auto; }
  .note { background: #fff8e6; border-left: 4px solid #E69F00; padding: .6rem 1rem;
          font-size: .85rem; }
  .stat { font-size: .9rem; color: #333; }
  nav { position: sticky; top: 0; background: #fff; padding: .7rem 0; border-bottom: 1px solid #ddd;
        z-index: 2; display: flex; flex-wrap: wrap; gap: .8rem; font-size: .85rem; }
  nav a { color: #075985; text-decoration: none; }
  .badges { display: flex; flex-wrap: wrap; gap: .45rem; margin: .8rem 0; }
  .badge { border-radius: 999px; padding: .24rem .55rem; font-size: .8rem; font-weight: 600; }
  .badge.ok { background: #d1fae5; color: #065f46; }
  .badge.warn { background: #fef3c7; color: #92400e; }
  .badge.info { background: #dbeafe; color: #1e40af; }
  .package { border: 1px solid #d9e2ec; border-radius: .5rem; padding: .8rem 1rem;
             background: #f8fafc; }
  .package ul { columns: 2; padding-left: 1.2rem; }
  .interactive-pca { max-width: 680px; width: 100%; border: 1px solid #ddd; background: white; }
  .interactive-pca circle { cursor: pointer; }
  .interactive-pca circle:focus { stroke: #111; stroke-width: 3; }
  .small { font-size: .83rem; color: #4b5563; }
</style></head><body>
<nav aria-label="Report navigation">
  <a href="#summary">Summary</a><a href="#cohort">Cohort &amp; QC</a>
  <a href="#pca">PCA</a><a href="#diversity">Diversity</a><a href="#introgression">Introgression</a>
  <a href="#sources">Sources</a><a href="#limitations">Limitations</a><a href="#downloads">Outputs</a><a href="#methods">Methods</a>
</nav>
<h1>{{ title }}</h1>
<p class="meta">Project <b>{{ project }}</b> &middot; generated {{ generated }}
   &middot; config digest <code>{{ digest }}</code></p>
<div class="badges">{% for badge in badges %}<span class="badge {{ badge.kind }}" title="{{ badge.detail }}">{{ badge.label }}</span>{% endfor %}</div>
<p class="note">Comparative evolutionary genomics only &mdash; not for veterinary
   diagnostics or conservation-management decisions. Results are exploratory unless their
   callable-site definitions, cohort design, model assumptions, and external validation meet
   a study's publication standard.</p>

<h2 id="summary">Executive summary</h2>
<div class="package"><p>{{ executive_summary }}</p><p class="small">{{ run_status }}</p></div>

<h2 id="cohort">Included cohort and QC</h2>
{{ cohort_table }}
{% if excluded_table %}<h3>Excluded samples and sites</h3>{{ excluded_table }}{% endif %}

<h2 id="pca">Principal component analysis</h2>
<img src="data:image/png;base64,{{ pca_img }}" alt="PCA scatter">
{% if pca_svg %}<p class="small">Interactive PCA: hover or focus a point to inspect its sample and coordinates.</p>{{ pca_svg }}{% endif %}

{% if admixture_img %}<h2>Admixture</h2>
<p class="stat">Backend: {{ admix_backend }} &middot; selected K = {{ admix_k }}
   (cross-validation).</p>
<img src="data:image/png;base64,{{ admixture_img }}" alt="admixture barplot">{% endif %}

{% if dendro_img %}<h2>Population clustering</h2>
<p class="stat">Best k = {{ cluster_k }} &middot; silhouette = {{ cluster_sil }}
   &middot; adjusted Rand index vs population labels = {{ cluster_ari }}.</p>
<img src="data:image/png;base64,{{ dendro_img }}" alt="dendrogram">{% endif %}

{% if tree_img %}<h2>Phylogeny (neighbor-joining)</h2>
<img src="data:image/png;base64,{{ tree_img }}" alt="NJ tree">{% endif %}

<h2>Population differentiation (F<sub>ST</sub>)</h2>
<img src="data:image/png;base64,{{ fst_img }}" alt="FST heatmap">
{{ fst_table }}

<h2 id="diversity">Genetic diversity</h2>
{% if diversity_note %}<p class="note">{{ diversity_note }}</p>{% endif %}
<img src="data:image/png;base64,{{ div_img }}" alt="diversity barplot">
{{ div_table }}

{% if dstats_table or f3_table %}<h2 id="introgression">Introgression</h2>
{% if dstats_table %}<h3>D-statistics (ABBA-BABA)</h3>
<p class="stat">|Z| &gt; 3 indicates significant gene flow (P2&harr;P3). Outgroup:
   {{ dstats_outgroup }}.</p>
{{ dstats_table }}{% endif %}
{% if fd_img %}<img src="data:image/png;base64,{{ fd_img }}" alt="fd scan">{% endif %}
{% if f3_table %}<h3>Outgroup f<sub>3</sub> (shared drift)</h3>{{ f3_table }}{% endif %}{% endif %}

{% if karyo_img %}<h2>Local ancestry</h2>
<p class="stat">Sources: {{ lai_sources }}.</p>
<img src="data:image/png;base64,{{ karyo_img }}" alt="karyogram">
{{ lai_table }}{% endif %}

{% if selection_table %}<h2>Selection scan</h2>
<p class="stat">Focal population: {{ selection_focal }} (PBS vs {{ selection_refs }}).
   Top windows by PBS:</p>
{{ selection_table }}{% endif %}

{% if demography_table %}<h2>Demographic summary</h2>{{ demography_table }}{% endif %}

{% if roh_img %}<h2>Runs of homozygosity</h2>
<img src="data:image/png;base64,{{ roh_img }}" alt="ROH barplot">{% endif %}

{% if map_img or ibd_img %}<h2>Geography</h2>
{% if mantel_line %}<p class="stat">{{ mantel_line }}</p>{% endif %}
{% if map_img %}<img src="data:image/png;base64,{{ map_img }}" alt="locality map">{% endif %}
{% if ibd_img %}<img src="data:image/png;base64,{{ ibd_img }}" alt="IBD">{% endif %}
{% if regional_table %}<h3>Regional ancestry</h3>{{ regional_table }}{% endif %}{% endif %}

{% if qc_table %}<h2>Sample quality control</h2>{{ qc_table }}{% endif %}
{% if harmonization_table %}<h2>Callset harmonization diagnostics</h2>{{ harmonization_table }}{% endif %}

{% if sources_table %}<h2 id="sources">Source accessions and citations</h2>
<p class="small">Cite the source study/accession listed here and preserve its access conditions. The local manifest records the exact URLs and checksums used.</p>
{{ sources_table }}{% endif %}

<h2 id="limitations">Statistical limitations and warnings</h2>
<ul>{% for limitation in limitations %}<li>{{ limitation }}</li>{% endfor %}</ul>

<h2 id="downloads">Analysis package and downloads</h2>
<div class="package"><p>The self-contained HTML report can be opened by double-clicking this file. The companion results folder is the artifact store referenced below; no genomic files are duplicated for the report.</p>
<ul>{% for output in downloads %}<li><a href="{{ output.href }}">{{ output.label }}</a> <span class="small">{{ output.kind }}</span></li>{% endfor %}</ul></div>

<h2 id="methods">Methods, provenance, and reproducibility</h2>
<div class="package"><p>{{ methods }}</p><p class="small">Run manifest: <code>{{ manifest_path }}</code></p></div>
</body></html>"""
)


class ReportConfig(StageConfig):
    title: str = "Canid population-genomics report"
    exploratory: bool = True


@STAGES.register("report")
class ReportStage(Stage):
    name = "report"
    config_model = ReportConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        R = ArtifactKind.ANALYSIS_RESULT
        return [
            ArtifactSpec(R, "pca"),
            ArtifactSpec(R, "fst"),
            ArtifactSpec(R, "diversity"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
            ArtifactSpec(R, "admixture", optional=True),
            ArtifactSpec(R, "cluster", optional=True),
            ArtifactSpec(R, "distance", optional=True),
            ArtifactSpec(R, "roh", optional=True),
            ArtifactSpec(R, "geography", optional=True),
            ArtifactSpec(R, "dstats", optional=True),
            ArtifactSpec(R, "f3", optional=True),
            ArtifactSpec(ArtifactKind.TREE, "nj", optional=True),
            ArtifactSpec(R, "local_ancestry", optional=True),
            ArtifactSpec(R, "selection", optional=True),
            ArtifactSpec(R, "demography", optional=True),
            ArtifactSpec(ArtifactKind.QC_TABLE, "sample_qc", optional=True),
            ArtifactSpec(ArtifactKind.QC_TABLE, "qc_exclusions", optional=True),
            ArtifactSpec(ArtifactKind.QC_TABLE, "harmonization_diagnostics", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "qc_sample_sheet", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.REPORT, "html")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: ReportConfig = self.config  # type: ignore[assignment]
        ds = ctx.datastore
        R = ArtifactKind.ANALYSIS_RESULT
        d = ds.stage_dir(self.name)

        pca_art = ds.get(R, "pca")
        fst_art = ds.get(R, "fst")
        div_art = ds.get(R, "diversity")
        sheet = (
            ds.get(ArtifactKind.SAMPLE_SHEET, "qc_sample_sheet").path
            if ds.has(ArtifactKind.SAMPLE_SHEET, "qc_sample_sheet")
            else ds.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path
        )
        diversity_meta = div_art.metadata
        limitations = _limitations(ds, diversity_meta)

        ctx_vars: dict[str, object] = {
            "title": cfg.title,
            "project": ctx.config.project_name,
            "generated": datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"),
            "digest": ctx.config.digest()[:16],
            "cohort_table": _cohort_table(sheet),
            "pca_img": _b64(figures.pca_scatter(
                pca_art.path, d / "pca.png",
                explained_variance=pca_art.metadata.get("explained_variance_ratio"))),
            "fst_img": _b64(figures.fst_heatmap(fst_art.path, d / "fst.png")),
            "fst_table": pd.read_csv(fst_art.path, index_col=0).to_html(border=0),
            "div_img": _b64(figures.diversity_bar(div_art.path, d / "diversity.png")),
            "div_table": pd.read_csv(div_art.path).to_html(index=False, border=0),
            "diversity_note": diversity_meta.get("limitations"),
            "badges": _badges(ds, diversity_meta, cfg.exploratory),
            "limitations": limitations,
            "executive_summary": _executive_summary(ds, diversity_meta),
            "run_status": _run_status(ctx),
            # `d` is a transactional staging directory; links must be computed from the
            # final promoted report directory so they keep working after double-clicking.
            "downloads": _downloads(ds, ctx.datastore.root / self.name),
            "manifest_path": str(ctx.run_dir / "manifest.json") if ctx.run_dir else "unavailable",
            "methods": _methods(ctx, diversity_meta),
            "pca_svg": _pca_svg(pca_art.path),
            "sources_table": _sources_table(ds),
        }

        self._add_admixture(ds, d, ctx_vars)
        self._add_clustering(ds, d, ctx_vars)
        self._add_phylogeny(ds, d, ctx_vars, sheet)
        self._add_introgression(ds, d, ctx_vars)
        self._add_local_ancestry(ds, d, ctx_vars)
        self._add_selection(ds, ctx_vars)
        self._add_demography(ds, ctx_vars)
        self._add_roh(ds, d, ctx_vars)
        self._add_geography(ds, d, ctx_vars)
        self._add_qc(ds, ctx_vars)
        self._add_harmonization(ds, ctx_vars)

        html = _TEMPLATE.render(**ctx_vars)
        out = d / "report.html"
        temporary = out.with_suffix(".html.tmp")
        temporary.write_text(html, encoding="utf-8")
        temporary.replace(out)
        art = Artifact(ArtifactKind.REPORT, "html", out, fmt=FileFormat.HTML,
                       produced_by=self.name)
        return StageResult(artifacts=[art], metrics={"report": str(out)})

    # -- optional sections -------------------------------------------------------------

    def _add_admixture(self, ds, d: Path, v: dict) -> None:
        if not ds.has(ArtifactKind.ANALYSIS_RESULT, "admixture"):
            return
        art = ds.get(ArtifactKind.ANALYSIS_RESULT, "admixture")
        v["admixture_img"] = _b64(figures.admixture_barplot(art.path, d / "admixture.png"))
        v["admix_backend"] = art.metadata.get("backend", "?")
        v["admix_k"] = art.metadata.get("best_k", "?")

    def _add_clustering(self, ds, d: Path, v: dict) -> None:
        if not (ds.has(ArtifactKind.ANALYSIS_RESULT, "cluster")
                and ds.has(ArtifactKind.ANALYSIS_RESULT, "distance")):
            return
        dist = ds.get(ArtifactKind.ANALYSIS_RESULT, "distance")
        meta = ds.get(ArtifactKind.ANALYSIS_RESULT, "cluster").metadata
        v["dendro_img"] = _b64(figures.cluster_dendrogram(dist.path, d / "dendro.png"))
        v["cluster_k"] = meta.get("best_k", "?")
        v["cluster_sil"] = meta.get("silhouette", "?")
        v["cluster_ari"] = meta.get("adjusted_rand_index", "?")

    def _add_phylogeny(self, ds, d: Path, v: dict, sheet: Path) -> None:
        if not ds.has(ArtifactKind.TREE, "nj"):
            return
        newick = ds.get(ArtifactKind.TREE, "nj").path
        v["tree_img"] = _b64(figures.nj_tree_figure(newick, d / "nj_tree.png",
                                                     tip_colors=_tip_colors(sheet)))

    def _add_introgression(self, ds, d: Path, v: dict) -> None:
        R = ArtifactKind.ANALYSIS_RESULT
        if ds.has(R, "dstats"):
            art = ds.get(R, "dstats")
            df = pd.read_csv(art.path)
            if not df.empty:
                v["dstats_table"] = df.to_html(index=False, border=0)
                v["dstats_outgroup"] = art.metadata.get("outgroup", "?")
            fd = art.path.parent / "fd_windows.csv"
            if fd.exists():
                v["fd_img"] = _b64(figures.fd_scan_plot(fd, d / "fd_scan.png"))
        if ds.has(R, "f3"):
            v["f3_table"] = pd.read_csv(ds.get(R, "f3").path, index_col=0).to_html(border=0)

    def _add_local_ancestry(self, ds, d: Path, v: dict) -> None:
        R = ArtifactKind.ANALYSIS_RESULT
        if not ds.has(R, "local_ancestry"):
            return
        art = ds.get(R, "local_ancestry")
        windows = art.path.parent / "local_ancestry_windows.csv"
        if windows.exists() and len(pd.read_csv(windows)):
            v["karyo_img"] = _b64(figures.karyogram(windows, d / "karyogram.png"))
        v["lai_sources"] = ", ".join(art.metadata.get("sources", []))
        v["lai_table"] = pd.read_csv(art.path).to_html(index=False, border=0)

    def _add_selection(self, ds, v: dict) -> None:
        R = ArtifactKind.ANALYSIS_RESULT
        if not ds.has(R, "selection"):
            return
        art = ds.get(R, "selection")
        df = pd.read_csv(art.path).sort_values("pbs", ascending=False).head(10)
        v["selection_table"] = df.to_html(index=False, border=0)
        v["selection_focal"] = art.metadata.get("focal", "?")
        v["selection_refs"] = ", ".join(art.metadata.get("references", []))

    def _add_demography(self, ds, v: dict) -> None:
        R = ArtifactKind.ANALYSIS_RESULT
        if ds.has(R, "demography"):
            v["demography_table"] = pd.read_csv(
                ds.get(R, "demography").path).to_html(index=False, border=0)

    def _add_roh(self, ds, d: Path, v: dict) -> None:
        if ds.has(ArtifactKind.ANALYSIS_RESULT, "roh"):
            art = ds.get(ArtifactKind.ANALYSIS_RESULT, "roh")
            v["roh_img"] = _b64(figures.roh_barplot(art.path, d / "roh.png"))

    def _add_geography(self, ds, d: Path, v: dict) -> None:
        if not ds.has(ArtifactKind.ANALYSIS_RESULT, "geography"):
            return
        art = ds.get(ArtifactKind.ANALYSIS_RESULT, "geography")
        geo_dir = art.path.parent
        localities = geo_dir / "localities.csv"
        ibd = geo_dir / "ibd_pairs.csv"
        regional = geo_dir / "regional_ancestry.csv"
        mantel_r = art.metadata.get("mantel_r")
        if localities.exists() and pd.read_csv(localities)["latitude"].notna().any():
            v["map_img"] = _b64(figures.locality_map(localities, d / "map.png"))
        if ibd.exists():
            v["ibd_img"] = _b64(figures.ibd_scatter(
                ibd, d / "ibd.png",
                mantel_r=float(mantel_r) if mantel_r is not None else None))
        if mantel_r is not None:
            v["mantel_line"] = (f"Isolation by distance: Mantel r = {mantel_r}, "
                                f"p = {art.metadata.get('mantel_p')} "
                                f"({art.metadata.get('n_localities')} localities).")
        if regional.exists():
            v["regional_table"] = pd.read_csv(regional).to_html(index=False, border=0)

    def _add_qc(self, ds, v: dict) -> None:
        if ds.has(ArtifactKind.QC_TABLE, "sample_qc"):
            qc = pd.read_csv(ds.get(ArtifactKind.QC_TABLE, "sample_qc").path)
            v["qc_table"] = qc.to_html(index=False, border=0)
        if ds.has(ArtifactKind.QC_TABLE, "qc_exclusions"):
            excluded = pd.read_csv(ds.get(ArtifactKind.QC_TABLE, "qc_exclusions").path)
            v["excluded_table"] = excluded.to_html(index=False, border=0)

    def _add_harmonization(self, ds, v: dict) -> None:
        if ds.has(ArtifactKind.QC_TABLE, "harmonization_diagnostics"):
            table = pd.read_csv(ds.get(ArtifactKind.QC_TABLE, "harmonization_diagnostics").path)
            v["harmonization_table"] = table.to_html(index=False, border=0)


def _cohort_table(sheet: Path) -> str:
    df = pd.read_csv(sheet)
    return (df.groupby(["taxon", "population"]).size().reset_index(name="n_samples")
            .to_html(index=False, border=0))


def _tip_colors(sheet: Path) -> dict[str, str]:
    """Map each sample_id to a color by its population (for tree tip labels)."""
    df = pd.read_csv(sheet, dtype=str)
    pops = sorted(df["population"].dropna().unique())
    palette = {p: figures._PALETTE[i % len(figures._PALETTE)] for i, p in enumerate(pops)}
    return {row.sample_id: palette.get(row.population, "#333333")
            for row in df.itertuples(index=False)}


def _b64(png_path: Path) -> str:
    return base64.b64encode(png_path.read_bytes()).decode("ascii")


def _badges(ds, diversity_metadata: dict, exploratory: bool) -> list[dict[str, str]]:
    badges = [{
        "kind": "warn" if exploratory else "ok",
        "label": "Exploratory" if exploratory else "Publication workflow configured",
        "detail": "Interpret alongside study design, callable-site definitions, and validation.",
    }]
    if diversity_metadata.get("scope") == "panel_relative":
        badges.append({
            "kind": "warn", "label": "Panel-relative diversity",
            "detail": "Selected SNPs are not a whole-genome callable-site denominator.",
        })
    else:
        badges.append({
            "kind": "ok", "label": "Callable-site denominator supplied",
            "detail": "Diversity output uses an explicit callable-site denominator.",
        })
    if ds.has(ArtifactKind.QC_TABLE, "qc_exclusions"):
        exclusions = pd.read_csv(ds.get(ArtifactKind.QC_TABLE, "qc_exclusions").path)
        if len(exclusions):
            badges.append({
                "kind": "warn", "label": f"QC excluded {len(exclusions)} records",
                "detail": "See the included exclusion table for every reason.",
            })
        else:
            badges.append({
                "kind": "ok", "label": "QC gate passed", "detail": "No samples or sites were excluded."
            })
    else:
        badges.append({
            "kind": "info", "label": "QC not in this recipe",
            "detail": "This report was produced without the enforced QC stage.",
        })
    return badges


def _limitations(ds, diversity_metadata: dict) -> list[str]:
    values = [
        "Exploratory analysis: this package is not a substitute for preregistered, "
        "independently replicated publication analyses.",
        str(diversity_metadata.get("limitations", "Diversity scope was not declared.")),
    ]
    if ds.has(ArtifactKind.ANALYSIS_RESULT, "demography"):
        demo = ds.get(ArtifactKind.ANALYSIS_RESULT, "demography")
        values.append(str(demo.metadata.get("limitations", "Demographic limitations unavailable.")))
    if ds.has(ArtifactKind.QC_TABLE, "harmonization_diagnostics"):
        values.append(
            "Harmonization is limited to pre-called VCFs on a common verified build; "
            "automatic cross-assembly liftover is intentionally not run."
        )
    return list(dict.fromkeys(values))


def _executive_summary(ds, diversity_metadata: dict) -> str:
    parts: list[str] = []
    if ds.has(ArtifactKind.GENOTYPES, "genotypes"):
        meta = ds.get(ArtifactKind.GENOTYPES, "genotypes").metadata
        parts.append(
            f"The analysis contains {meta.get('n_samples', '?')} included samples and "
            f"{meta.get('n_variants', '?')} retained variants."
        )
    if ds.has(ArtifactKind.ANALYSIS_RESULT, "dstats"):
        parts.append("D-statistics were calculated with uncertainty reported in the introgression section.")
    if ds.has(ArtifactKind.ANALYSIS_RESULT, "local_ancestry"):
        parts.append("Chromosome-reset local ancestry calls are included where source panels were configured.")
    if diversity_metadata.get("scope") == "panel_relative":
        parts.append("Diversity values are panel-relative, not per-base whole-genome estimates.")
    return " ".join(parts) or "The configured analyses completed; inspect each section and its limits."


def _run_status(ctx: RunContext) -> str:
    if ctx.run_dir is None:
        return "Run manifest location unavailable."
    return f"Provenance, parameters, cache decisions, and recovery instructions are recorded in {ctx.run_dir}."


def _downloads(ds, report_dir: Path) -> list[dict[str, str]]:
    outputs: list[dict[str, str]] = []
    for artifact in sorted(ds.all(), key=lambda a: (a.kind.value, a.role)):
        try:
            href = os.path.relpath(artifact.path, report_dir).replace("\\", "/")
        except ValueError:  # different Windows drive; retain a useful local path label
            href = str(artifact.path)
        outputs.append({
            "href": href,
            "label": f"{artifact.role}: {artifact.path.name}",
            "kind": artifact.kind.value,
        })
    return outputs


def _methods(ctx: RunContext, diversity_metadata: dict) -> str:
    scope = diversity_metadata.get("scope", "not declared")
    return (
        "Configured CANIS stages exchange typed, checksummed artifacts. QC removes failed "
        "samples/sites before materialization; cache fingerprints bind inputs, configuration, and "
        f"stage code. Diversity scope is {scope}; parameters are identified by config digest "
        f"{ctx.config.digest()[:16]}. Source accessions and source paths are retained in the "
        "artifact manifest."
    )


def _sources_table(ds) -> str:
    rows: list[dict[str, object]] = []
    for artifact in ds.find(ArtifactKind.CALLSET):
        metadata = artifact.metadata
        rows.append({
            "artifact_role": artifact.role,
            "dataset_or_accession": metadata.get("dataset_id", "not supplied"),
            "reference_build": metadata.get("reference_build", metadata.get("reference_id", "unknown")),
            "source": metadata.get("source", metadata.get("vcf_url", str(artifact.path))),
            "manifest_or_checksum": metadata.get("manifest", metadata.get("sha256", "see run manifest")),
        })
    return pd.DataFrame(rows).to_html(index=False, border=0) if rows else ""


def _pca_svg(pca_csv: Path) -> str:
    """A dependency-free inline PCA scatter with native hover/focus tooltips."""
    df = pd.read_csv(pca_csv)
    if not {"PC1", "PC2"}.issubset(df.columns) or df.empty:
        return ""
    width, height, pad = 680, 400, 42
    x = pd.to_numeric(df["PC1"], errors="coerce")
    y = pd.to_numeric(df["PC2"], errors="coerce")
    valid = x.notna() & y.notna()
    if not valid.any():
        return ""
    x, y, shown = x[valid], y[valid], df.loc[valid]
    xlow, xhigh = float(x.min()), float(x.max())
    ylow, yhigh = float(y.min()), float(y.max())
    xspan, yspan = max(xhigh - xlow, 1e-12), max(yhigh - ylow, 1e-12)
    populations = shown.get("population", pd.Series("unknown", index=shown.index)).astype(str)
    palette = {p: figures._PALETTE[i % len(figures._PALETTE)]
               for i, p in enumerate(sorted(populations.unique()))}
    circles: list[str] = []
    for row, xv, yv, pop in zip(shown.itertuples(index=False), x, y, populations, strict=True):
        cx = pad + (float(xv) - xlow) / xspan * (width - 2 * pad)
        cy = height - pad - (float(yv) - ylow) / yspan * (height - 2 * pad)
        sample = html.escape(str(getattr(row, "sample_id", "sample")))
        tooltip = html.escape(f"{sample} | {pop} | PC1={float(xv):.4g}, PC2={float(yv):.4g}")
        circles.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="5" fill="{palette[str(pop)]}" '
            f'tabindex="0"><title>{tooltip}</title></circle>'
        )
    legend = "".join(
        f'<text x="{pad + i * 125}" y="{height - 8}" fill="{color}" font-size="12">'
        f'● {html.escape(pop)}</text>'
        for i, (pop, color) in enumerate(palette.items())
    )
    return (
        f'<svg class="interactive-pca" viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Interactive PCA scatterplot"><line x1="42" y1="358" x2="638" y2="358" '
        'stroke="#555"/><line x1="42" y1="42" x2="42" y2="358" stroke="#555"/>'
        f'<text x="300" y="392" font-size="13">PC1</text><text x="8" y="30" font-size="13">PC2</text>'
        f'{"".join(circles)}{legend}</svg>'
    )
