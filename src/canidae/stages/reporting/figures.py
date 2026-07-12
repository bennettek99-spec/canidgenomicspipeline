"""Publication-style figures for the population-genomics report.

All plotting uses the non-interactive Agg backend and a single shared theme
(colorblind-safe palette, vector-friendly) so figures look consistent across reports and
render headlessly in CI. Each function reads a result CSV and writes one PNG.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Okabe-Ito colorblind-safe qualitative palette.
_PALETTE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#999999",
]
_RC = {
    "figure.dpi": 120,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "font.size": 11,
}


def _color_map(labels: list[str]) -> dict[str, str]:
    uniq = sorted(dict.fromkeys(labels))
    return {lab: _PALETTE[i % len(_PALETTE)] for i, lab in enumerate(uniq)}


def pca_scatter(pca_csv: Path, out_png: Path, *,
                explained_variance: list[float] | None = None) -> Path:
    df = pd.read_csv(pca_csv)
    colors = _color_map(df["population"].astype(str).tolist())
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        for pop, sub in df.groupby("population"):
            ax.scatter(sub["PC1"], sub["PC2"], s=60, alpha=0.85,
                       color=colors[str(pop)], edgecolor="white", linewidth=0.5,
                       label=str(pop))
        xlab, ylab = "PC1", "PC2"
        if explained_variance and len(explained_variance) >= 2:
            xlab = f"PC1 ({explained_variance[0] * 100:.1f}%)"
            ylab = f"PC2 ({explained_variance[1] * 100:.1f}%)"
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title("Principal component analysis")
        ax.legend(title="population", fontsize=9, frameon=False)
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def fst_heatmap(fst_csv: Path, out_png: Path) -> Path:
    df = pd.read_csv(fst_csv, index_col=0)
    pops = list(df.index)
    mat = df.to_numpy(dtype=float)
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(1.2 + 0.8 * len(pops), 1.0 + 0.8 * len(pops)))
        im = ax.imshow(mat, cmap="viridis")
        ax.set_xticks(range(len(pops)), pops, rotation=45, ha="right")
        ax.set_yticks(range(len(pops)), pops)
        ax.grid(False)
        for i in range(len(pops)):
            for j in range(len(pops)):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                            color="white" if val < np.nanmax(mat) * 0.6 else "black",
                            fontsize=8)
        ax.set_title("Pairwise $F_{ST}$ (Hudson)")
        fig.colorbar(im, ax=ax, shrink=0.8, label="$F_{ST}$")
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def admixture_barplot(q_csv: Path, out_png: Path) -> Path:
    df = pd.read_csv(q_csv)
    q_cols = [c for c in df.columns if c.startswith("Q")]
    df = df.sort_values(["population", q_cols[0]]).reset_index(drop=True)
    x = np.arange(len(df))
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(max(6.0, 0.32 * len(df)), 4.0))
        bottom = np.zeros(len(df))
        for i, col in enumerate(q_cols):
            ax.bar(x, df[col], bottom=bottom, width=1.0, color=_PALETTE[i % len(_PALETTE)],
                   edgecolor="white", linewidth=0.15, label=col)
            bottom += df[col].to_numpy()
        ax.set_ylim(0, 1)
        ax.set_ylabel("ancestry proportion")
        ax.set_title(f"Admixture (K={len(q_cols)})")
        ax.grid(False)
        # population group labels + separators
        centers, bounds = [], []
        for pop in dict.fromkeys(df["population"]):
            idx = np.where(df["population"].to_numpy() == pop)[0]
            centers.append((pop, idx.mean()))
            bounds.append(idx.max() + 0.5)
        ax.set_xticks([c for _, c in centers], [p for p, _ in centers], rotation=30,
                      ha="right")
        for b in bounds[:-1]:
            ax.axvline(b, color="black", linewidth=0.8)
        ax.set_xlim(-0.5, len(df) - 0.5)
        ax.legend(title="component", fontsize=8, frameon=False, ncol=len(q_cols))
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def cluster_dendrogram(distance_csv: Path, out_png: Path) -> Path:
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform

    dist = pd.read_csv(distance_csv, index_col=0)
    z = linkage(squareform(dist.to_numpy(dtype=float), checks=False), method="average")
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(max(6.0, 0.3 * len(dist)), 4.0))
        dendrogram(z, labels=list(dist.index), ax=ax, leaf_rotation=90,
                   color_threshold=0.7 * z[:, 2].max())
        ax.set_ylabel("genetic distance")
        ax.set_title("Hierarchical clustering")
        ax.grid(False)
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def roh_barplot(roh_csv: Path, out_png: Path) -> Path:
    df = pd.read_csv(roh_csv).sort_values("froh", ascending=False)
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(max(6.0, 0.3 * len(df)), 4.0))
        ax.bar(df["sample_id"].astype(str), df["froh"], color=_PALETTE[0])
        ax.set_ylabel(r"$F_{ROH}$ (ROH fraction of genome)")
        ax.set_title("Runs of homozygosity")
        ax.tick_params(axis="x", rotation=90, labelsize=7)
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def locality_map(localities_csv: Path, out_png: Path) -> Path:
    df = pd.read_csv(localities_csv).dropna(subset=["latitude", "longitude"])
    colors = _color_map(df["population"].astype(str).tolist())
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(6.8, 4.6))
        for pop, sub in df.groupby("population"):
            ax.scatter(sub["longitude"], sub["latitude"], s=70, alpha=0.85,
                       color=colors[str(pop)], edgecolor="white", linewidth=0.5,
                       label=str(pop))
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        ax.set_title("Sample localities")
        ax.legend(title="population", fontsize=8, frameon=False)
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def ibd_scatter(ibd_csv: Path, out_png: Path, *, mantel_r: float | None = None) -> Path:
    df = pd.read_csv(ibd_csv)
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(6.0, 4.4))
        ax.scatter(df["geo_km"], df["genetic_distance"], s=18, alpha=0.5,
                   color=_PALETTE[0])
        if len(df) >= 2 and df["geo_km"].nunique() > 1:
            slope, intercept = np.polyfit(df["geo_km"], df["genetic_distance"], 1)
            xs = np.array([df["geo_km"].min(), df["geo_km"].max()])
            ax.plot(xs, slope * xs + intercept, color=_PALETTE[1], linewidth=1.5)
        ax.set_xlabel("geographic distance (km)")
        ax.set_ylabel("genetic distance")
        title = "Isolation by distance"
        if mantel_r is not None and mantel_r == mantel_r:  # not NaN
            title += f"  (Mantel r = {mantel_r:.3f})"
        ax.set_title(title)
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def nj_tree_figure(newick_path: Path, out_png: Path,
                   tip_colors: dict[str, str] | None = None) -> Path:
    from Bio import Phylo

    tree = Phylo.read(str(newick_path), "newick")
    n_tips = tree.count_terminals()
    with plt.rc_context({**_RC, "axes.grid": False}):
        fig, ax = plt.subplots(figsize=(7.0, max(3.0, 0.26 * n_tips)))
        Phylo.draw(
            tree, do_show=False, axes=ax,
            label_colors=(tip_colors or {}),
            branch_labels=lambda c: (c.name if (c.name and not c.is_terminal()) else None),
        )
        ax.set_title("Neighbor-joining tree (internal labels = bootstrap %)")
        ax.set_xlabel("genetic distance")
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def fd_scan_plot(fd_csv: Path, out_png: Path) -> Path:
    df = pd.read_csv(fd_csv)
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(7.0, 3.4))
        x = np.arange(len(df))
        ax.plot(x, df["f_d"], color=_PALETTE[1], linewidth=1.0, marker="o", markersize=3)
        ax.axhline(0.0, color="#888", linewidth=0.8)
        ax.set_ylabel(r"$f_d$")
        ax.set_xlabel("genomic window")
        ax.set_title("Windowed introgression scan ($f_d$)")
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def karyogram(windows_csv: Path, out_png: Path) -> Path:
    from matplotlib.patches import Patch

    df = pd.read_csv(windows_csv)
    contigs = sorted(df["chrom"].astype(str).unique())
    offset, cum = {}, 0.0
    for c in contigs:
        offset[c] = cum
        cum += float(df.loc[df["chrom"].astype(str) == c, "end"].max()) * 1.05
    samples = list(dict.fromkeys(df["sample_id"]))
    ancestries = sorted(df["ancestry"].astype(str).unique())
    colors = {a: _PALETTE[i % len(_PALETTE)] for i, a in enumerate(ancestries)}
    with plt.rc_context({**_RC, "axes.grid": False}):
        fig, ax = plt.subplots(figsize=(8.0, max(2.5, 0.32 * len(samples))))
        for i, s in enumerate(samples):
            sub = df[df["sample_id"] == s]
            bars = [(offset[str(r.chrom)] + r.start, r.end - r.start)
                    for r in sub.itertuples(index=False)]
            cols = [colors[str(r.ancestry)] for r in sub.itertuples(index=False)]
            ax.broken_barh(bars, (i - 0.4, 0.8), facecolors=cols)
        ax.set_yticks(range(len(samples)), samples)
        ax.set_xlabel("genome position (bp, contigs concatenated)")
        ax.set_title("Local ancestry (karyogram)")
        ax.legend(handles=[Patch(color=colors[a], label=a) for a in ancestries],
                  fontsize=8, frameon=False, title="ancestry")
        fig.savefig(out_png)
        plt.close(fig)
    return out_png


def diversity_bar(diversity_csv: Path, out_png: Path) -> Path:
    df = pd.read_csv(diversity_csv).sort_values("pi", ascending=False)
    is_panel = "diversity_scope" in df and (df["diversity_scope"] == "panel_relative").any()
    colors = _color_map(df["population"].astype(str).tolist())
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        bars = ax.bar(df["population"].astype(str), df["pi"],
                      color=[colors[str(p)] for p in df["population"]])
        ax.set_ylabel(
            "panel-relative diversity (selected sites)" if is_panel
            else r"nucleotide diversity $\pi$ (per callable site)"
        )
        ax.set_xlabel("population")
        ax.set_title("Within-population diversity" + (" (panel-relative)" if is_panel else ""))
        ax.tick_params(axis="x", rotation=30)
        for rect, val in zip(bars, df["pi"], strict=False):
            ax.text(rect.get_x() + rect.get_width() / 2, val, f"{val:.3g}",
                    ha="center", va="bottom", fontsize=8)
        fig.savefig(out_png)
        plt.close(fig)
    return out_png
