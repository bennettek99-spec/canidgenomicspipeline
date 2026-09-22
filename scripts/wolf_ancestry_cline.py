"""Red wolf, Great Lakes wolf and eastern wolf: a gray-wolf-to-coyote ancestry cline.

Reads the window panel written by ``scripts/fetch_wolf_cline_panel.py`` and asks one
question with several independent statistics: how much of each target genome traces
to the gray-wolf lineage versus the coyote lineage?

* f4-ratio  alpha = f4(Eurasian wolf, jackal; X, coyote) / f4(Eurasian wolf, jackal;
  western gray wolf, coyote) -> gray-wolf ancestry fraction, chromosome jackknife SE.
* D(X, western gray wolf; coyote, jackal) -> does X share excess drift with coyotes?
* PCA and unsupervised two-component sNMF of the North American canids.
* Heterozygosity of the high-coverage (>= 15x) genomes.

Every f4 / D statistic is computed twice: on GQ >= 20 calls and on all hard calls.
Low-coverage genomes lose most homozygous-reference calls at GQ >= 20, so agreement
between the two is the robustness check. Reference controls (western wolves, coyotes)
are scored leave-one-out, so each is never part of its own reference group.

Outputs (small, committed) go to ``results/wolf_ancestry_cline/``.

Run:
    python scripts/wolf_ancestry_cline.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from canidae.analysis.f_statistics import (
    JackknifeEstimate,
    d_statistic,
    f4_ratio,
    population_frequencies,
)
from canidae.stages.popgen.ancestry_nmf import weighted_nmf

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "data" / "wolf_cline" / "panel.npz"
SAMPLES_CSV = ROOT / "configs" / "examples" / "wolf_cline_samples.csv"
OUT = ROOT / "results" / "wolf_ancestry_cline"

MIN_GQ = 20
HIGH_COVERAGE = 15.0
TARGET_GROUPS = ["great_lakes_wolf", "eastern_wolf", "red_wolf"]
NA_ORDER = ["western_gray_wolf", "great_lakes_wolf", "eastern_wolf", "red_wolf", "coyote"]
GROUP_ORDER = ["eurasian_wolf", *NA_ORDER, "golden_jackal"]
LABELS = {
    "eurasian_wolf": "Eurasian gray wolf",
    "western_gray_wolf": "Western gray wolf",
    "great_lakes_wolf": "Great Lakes wolf",
    "eastern_wolf": "Eastern wolf",
    "algonquin": "Eastern wolf (Algonquin only)",
    "red_wolf": "Red wolf",
    "coyote": "Coyote",
    "golden_jackal": "Golden jackal",
}
COLORS = {
    "eurasian_wolf": "#6b7a8f",
    "western_gray_wolf": "#2f5d8a",
    "great_lakes_wolf": "#2a9d8f",
    "eastern_wolf": "#8a6fb3",
    "algonquin": "#8a6fb3",
    "red_wolf": "#c8553d",
    "coyote": "#d9a13b",
    "golden_jackal": "#8c8c8c",
}
WOLF_COLOR, COYOTE_COLOR = "#2f5d8a", "#d9a13b"
CALLSETS = {"gq20": f"GQ ≥ {MIN_GQ}", "all": "all hard calls"}


def load_samples() -> dict[str, dict[str, str]]:
    with SAMPLES_CSV.open(newline="", encoding="utf-8") as handle:
        return {row["sample_id"]: row for row in csv.DictReader(handle)}


def short(sample: str) -> str:
    return sample.replace("AlgonquinWolf", "Algonquin")


def pca(genotypes: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Patterson-normalized PCA (sites x samples in); missing calls mean-imputed."""
    g = genotypes.astype(float)
    g[g < 0] = np.nan
    mean = np.nanmean(g, axis=1, keepdims=True)
    p = mean / 2
    x = np.nan_to_num((g - mean) / np.sqrt(p * (1 - p))).T
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    explained = s**2 / np.sum(s**2)
    return u[:, :n_components] * s[:n_components], explained[:n_components]


def snmf(genotypes: np.ndarray, k: int, restarts: int = 5) -> np.ndarray:
    """Best-of-restarts masked sNMF ancestry proportions (samples x k)."""
    x = genotypes.T.astype(float)
    observed = (x >= 0).astype(float)
    x = np.where(observed > 0, x, 0.0)
    best_err, best_q = np.inf, np.zeros((x.shape[0], k))
    for seed in range(restarts):
        w, h = weighted_nmf(x, observed, k, n_iter=500, seed=seed)
        err = float(np.sum(observed * (x - w @ h) ** 2))
        if err < best_err:
            best_err, best_q = err, w / w.sum(axis=1, keepdims=True)
    return best_q


def estimate_fields(est: JackknifeEstimate, suffix: str) -> dict[str, Any]:
    return {
        f"estimate_{suffix}": round(est.estimate, 4),
        f"se_{suffix}": round(est.se, 4),
        f"z_{suffix}": round(est.z, 2),
        f"n_sites_{suffix}": est.n_sites,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = np.load(PANEL)
    samples = [str(s) for s in panel["samples"]]
    meta = load_samples()
    if set(samples) != set(meta):
        raise SystemExit("panel samples do not match the sample sheet; re-run the fetch")

    raw = panel["gt"].astype(np.int8)
    filtered = raw.copy()
    filtered[panel["gq"] < MIN_GQ] = -1
    callsets = {"gq20": filtered, "all": raw}
    blocks = panel["chrom"].astype(int)  # 38 autosome jackknife blocks
    groups: dict[str, list[str]] = {g: [] for g in GROUP_ORDER}
    for sample in samples:
        groups[meta[sample]["population"]].append(sample)
    col = {s: i for i, s in enumerate(samples)}
    spans = panel["spans"]
    total_span = int(np.sum(np.maximum(spans[:, 3] - spans[:, 2], 0)))

    # Units scored: target groups (plus the Algonquin-only eastern wolves) and every
    # target/control individual. Controls are scored leave-one-out.
    units: list[tuple[str, str, str, list[str]]] = []
    for group in TARGET_GROUPS:
        units.append(("group", group, group, groups[group]))
        if group == "eastern_wolf":
            algonquin = [s for s in groups[group] if s.startswith("Algonquin")]
            units.append(("group", "algonquin", group, algonquin))
    for group in NA_ORDER:
        units.extend(("individual", s, group, [s]) for s in groups[group])

    ancestry_rows: list[dict[str, Any]] = []
    d_rows: list[dict[str, Any]] = []
    for unit, unit_id, group, members in units:
        wolves = [s for s in groups["western_gray_wolf"] if s not in members]
        coyotes = [s for s in groups["coyote"] if s not in members]
        pops = {
            "A": groups["eurasian_wolf"],
            "O": groups["golden_jackal"],
            "B": wolves,
            "C": coyotes,
            "X": members,
        }
        a_row: dict[str, Any] = {"unit": unit, "id": unit_id, "group": group}
        d_row: dict[str, Any] = {"unit": unit, "id": unit_id, "group": group}
        for key, gt in callsets.items():
            freqs = population_frequencies(gt, samples, pops)
            a_row |= estimate_fields(f4_ratio(freqs, "A", "O", "X", "B", "C", blocks), key)
            # D(X, W; C, O) > 0: X shares more derived alleles with coyotes than western wolves do.
            d_row |= estimate_fields(d_statistic(freqs, "X", "B", "C", "O", blocks), key)
        ancestry_rows.append(a_row)
        if group in TARGET_GROUPS:
            d_rows.append(d_row)

    # Between-target contrasts on the high-quality calls.
    freqs = population_frequencies(
        filtered,
        samples,
        {
            # Algonquin only: the Quebec genome is coyote-like (see README).
            "E": [s for s in groups["eastern_wolf"] if s.startswith("Algonquin")],
            "G": groups["great_lakes_wolf"],
            "R": groups["red_wolf"],
            "C": groups["coyote"],
            "O": groups["golden_jackal"],
        },
    )
    contrasts = {
        "D(red_wolf, algonquin; coyote, jackal)": d_statistic(freqs, "R", "E", "C", "O", blocks),
        "D(algonquin, great_lakes_wolf; coyote, jackal)": d_statistic(
            freqs, "E", "G", "C", "O", blocks
        ),
        "D(red_wolf, great_lakes_wolf; coyote, jackal)": d_statistic(
            freqs, "R", "G", "C", "O", blocks
        ),
    }

    # PCA + sNMF on North American canids: all hard calls, >= 90% call rate, MAF >= 5%,
    # thinned to one SNP per 2 kb to damp local LD.
    na_samples = [s for g in NA_ORDER for s in groups[g]]
    na = raw[:, [col[s] for s in na_samples]]
    called = na >= 0
    alt = np.where(called, na, 0).sum(1) / np.maximum(2 * called.sum(1), 1)
    keep = (called.mean(1) >= 0.9) & (np.minimum(alt, 1 - alt) >= 0.05)
    bins = blocks.astype(np.int64) * 10**6 + panel["pos"] // 2000
    kept_bins = bins[keep]
    thin = np.zeros_like(keep)
    thin[np.flatnonzero(keep)[np.r_[True, kept_bins[1:] != kept_bins[:-1]]]] = True
    pcs, explained = pca(na[thin])
    q2 = snmf(na[thin], 2)
    wolf_comp = int(
        np.argmax(q2[[na_samples.index(s) for s in groups["western_gray_wolf"]]].mean(0))
    )
    snmf_wolf = q2[:, wolf_comp]

    # Heterozygosity: fraction of called panel sites that are heterozygous, scaled by
    # panel-site density per bp of window span. Reliable only for high coverage.
    density = raw.shape[0] / total_span
    het_rows: list[dict[str, Any]] = []
    for s in samples:
        calls = raw[:, col[s]]
        het = float((calls == 1).sum() / (calls >= 0).sum())
        coverage = float(meta[s]["coverage"])
        het_rows.append(
            {
                "sample_id": s,
                "group": meta[s]["population"],
                "locality": meta[s]["locality"],
                "coverage": coverage,
                "missing_all_calls": round(float((calls < 0).mean()), 4),
                "missing_gq20": round(float((filtered[:, col[s]] < 0).mean()), 4),
                "het_per_kb": round(het * density * 1000, 3),
                "high_coverage": coverage >= HIGH_COVERAGE,
            }
        )

    def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(OUT / "f4_ratio_wolf_ancestry.csv", ancestry_rows)
    write_csv(OUT / "d_statistics.csv", d_rows)
    write_csv(OUT / "sample_qc_heterozygosity.csv", het_rows)
    write_csv(
        OUT / "pca_snmf.csv",
        [
            {
                "sample_id": s,
                "group": meta[s]["population"],
                "PC1": round(float(pcs[i, 0]), 4),
                "PC2": round(float(pcs[i, 1]), 4),
                "snmf_k2_wolf": round(float(snmf_wolf[i]), 4),
                "snmf_k2_coyote": round(float(1 - snmf_wolf[i]), 4),
            }
            for i, s in enumerate(na_samples)
        ],
    )
    summary: dict[str, Any] = {
        "source": str(panel["source"]),
        "reference": "CanFam3.1",
        "bytes_transferred": int(panel["bytes_transferred"]),
        "n_windows": len(spans),
        "window_span_bp": total_span,
        "n_snps": int(raw.shape[0]),
        "jackknife_blocks": "38 autosomes",
        "f4_ratio_model": "alpha = f4(eurasian_wolf, golden_jackal; X, coyote) / "
        "f4(eurasian_wolf, golden_jackal; western_gray_wolf, coyote)",
        "wolf_ancestry": {
            r["id"]: {
                "gq20": [r["estimate_gq20"], r["se_gq20"]],
                "all_calls": [r["estimate_all"], r["se_all"]],
            }
            for r in ancestry_rows
        },
        "d_x_westernwolf_coyote_jackal": {
            r["id"]: {
                "gq20": [r["estimate_gq20"], r["z_gq20"]],
                "all_calls": [r["estimate_all"], r["z_all"]],
            }
            for r in d_rows
        },
        "target_contrasts_gq20": {
            name: {"D": round(e.estimate, 4), "se": round(e.se, 4), "z": round(e.z, 2)}
            for name, e in contrasts.items()
        },
        "pca": {
            "n_sites": int(thin.sum()),
            "explained_variance_ratio": [round(float(v), 4) for v in explained],
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plot_figure(
        ancestry_rows, d_rows, het_rows, na_samples, meta, pcs, explained, snmf_wolf, summary
    )
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "wolf_ancestry",
                    "d_x_westernwolf_coyote_jackal",
                    "target_contrasts_gq20",
                    "pca",
                )
            },
            indent=1,
        )
    )


def plot_figure(
    ancestry_rows: list[dict[str, Any]],
    d_rows: list[dict[str, Any]],
    het_rows: list[dict[str, Any]],
    na_samples: list[str],
    meta: dict[str, dict[str, str]],
    pcs: np.ndarray,
    explained: np.ndarray,
    snmf_wolf: np.ndarray,
    summary: dict[str, Any],
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.titlesize": 10,
            "axes.titlelocation": "left",
            "axes.titlepad": 10,
        }
    )
    fig = plt.figure(figsize=(14, 12))
    grid = fig.add_gridspec(3, 3, height_ratios=[1.1, 1, 1], hspace=0.62, wspace=0.34)
    indiv = [r for r in ancestry_rows if r["unit"] == "individual"]
    indiv.sort(key=lambda r: (NA_ORDER.index(r["group"]), -r["estimate_gq20"]))

    def color_ticks(ax: Any, ids: list[str]) -> None:
        for tick, sample in zip(ax.get_xticklabels(), ids, strict=True):
            tick.set_color(COLORS[meta[sample]["population"]])

    # (a) Per-genome ancestry bars.
    ax = fig.add_subplot(grid[0, :2])
    x = np.arange(len(indiv))
    wolf = np.clip([r["estimate_gq20"] for r in indiv], 0, 1)
    ax.bar(x, wolf, color=WOLF_COLOR, width=0.82, label="Gray-wolf lineage")
    ax.bar(x, 1 - wolf, bottom=wolf, color=COYOTE_COLOR, width=0.82, label="Coyote lineage")
    ax.errorbar(
        x,
        [r["estimate_gq20"] for r in indiv],
        yerr=[1.96 * r["se_gq20"] for r in indiv],
        fmt="none",
        ecolor="#1d1d1d",
        elinewidth=1,
        capsize=2,
    )
    ids = [r["id"] for r in indiv]
    ax.set_xticks(x, [short(s) for s in ids], rotation=50, ha="right", fontsize=7.5)
    color_ticks(ax, ids)
    edges = [i for i in range(1, len(indiv)) if indiv[i]["group"] != indiv[i - 1]["group"]]
    for edge in edges:
        ax.axvline(edge - 0.5, color="white", lw=2.5)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(-0.6, len(indiv) - 0.4)
    ax.set_ylabel("Ancestry fraction (f4-ratio)")
    ax.set_title("a  Gray-wolf vs coyote ancestry of each genome (±95% CI)")
    ax.legend(loc="lower left", frameon=True, framealpha=0.9, fontsize=8)

    # (b) The cline: group estimates with individuals behind them.
    ax = fig.add_subplot(grid[0, 2])
    rows = [
        "western_gray_wolf",
        "great_lakes_wolf",
        "algonquin",
        "eastern_wolf",
        "red_wolf",
        "coyote",
    ]
    groups_by_id = {r["id"]: r for r in ancestry_rows if r["unit"] == "group"}
    ys = np.arange(len(rows))[::-1]
    for y, row in zip(ys, rows, strict=True):
        members = (
            [r for r in indiv if r["group"] == row]
            if row != "algonquin"
            else [r for r in indiv if r["id"].startswith("Algonquin")]
        )
        ax.scatter(
            [r["estimate_gq20"] for r in members],
            [y] * len(members),
            s=18,
            color=COLORS[row],
            alpha=0.45,
            zorder=2,
        )
        if row in groups_by_id:
            g = groups_by_id[row]
            ax.errorbar(
                g["estimate_gq20"],
                y,
                xerr=1.96 * g["se_gq20"],
                fmt="D",
                ms=8,
                color=COLORS[row],
                capsize=3,
                zorder=3,
                mec="white",
            )
            ax.annotate(
                f"{g['estimate_gq20'] * 100:.0f}% wolf",
                (g["estimate_gq20"], y),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=8,
            )
    ax.set_yticks(
        ys, [LABELS[r].replace(" (Algonquin only)", "\n(Algonquin only)") for r in rows], fontsize=8
    )
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.axvline(0, color="#cccccc", lw=0.8, ls=":")
    ax.axvline(1, color="#cccccc", lw=0.8, ls=":")
    ax.set_xlabel("Gray-wolf ancestry")
    ax.set_title("b  The wolf → coyote cline")
    ax.text(
        0.99,
        0.01,
        "◆ group estimate  ● individual\ncontrols scored leave-one-out",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color="#555555",
    )

    # (c) PCA.
    ax = fig.add_subplot(grid[1, 0])
    for group in NA_ORDER:
        idx = [i for i, s in enumerate(na_samples) if meta[s]["population"] == group]
        ax.scatter(
            pcs[idx, 0],
            pcs[idx, 1],
            s=52,
            color=COLORS[group],
            label=LABELS[group],
            edgecolor="white",
            lw=0.7,
            zorder=3,
        )
    for i, s in enumerate(na_samples):
        if s in ("QuebecWolf", "Wolf40"):
            ax.annotate(
                short(s),
                (pcs[i, 0], pcs[i, 1]),
                xytext=(5, 4),
                textcoords="offset points",
                fontsize=7,
            )
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax.set_title(f"c  PCA ({summary['pca']['n_sites']:,} LD-thinned SNPs)")
    ax.legend(fontsize=7, frameon=False, loc="best")

    # (d) sNMF K=2 alongside the f4-ratio.
    ax = fig.add_subplot(grid[1, 1:])
    order = [na_samples.index(i) for i in ids]
    xs = np.arange(len(order))
    q = snmf_wolf[order]
    ax.bar(xs, q, color=WOLF_COLOR, width=0.82)
    ax.bar(xs, 1 - q, bottom=q, color=COYOTE_COLOR, width=0.82)
    ax.scatter(
        xs,
        np.clip([r["estimate_gq20"] for r in indiv], 0, 1),
        marker="_",
        s=260,
        color="white",
        linewidths=2.2,
        zorder=3,
        label="f4-ratio estimate (panel a)",
    )
    ax.set_xticks(xs, [short(s) for s in ids], rotation=50, ha="right", fontsize=7.5)
    color_ticks(ax, ids)
    for edge in edges:
        ax.axvline(edge - 0.5, color="white", lw=2.5)
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.6, len(ids) - 0.4)
    ax.set_ylabel("Ancestry (Q)")
    ax.set_title("d  Independent check: unsupervised sNMF, K = 2 (no reference groups)")
    legend = ax.legend(
        loc="lower left", fontsize=7.5, frameon=True, facecolor="#555555", framealpha=0.85
    )
    for text in legend.get_texts():
        text.set_color("white")

    # (e) D-statistics.
    ax = fig.add_subplot(grid[2, 0])
    drows = [r for r in d_rows if r["unit"] == "individual"]
    ys = np.arange(len(drows))[::-1]
    for y, r in zip(ys, drows, strict=True):
        ax.errorbar(
            r["estimate_gq20"],
            y,
            xerr=1.96 * r["se_gq20"],
            fmt="o",
            color=COLORS[r["group"]],
            ms=6,
            capsize=2,
        )
        ax.annotate(
            f"Z={r['z_gq20']:.1f}",
            (r["estimate_gq20"] + 1.96 * r["se_gq20"], y),
            xytext=(4, -3),
            textcoords="offset points",
            fontsize=6.5,
            color="#555555",
        )
    ax.axvline(0, color="#999999", lw=0.8)
    ax.set_yticks(ys, [short(r["id"]) for r in drows], fontsize=7.5)
    ax.set_xlim(-0.02, 0.42)
    ax.set_xlabel("D(X, western wolf; coyote, jackal)")
    ax.set_title("e  Excess allele sharing with coyotes")

    # (f) Call-quality sensitivity.
    ax = fig.add_subplot(grid[2, 1])
    for r in indiv:
        low = float(meta[r["id"]]["coverage"]) < 10
        ax.scatter(
            r["estimate_gq20"],
            r["estimate_all"],
            s=40,
            color=COLORS[r["group"]],
            marker="o" if not low else "^",
            edgecolor="white",
            lw=0.5,
            zorder=3,
        )
    ax.plot([-0.1, 1.1], [-0.1, 1.1], color="#999999", lw=0.8, ls="--")
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.12, 1.12)
    ax.set_xlabel(f"Gray-wolf ancestry, GQ ≥ {MIN_GQ}")
    ax.set_ylabel("Gray-wolf ancestry, all hard calls")
    ax.set_title("f  Robustness to genotype-quality filter")
    ax.text(
        0.98,
        0.03,
        "● ≥10x coverage   ▲ <10x",
        transform=ax.transAxes,
        ha="right",
        fontsize=7,
        color="#555555",
    )

    # (g) Heterozygosity, high-coverage genomes only.
    ax = fig.add_subplot(grid[2, 2])
    hrows = [r for r in het_rows if r["high_coverage"] and r["group"] != "golden_jackal"]
    hrows.sort(key=lambda r: (GROUP_ORDER.index(r["group"]), r["sample_id"]))
    ys = np.arange(len(hrows))[::-1]
    ax.barh(
        ys, [r["het_per_kb"] for r in hrows], color=[COLORS[r["group"]] for r in hrows], height=0.72
    )
    places = [
        r["locality"].split(" (")[0].split(",")[0].replace(" National Park", "") for r in hrows
    ]
    ax.set_yticks(
        ys, [f"{r['sample_id']} · {p}" for r, p in zip(hrows, places, strict=True)], fontsize=7
    )
    for y, r in zip(ys, hrows, strict=True):
        ax.text(r["het_per_kb"] + 0.02, y, f"{r['het_per_kb']:.2f}", va="center", fontsize=6.5)
    ax.set_xlabel("Heterozygous sites per kb")
    ax.set_title(f"g  Heterozygosity (genomes ≥ {HIGH_COVERAGE:.0f}x)")

    fig.suptitle(
        "Red wolf, Great Lakes wolf and eastern wolf on the gray-wolf ↔ coyote continuum",
        fontsize=14,
        fontweight="bold",
        x=0.01,
        ha="left",
        y=0.985,
    )
    fig.text(
        0.01,
        0.958,
        f"NHGRI 722-genome WGS callset (CanFam3.1) · {summary['n_snps']:,} PASS biallelic "
        f"SNPs from {summary['n_windows']} unascertained windows on 38 autosomes · "
        "chromosome block jackknife\n"
        "f4-ratio = f4(Eurasian wolves, golden jackal; X, coyotes) / "
        "f4(Eurasian wolves, golden jackal; western gray wolves, coyotes)",
        fontsize=8.5,
        color="#555555",
        va="top",
    )
    fig.savefig(OUT / "wolf_ancestry_cline.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "wolf_ancestry_cline.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
