"""Robustness tests for the wolf ancestry cline (companion to wolf_ancestry_cline.py).

All statistics use allele frequencies estimated by EM from genotype likelihoods.

1. Dog ancestry. Dogs sit inside the gray-wolf clade, so any dog ancestry would be
   counted as "gray wolf" by the f4-ratio. qpAdm fits each target as
   western gray wolf + dog + coyote, and as western gray wolf + coyote alone.
2. Where did the coyote ancestry come from? Recent hybridization predicts that a
   target's coyote component matches the coyotes living near it; an old, separate
   lineage predicts it matches no living coyote better than another. qpAdm is run
   with each coyote genome as the source and the other two as outgroups, alongside
   D(coyote i, coyote j; X, jackal).
3. Reference rotation. The f4-ratio is recomputed under every combination of
   Eurasian reference, North American wolf reference, coyote reference and outgroup.

Outputs go to ``results/wolf_ancestry_cline/robustness/``.

Run (after scripts/fetch_wolf_cline_panel.py):
    python scripts/wolf_cline_robustness.py
"""

from __future__ import annotations

import csv
import itertools
import json
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wolf_ancestry_cline import COLORS, LABELS, OUT, PANEL, load_samples, short

from canidae.analysis.f_statistics import (
    QpAdmResult,
    d_statistic,
    f4_ratio,
    gl_population_frequencies,
    qpadm,
)

ROBUST = OUT / "robustness"
WOLF_COLOR, DOG_COLOR, COYOTE_COLOR = "#2f5d8a", "#b07aa1", "#d9a13b"
RIGHTS = ["golden_jackal", "eurasian_wolf", "asian_wolf", "village_dog", "dhole"]
COYOTES = {"Coyote01": "California", "Coyote02": "Alabama", "Coyote03": "Midwest"}


def main() -> None:
    ROBUST.mkdir(parents=True, exist_ok=True)
    panel = np.load(PANEL)
    samples = [str(s) for s in panel["samples"]]
    meta = load_samples()
    blocks = panel["chrom"].astype(int)
    groups: dict[str, list[str]] = {}
    for sample in samples:
        groups.setdefault(meta[sample]["population"], []).append(sample)
    algonquin = [s for s in groups["eastern_wolf"] if s.startswith("Algonquin")]

    targets: dict[str, list[str]] = {
        "great_lakes_wolf": groups["great_lakes_wolf"],
        "algonquin": algonquin,
        "red_wolf": groups["red_wolf"],
    }
    individuals = [*groups["great_lakes_wolf"], *groups["eastern_wolf"], *groups["red_wolf"]]
    targets |= {s: [s] for s in individuals}

    # Controls with a known answer: a western wolf scored leave-one-out, Mexican wolves.
    controls: dict[str, list[str]] = {"Wolf28": ["Wolf28"], "mexican_wolf": groups["mexican_wolf"]}

    pops: dict[str, list[str]] = {**groups, **targets, **controls, **{c: [c] for c in COYOTES}}
    pops["europe_wolf"] = ["Wolf03", "Wolf06", "Wolf24", "Wolf27"]
    pops["asia_wolf_all"] = ["Wolf01", "Wolf02", *groups["asian_wolf"]]
    pops["yellowstone"] = ["Wolf28", "Wolf29", "Wolf30"]
    pops["AlaskanWolf"] = ["AlaskanWolf"]
    freqs = gl_population_frequencies(panel["pl"], samples, pops)

    def without(source: str, target: str) -> str:
        """Name of ``source`` with ``target``'s genomes removed (leave-one-out)."""
        overlap = set(pops[target]) & set(pops[source])
        if not overlap:
            return source
        name = f"{source}-minus-{target}"
        if name not in freqs:
            pops[name] = [m for m in pops[source] if m not in overlap]
            freqs.update(gl_population_frequencies(panel["pl"], samples, {name: pops[name]}))
        return name

    def label(target: str) -> str:
        if target == "Wolf28":
            return "Wolf28 (Yellowstone, control)"
        if target == "mexican_wolf":
            return "Mexican wolves (control)"
        if target in COYOTES:
            return f"{COYOTES[target]} coyote"
        return LABELS.get(target, short(target))

    def color(target: str) -> str:
        return COLORS[target] if target in COLORS else COLORS[meta[target]["population"]]

    def safe_qpadm(target: str, sources: list[str], rights: list[str]) -> QpAdmResult | None:
        try:
            return qpadm(freqs, target, sources, rights, blocks)
        except ValueError:
            return None

    # ---- 1. Dog ancestry ----------------------------------------------------
    dog_rows: list[dict[str, Any]] = []
    for target in [*controls, *targets]:
        wolf_source = without("western_gray_wolf", target)
        three = safe_qpadm(target, [wolf_source, "dog", "coyote"], RIGHTS)
        two = safe_qpadm(target, [wolf_source, "coyote"], RIGHTS)
        assert three is not None and two is not None
        dog_rows.append(
            {
                "target": target,
                "wolf": round(float(three.weights[0]), 4),
                "wolf_se": round(float(three.se[0]), 4),
                "dog": round(float(three.weights[1]), 4),
                "dog_se": round(float(three.se[1]), 4),
                "dog_z": round(float(three.weights[1] / three.se[1]), 2),
                "coyote": round(float(three.weights[2]), 4),
                "coyote_se": round(float(three.se[2]), 4),
                "p_wolf_dog_coyote": round(three.p_value, 4),
                "p_wolf_coyote": round(two.p_value, 4),
                "two_way_wolf": round(float(two.weights[0]), 4),
                "two_way_wolf_se": round(float(two.se[0]), 4),
                "n_sites": three.n_sites,
            }
        )

    # ---- 2. Coyote source ---------------------------------------------------
    source_rows: list[dict[str, Any]] = []
    for target in [*controls, *targets]:
        for coyote in COYOTES:
            others = [c for c in COYOTES if c != coyote]
            wolf_source = without("western_gray_wolf", target)
            fit = safe_qpadm(target, [wolf_source, coyote], [*RIGHTS, *others])
            assert fit is not None
            source_rows.append(
                {
                    "target": target,
                    "coyote_source": coyote,
                    "coyote_locality": COYOTES[coyote],
                    "wolf": round(float(fit.weights[0]), 4),
                    "wolf_se": round(float(fit.se[0]), 4),
                    "chisq": round(fit.chisq, 2),
                    "dof": fit.dof,
                    "p_value": round(fit.p_value, 4),
                }
            )
    # Is any coyote genome itself admixed? Model each as gray wolf + one other coyote.
    coyote_target_rows: list[dict[str, Any]] = []
    for coyote in COYOTES:
        for source in COYOTES:
            if source == coyote:
                continue
            other = next(c for c in COYOTES if c not in (coyote, source))
            fit = safe_qpadm(coyote, ["western_gray_wolf", source], [*RIGHTS, other])
            assert fit is not None
            coyote_target_rows.append(
                {
                    "coyote": coyote,
                    "locality": COYOTES[coyote],
                    "coyote_source": source,
                    "wolf_weight": round(float(fit.weights[0]), 4),
                    "wolf_se": round(float(fit.se[0]), 4),
                    "p_value": round(fit.p_value, 4),
                }
            )

    d_rows: list[dict[str, Any]] = []
    d_targets = [
        "western_gray_wolf",
        "eurasian_wolf",
        "great_lakes_wolf",
        "algonquin",
        "red_wolf",
        "QuebecWolf",
    ]
    for c1, c2 in itertools.combinations(COYOTES, 2):
        for target in d_targets:
            est = d_statistic(freqs, c1, c2, target, "golden_jackal", blocks)
            d_rows.append(
                {
                    "test": f"D({COYOTES[c1]}, {COYOTES[c2]}; X, jackal)",
                    "coyote_1": c1,
                    "coyote_2": c2,
                    "target": target,
                    "D": round(est.estimate, 4),
                    "se": round(est.se, 4),
                    "z": round(est.z, 2),
                }
            )

    # ---- 3. Reference rotation ----------------------------------------------
    a_sets = {"all Eurasian": "eurasian_wolf", "Europe": "europe_wolf", "Asia": "asia_wolf_all"}
    b_sets = {
        "Yellowstone + Alaska": "western_gray_wolf",
        "Yellowstone": "yellowstone",
        "Alaska": "AlaskanWolf",
        "Mexican": "mexican_wolf",
    }
    c_sets = {"all coyotes": "coyote", **{f"{v} coyote": k for k, v in COYOTES.items()}}
    o_sets = {"golden jackal": "golden_jackal", "dhole": "dhole"}
    rotation_rows: list[dict[str, Any]] = []
    rotation_targets = ["great_lakes_wolf", "algonquin", "red_wolf", "QuebecWolf"]
    for target in rotation_targets:
        for (a_name, a), (b_name, b), (c_name, c), (o_name, o) in itertools.product(
            a_sets.items(), b_sets.items(), c_sets.items(), o_sets.items()
        ):
            est = f4_ratio(freqs, a, o, target, b, c, blocks)
            rotation_rows.append(
                {
                    "target": target,
                    "A": a_name,
                    "B": b_name,
                    "C": c_name,
                    "O": o_name,
                    "wolf_ancestry": round(est.estimate, 4),
                    "se": round(est.se, 4),
                    "baseline": (a, b, c, o)
                    == ("eurasian_wolf", "western_gray_wolf", "coyote", "golden_jackal"),
                }
            )

    # ---- Write --------------------------------------------------------------
    def write_csv(name: str, rows: list[dict[str, Any]]) -> None:
        with (ROBUST / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    write_csv("qpadm_dog_ancestry.csv", dog_rows)
    write_csv("qpadm_coyote_source.csv", source_rows)
    write_csv("qpadm_coyote_as_target.csv", coyote_target_rows)
    write_csv("d_coyote_pairs.csv", d_rows)
    write_csv("reference_rotation.csv", rotation_rows)
    rotation_summary = {}
    for target in rotation_targets:
        values = np.array([r["wolf_ancestry"] for r in rotation_rows if r["target"] == target])
        by_factor = {}
        for factor in ("A", "B", "C", "O"):
            levels = sorted({r[factor] for r in rotation_rows})
            by_factor[factor] = {
                level: round(
                    float(
                        np.median(
                            [
                                r["wolf_ancestry"]
                                for r in rotation_rows
                                if r["target"] == target and r[factor] == level
                            ]
                        )
                    ),
                    4,
                )
                for level in levels
            }
        rotation_summary[target] = {
            "n_models": len(values),
            "median": round(float(np.median(values)), 4),
            "min": round(float(values.min()), 4),
            "max": round(float(values.max()), 4),
            "median_by_factor": by_factor,
        }
    # How often does the ordering Great Lakes > Algonquin > red wolf hold across rotations?
    keyed: dict[tuple[str, str, str, str], dict[str, float]] = {}
    for r in rotation_rows:
        keyed.setdefault((r["A"], r["B"], r["C"], r["O"]), {})[r["target"]] = r["wolf_ancestry"]
    ordered = [
        v["great_lakes_wolf"] > v["algonquin"] > v["red_wolf"] > v["QuebecWolf"]
        for v in keyed.values()
    ]
    summary = {
        "ordering_holds_in_models": f"{sum(ordered)}/{len(ordered)}",
        "coyote_as_target": coyote_target_rows,
        "frequencies": "EM from genotype likelihoods",
        "qpadm_rights": RIGHTS,
        "dog_ancestry": {r["target"]: r for r in dog_rows},
        "coyote_source_p_values": {
            t: {r["coyote_locality"]: r["p_value"] for r in source_rows if r["target"] == t}
            for t in [*controls, *targets]
        },
        "reference_rotation": rotation_summary,
    }
    (ROBUST / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    plot(dog_rows, source_rows, d_rows, rotation_rows, label, color)
    print(json.dumps(summary, indent=1))


def plot(
    dog_rows: list[dict[str, Any]],
    source_rows: list[dict[str, Any]],
    d_rows: list[dict[str, Any]],
    rotation_rows: list[dict[str, Any]],
    label: Any,
    color: Any,
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
    fig = plt.figure(figsize=(14, 10.5))
    grid = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.3)

    # (a) Three-way qpAdm.
    ax = fig.add_subplot(grid[0, 0])
    rows = dog_rows
    x = np.arange(len(rows))
    wolf = np.array([r["wolf"] for r in rows])
    dog = np.array([r["dog"] for r in rows])
    coyote = np.array([r["coyote"] for r in rows])
    ax.bar(x, wolf, color=WOLF_COLOR, width=0.8, label="Gray wolf")
    ax.bar(x, dog, bottom=wolf, color=DOG_COLOR, width=0.8, label="Dog")
    ax.bar(x, coyote, bottom=wolf + dog, color=COYOTE_COLOR, width=0.8, label="Coyote")
    ax.errorbar(
        x,
        wolf + dog,
        yerr=[1.96 * r["dog_se"] for r in rows],
        fmt="none",
        ecolor="#1d1d1d",
        elinewidth=1,
        capsize=2,
    )
    for xi, r in zip(x, rows, strict=True):
        ax.text(xi, 1.02, f"{r['dog'] * 100:+.0f}%", ha="center", fontsize=6.5, color=DOG_COLOR)
    ax.set_xticks(
        x,
        [label(r["target"]).replace(" (Algonquin only)", "\n(Algonquin)") for r in rows],
        rotation=50,
        ha="right",
        fontsize=7.5,
    )
    for tick, r in zip(ax.get_xticklabels(), rows, strict=True):
        tick.set_color(color(r["target"]))
    ax.axvline(1.5, color="#999999", lw=0.8, ls=":")  # controls | targets
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("qpAdm mixture weight")
    ax.set_title("a  Three-way fit: gray wolf + dog + coyote (dog weight ±95% CI above bars)")
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.9)

    # (b) Coyote-source p-values.
    ax = fig.add_subplot(grid[0, 1])
    tlist = list(dict.fromkeys(r["target"] for r in source_rows))
    clist = list(COYOTES)
    pmat = np.array(
        [
            [
                next(
                    r["p_value"]
                    for r in source_rows
                    if r["target"] == t and r["coyote_source"] == c
                )
                for c in clist
            ]
            for t in tlist
        ]
    )
    image = ax.imshow(
        np.log10(np.maximum(pmat, 1e-6)), cmap="RdYlGn", vmin=-4, vmax=0, aspect="auto"
    )
    for i in range(len(tlist)):
        for j in range(len(clist)):
            p = pmat[i, j]
            ax.text(
                j, i, f"{p:.2g}" if p >= 1e-4 else "<1e-4", ha="center", va="center", fontsize=7.5
            )
    ax.set_xticks(range(len(clist)), [f"{COYOTES[c]}\n({c})" for c in clist], fontsize=8)
    ax.set_yticks(
        range(len(tlist)),
        [label(t).replace(" (Algonquin only)", " (Algonquin)") for t in tlist],
        fontsize=7.5,
    )
    for tick, t in zip(ax.get_yticklabels(), tlist, strict=True):
        tick.set_color(color(t))
    ax.axhline(1.5, color="white", lw=3)  # controls | targets
    ax.spines[:].set_visible(False)
    ax.set_xlabel("Coyote genome used as the coyote source (other two become outgroups)")
    ax.set_title("b  Which coyote fits? qpAdm p-value (wolf + that coyote)")
    cbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("log10 p  (green = model fits)", fontsize=7.5)

    # (c) D(coyote i, coyote j; X, jackal).
    ax = fig.add_subplot(grid[1, 0])
    tests = list(dict.fromkeys(r["test"] for r in d_rows))
    d_targets = list(dict.fromkeys(r["target"] for r in d_rows))
    offsets = np.linspace(-0.3, 0.3, len(d_targets))
    for ti, test in enumerate(tests):
        for off, target in zip(offsets, d_targets, strict=True):
            r = next(r for r in d_rows if r["test"] == test and r["target"] == target)
            ax.errorbar(
                ti + off,
                r["D"],
                yerr=1.96 * r["se"],
                fmt="o",
                ms=5,
                capsize=2,
                color=color(target),
                label=label(target) if ti == 0 else None,
            )
    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(
        range(len(tests)), [t.replace("; X, jackal)", ";\nX, jackal)") for t in tests], fontsize=8
    )
    ax.set_ylabel("D  (>0: X closer to the first coyote)")
    ax.set_title("c  Is X's coyote ancestry closer to one living coyote?")
    ax.legend(fontsize=7, frameon=False, ncols=2, loc="upper left")

    # (d) Reference rotation.
    ax = fig.add_subplot(grid[1, 1])
    rtargets = list(dict.fromkeys(r["target"] for r in rotation_rows))
    rng = np.random.default_rng(0)
    for yi, target in enumerate(rtargets[::-1]):
        rows_t = [r for r in rotation_rows if r["target"] == target]
        vals = np.array([r["wolf_ancestry"] for r in rows_t])
        jitter = rng.uniform(-0.18, 0.18, len(vals))
        ax.scatter(vals, yi + jitter, s=12, color=color(target), alpha=0.45, lw=0)
        base = next(r for r in rows_t if r["baseline"])
        ax.scatter(
            base["wolf_ancestry"],
            yi,
            marker="D",
            s=60,
            color=color(target),
            edgecolor="#1d1d1d",
            lw=0.8,
            zorder=3,
        )
        ax.text(
            vals.max() + 0.03,
            yi,
            f"{np.median(vals) * 100:.0f}%  ({vals.min() * 100:.0f} to {vals.max() * 100:.0f}%)",
            va="center",
            fontsize=7.5,
        )
    ax.set_yticks(
        range(len(rtargets)),
        [label(t).replace(" (Algonquin only)", "\n(Algonquin)") for t in rtargets[::-1]],
        fontsize=8,
    )
    ax.set_xlim(-0.15, 1.25)
    ax.set_xlabel("Gray-wolf ancestry (f4-ratio)")
    n_models = len([r for r in rotation_rows if r["target"] == rtargets[0]])
    ax.set_title(f"d  Every reference choice: {n_models} f4-ratio models per target")
    ax.text(
        0.99,
        0.02,
        "◆ baseline   ● alternative A/B/C/outgroup\nmedian (range)",
        transform=ax.transAxes,
        ha="right",
        fontsize=7,
        color="#555555",
    )

    fig.suptitle(
        "Robustness of the wolf ancestry cline",
        fontsize=14,
        fontweight="bold",
        x=0.01,
        ha="left",
        y=0.99,
    )
    fig.text(
        0.01,
        0.955,
        "Genotype-likelihood allele frequencies · qpAdm rights: golden jackal "
        "(base), Eurasian wolves, Asian wolves, village dogs, dhole · chromosome block "
        "jackknife",
        fontsize=8.5,
        color="#555555",
    )
    fig.savefig(ROBUST / "wolf_cline_robustness.png", dpi=160, bbox_inches="tight")
    fig.savefig(ROBUST / "wolf_cline_robustness.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
