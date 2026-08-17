# ruff: noqa: E501
"""Render a portable HTML report for the public NYC coydog validation."""

from __future__ import annotations

import argparse
import csv
import html
import json
from datetime import UTC, datetime
from pathlib import Path

SAMPLES = ("T211", "NY04", "NY05", "NY01")


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def fmt_fraction(value: object | None) -> str:
    return "Not estimable" if value is None else f"{float(value) * 100:.1f}%"


def table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(value)}</td>" for value in row) + "</tr>"
        for row in rows
    )
    return f"<div class='scroll'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def bar_chart(
    strict: dict[str, dict[str, object]], sensitivity: dict[str, dict[str, object]]
) -> str:
    """Return a compact inline SVG comparing strict and sensitivity estimates."""
    width, height, left, bottom = 740, 330, 65, 52
    plot_height = height - bottom - 30
    group_width = (width - left - 35) / len(SAMPLES)
    svg = [
        f"<svg viewBox='0 0 {width} {height}' role='img' aria-label='Dog ancestry estimates by sample'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
    ]
    for tick in range(6):
        value = tick / 5
        y = 30 + plot_height * (1 - value)
        svg.append(f"<line x1='{left}' y1='{y:.1f}' x2='{width - 25}' y2='{y:.1f}' class='grid'/>")
        svg.append(
            f"<text x='{left - 10}' y='{y + 4:.1f}' text-anchor='end' class='axis'>{value:.0%}</text>"
        )
    for index, sample in enumerate(SAMPLES):
        x = left + index * group_width + group_width / 2
        svg.append(
            f"<text x='{x:.1f}' y='{height - 20}' text-anchor='middle' class='label'>{sample}</text>"
        )
        for offset, (name, color, dataset) in enumerate(
            (
                ("Strict depth >=8", "#2563eb", strict),
                ("Sensitivity depth >=4", "#f59e0b", sensitivity),
            )
        ):
            value = dataset[sample]["dog_fraction"]
            if value is None:
                svg.append(
                    f"<text x='{x + (-14 if offset == 0 else 14):.1f}' y='{height - 42}' class='missing'>NA</text>"
                )
                continue
            bar_height = plot_height * float(value)
            bar_x = x + (-21 if offset == 0 else 3)
            bar_y = 30 + plot_height - bar_height
            svg.append(
                f"<rect x='{bar_x:.1f}' y='{bar_y:.1f}' width='18' height='{bar_height:.1f}' "
                f"rx='3' fill='{color}'><title>{sample}: {name} {float(value):.1%}</title></rect>"
            )
    svg.extend(
        [
            "<rect x='75' y='9' width='12' height='12' rx='2' fill='#2563eb'/>",
            "<text x='93' y='20' class='legend'>Strict calls (depth >=8)</text>",
            "<rect x='260' y='9' width='12' height='12' rx='2' fill='#f59e0b'/>",
            "<text x='278' y='20' class='legend'>Sensitivity (depth >=4)</text>",
            "</svg>",
        ]
    )
    return "".join(svg)


def calibration_chart(rows: list[dict[str, str]]) -> str:
    width, height, left, bottom = 740, 250, 65, 44
    plot_height = height - bottom - 26
    groups = {"coyote": [], "dog": []}
    for row in rows:
        groups[row["reference_group"]].append(float(row["leave_one_out_dog_fraction"]))
    svg = [
        f"<svg viewBox='0 0 {width} {height}' role='img' aria-label='Leave-one-out reference calibration'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
    ]
    for tick in range(6):
        value = tick / 5
        y = 26 + plot_height * (1 - value)
        svg.append(f"<line x1='{left}' y1='{y:.1f}' x2='{width - 25}' y2='{y:.1f}' class='grid'/>")
        svg.append(
            f"<text x='{left - 10}' y='{y + 4:.1f}' text-anchor='end' class='axis'>{value:.0%}</text>"
        )
    for group_index, (group, color) in enumerate((("coyote", "#64748b"), ("dog", "#16a34a"))):
        center = 245 + group_index * 260
        for index, value in enumerate(groups[group]):
            x = center + (index - (len(groups[group]) - 1) / 2) * 15
            y = 26 + plot_height * (1 - value)
            svg.append(
                f"<circle cx='{x:.1f}' cy='{y:.1f}' r='5' fill='{color}' opacity='.85'><title>{group}: {value:.1%}</title></circle>"
            )
        svg.append(
            f"<text x='{center}' y='{height - 14}' text-anchor='middle' class='label'>{group.title()} references</text>"
        )
    svg.append("</svg>")
    return "".join(svg)


def breed_html(breed: dict[str, object]) -> str:
    """Render the dog-parent breed-assignment section from the manifest."""
    calibration = breed["calibration"]
    results = breed["results"]
    sources = breed["sources"]
    rows = []
    for sample in SAMPLES:
        entry = results.get(sample)
        if entry is None:
            rows.append([sample, "Not estimable (no dog-fraction estimate)", "-", "-", "-"])
            continue
        rows.append(
            [
                sample,
                fmt_fraction(entry["dog_fraction"]),
                str(entry["best_breed"]),
                f"{entry['single_breed_gap']:+.2f}",
                "No" if not entry["single_breed_supported"] else "Yes",
            ]
        )
    accuracy = calibration["top1_accuracy"]
    return f"""
<h2>Dog-parent breed assignment</h2>
<p>The 722-genome NHGRI panel (~144 breeds) shares the same VCF records, so all {sources['n_wgs_samples']} samples were recovered at the 252 bridge loci with one repeat of the indexed byte-range transfer. {sources['n_candidate_groups']} candidate groups were compared against the NYC animals' dog component.</p>
<section class="callout warn"><strong>Result: no single breed matches.</strong> For every NYC animal the pooled any-dog panel fits as well as or better than every named breed or village-dog population (all likelihood gaps within 0.3 log units of zero). The dog ancestry is most consistent with a mixed-breed or unregistered dog parent, and/or a parent breed below the resolution of 252 loci.</section>
<section class="callout"><strong>Calibration: {accuracy:.0%} leave-one-out top-1 accuracy ({calibration['n_tested']} animals, {sources['n_candidate_groups']} groups).</strong> Large distinctive breeds assign reliably (e.g. Bernese Mountain Dog 18/18) but most small panels (3-8 members) do not resolve at this marker count. Single-breed calls are therefore reported as rankings only, never as identification.</section>
{table(["Sample", "Dog fraction", "Best-fitting group", "Gap vs any-dog (LL)", "Single breed supported"], rows)}
<p class="small">Depth-4 sensitivity (including NY04, which lacks a strict estimate) gives the same conclusion: best-fitting groups are village-dog pools with gaps near zero.</p>"""


def render(
    strict: dict[str, object],
    sensitivity: dict[str, object],
    calibration: list[dict[str, str]],
    breed: dict[str, object] | None = None,
) -> str:
    strict_rows = {row["sample_id"]: row for row in strict["results"]}  # type: ignore[index]
    low_rows = {row["sample_id"]: row for row in sensitivity["results"]}  # type: ignore[index]
    result_rows = []
    for sample in SAMPLES:
        high = strict_rows[sample]
        low = low_rows[sample]
        result_rows.append(
            [
                sample,
                fmt_fraction(high["dog_fraction"]),
                str(high["diagnostic_sites_called"]),
                fmt_fraction(low["dog_fraction"]),
                str(low["diagnostic_sites_called"]),
            ]
        )
    calibration_rows = [
        [
            row["sample_id"],
            row["reference_group"],
            fmt_fraction(row["leave_one_out_dog_fraction"]),
            row["diagnostic_sites"],
        ]
        for row in calibration
    ]
    generated = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    breed_section = breed_html(breed) if breed else ""
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>NYC Coydog Validation Report</title><style>
:root {{ color-scheme: light; --ink:#172033; --muted:#64748b; --line:#dbe3ee; --blue:#2563eb; --amber:#b45309; --red:#b42318; --green:#15803d; }}
body {{ font-family: Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif; margin:0; color:var(--ink); background:#f8fafc; }}
main {{ max-width:1040px; margin:auto; padding:34px 22px 56px; }} h1 {{ margin:0 0 7px; font-size:clamp(2rem,5vw,3.15rem); letter-spacing:-.04em; }} h2 {{ margin:38px 0 12px; font-size:1.35rem; }} h3 {{ margin:20px 0 8px; font-size:1rem; }}
.eyebrow {{ color:var(--blue); font-size:.78rem; text-transform:uppercase; letter-spacing:.12em; font-weight:750; }} .meta,.small {{ color:var(--muted); }} .meta {{ margin:0; }}
.grid2 {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:15px; margin:23px 0; }} .card,.callout {{ background:#fff; border:1px solid var(--line); border-radius:14px; padding:18px; }} .card strong {{ display:block; font-size:1.65rem; margin-bottom:3px; }}
.callout {{ border-left:5px solid var(--red); background:#fff8f7; }} .callout strong {{ color:var(--red); }} .ok {{ border-left-color:var(--green); background:#f4fbf6; }} .warn {{ border-left-color:var(--amber); background:#fffbeb; }}
table {{ width:100%; border-collapse:collapse; font-size:.91rem; }} th,td {{ padding:10px 9px; border-bottom:1px solid var(--line); text-align:left; white-space:nowrap; }} th {{ font-size:.77rem; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); }} .scroll {{ overflow:auto; background:#fff; border:1px solid var(--line); border-radius:12px; padding:0 11px; }}
svg {{ display:block; width:100%; height:auto; background:#fff; border:1px solid var(--line); border-radius:12px; }} .grid {{ stroke:#e2e8f0; }} .axis,.label,.legend,.missing {{ fill:#64748b; font:12px system-ui,sans-serif; }} .label {{ fill:#334155; font-weight:650; }} .missing {{ fill:#b42318; font-weight:700; }}
code {{ background:#eff6ff; padding:2px 5px; border-radius:4px; }} a {{ color:var(--blue); }} li {{ margin:8px 0; }} footer {{ color:var(--muted); border-top:1px solid var(--line); margin-top:38px; padding-top:16px; font-size:.85rem; }}
@media(max-width:640px) {{ .grid2 {{ grid-template-columns:1fr; }} main {{ padding:24px 14px 40px; }} }}
</style></head><body><main>
<div class="eyebrow">CANIS · public-data reanalysis</div><h1>NYC coydog validation</h1>
<p class="meta">Generated {generated} · PRJNA857904 RADseq samples cross-checked against PRJNA448733 WGS references</p>
<section class="callout"><strong>Result: not validated as a standalone quantitative ancestry estimate.</strong><br>
NY01 has an independent dog-ancestry signal near the published F1 value, but the strict panel cannot estimate NY04 and the coyote mother exceeds the prespecified calibration threshold. The data support a qualitative hybrid gradient, not a publishable new ancestry proportion.</section>
<div class="grid2"><section class="card"><strong>40.2%</strong><span>NY01 dog fraction, strict calls (published estimate: 46%)</span></section><section class="card"><strong>252</strong><span>Shared bridge loci retained across public WGS and RADseq data</span></section></div>
<h2>Observed ancestry signal</h2><p>The blue bars are the primary analysis: minimum depth 8, base quality 20, mapping quality 30. Orange bars lower only the depth threshold to 4 and are a sensitivity check, not the reported estimate.</p>
{bar_chart(strict_rows, low_rows)}
{table(["Sample", "Strict dog fraction", "Strict diagnostic calls", "Depth-4 sensitivity", "Sensitivity calls"], result_rows)}
<section class="callout warn"><strong>How to read this.</strong> Strict calls reproduce the expected ordering where data are sufficient (T211 &lt; NY05 &lt; NY01), but NY04 has only 7 strict diagnostic calls. At depth 4, all samples shift upward—including T211—so those values are not treated as calibrated ancestry estimates.</section>
<h2>Reference-panel calibration</h2><p>Leave-one-out scoring cleanly separates the WGS reference animals: coyote controls average 8.7% and dog controls average 95.1% dog fraction. That demonstrates the reference panel is discriminative; the weak point is sparse cross-platform observation in the NYC RADseq samples.</p>
{calibration_chart(calibration)}
{table(["Reference sample", "Group", "Leave-one-out dog fraction", "Diagnostic loci"], calibration_rows)}
{breed_section}
<h2>Methods and sources</h2><ul><li>Reference panel: 3 coyotes and 10 dogs drawn from public WGS BioProject <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA448733">PRJNA448733</a>.</li><li>Query samples: T211, NY04, NY05, and NY01 from aligned RADseq records in <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA857904">PRJNA857904</a>.</li><li>Only 267 pre-existing WGS/RADseq bridge loci were considered; 252 had complete reference calls. A least-squares dog-fraction score used loci with dog-coyote frequency separation of at least 0.40.</li><li>WGS extraction used 337 indexed byte ranges (670.2 MB actual) under a 9 GB transfer ceiling. The 149 MB NYC archive size was included in the conservative preflight.</li><li>Breed assignment: the same 337 byte ranges were re-fetched once to recover all 722 NHGRI panel samples at the 252 bridge loci (667.9 MB); 44 candidate groups (39 breeds plus pooled village-dog populations) were scored with Laplace-smoothed allele-frequency panels under a mixture likelihood.</li></ul>
<h2>Limitations</h2><ul><li>This is a sparse cross-platform bridge panel, not a jointly called genome-wide VCF.</li><li>The dog reference set is broad rather than specifically sampled from NYC dogs.</li><li>The result validates a qualitative pattern only; it must not be used as a veterinary, legal, conservation, or management decision.</li></ul>
<h2>Companion outputs</h2><ul><li><a href="validation.json">Strict validation manifest</a></li><li><a href="ancestry_validation.csv">Strict ancestry table</a></li><li><a href="reference_leave_one_out.csv">Reference calibration table</a></li><li><a href="../nyc_coydog_validation_depth4/validation.json">Depth-4 sensitivity manifest</a></li><li><a href="../nyc_coydog_validation_depth4/ancestry_validation.csv">Depth-4 sensitivity table</a></li><li><a href="../nyc_coydog_breeds/breed_assignment.json">Breed assignment manifest</a></li><li><a href="../nyc_coydog_breeds/breed_scores.csv">Breed likelihood rankings</a></li><li><a href="../nyc_coydog_breeds/leave_one_out.csv">Breed-panel leave-one-out calibration</a></li></ul>
<footer>Self-contained HTML report. The conclusion intentionally remains conservative: it corroborates the published NYC coydog signal but does not produce a replacement whole-genome ancestry estimate.</footer>
</main></body></html>"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--strict", type=Path, default=Path("data/nyc_coydog_validation/validation.json")
    )
    parser.add_argument(
        "--sensitivity",
        type=Path,
        default=Path("data/nyc_coydog_validation_depth4/validation.json"),
    )
    parser.add_argument(
        "--calibration",
        type=Path,
        default=Path("data/nyc_coydog_validation/reference_leave_one_out.csv"),
    )
    parser.add_argument("--out", type=Path, default=Path("data/nyc_coydog_validation/report.html"))
    parser.add_argument(
        "--breed-json",
        type=Path,
        default=Path("data/nyc_coydog_breeds/breed_assignment.json"),
        help="omit to skip the breed-assignment section",
    )
    args = parser.parse_args()
    strict, sensitivity = read_json(args.strict), read_json(args.sensitivity)
    breed = read_json(args.breed_json) if args.breed_json and args.breed_json.exists() else None
    with args.calibration.open(newline="", encoding="utf-8") as handle:
        calibration = list(csv.DictReader(handle))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render(strict, sensitivity, calibration, breed), encoding="utf-8")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
