"""Combine the wolf-cline figures and ten-group results into one PNG sheet."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "canidae-matplotlib"))
from matplotlib.font_manager import FontProperties, findfont

ROOT = Path(__file__).resolve().parents[1]
CLINE = ROOT / "results" / "wolf_ancestry_cline"
STUDY = ROOT / "results" / "north_american_wolf_admixture"
OUTPUT = STUDY / "wolf_study_overview.png"

WIDTH, HEIGHT = 4400, 3670
INK = "#18304a"
MUTED = "#53677b"
LIGHT = "#eef3f7"
TEAL = "#209c91"
PURPLE = "#9271b8"
RED = "#ce604d"
GOLD = "#d6a03b"


def main() -> None:
    def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
        path = findfont(FontProperties(family="DejaVu Sans", weight="bold" if bold else "normal"))
        return ImageFont.truetype(path, size)

    im = Image.new("RGB", (WIDTH, HEIGHT), "#f4f7fa")
    draw = ImageDraw.Draw(im)
    draw.rectangle((0, 0, WIDTH, 215), fill=INK)
    draw.text((82, 39), "NORTH AMERICAN WOLF ANCESTRY", fill="white", font=font(65, bold=True))
    draw.text(
        (84, 130),
        "Cline, robustness checks, paired contrasts, and the ten-group evidence audit",
        fill="#dce8f1",
        font=font(34),
    )

    def figure_panel(path: Path, box: tuple[int, int, int, int], label: str) -> None:
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=24, fill="white", outline="#d9e3eb", width=3)
        draw.text((x0 + 30, y0 + 21), label, fill=INK, font=font(34, bold=True))
        with Image.open(path) as original:
            img = ImageOps.contain(original.convert("RGB"), (x1 - x0 - 38, y1 - y0 - 105))
            im.paste(img, (x0 + (x1 - x0 - img.width) // 2, y0 + 82))

    figure_panel(
        CLINE / "wolf_ancestry_cline.png", (65, 245, 2168, 2120), "01  Genome-wide ancestry cline"
    )
    figure_panel(
        CLINE / "robustness" / "wolf_cline_robustness.png",
        (2220, 245, 4335, 2120),
        "02  Admixture and reference robustness",
    )
    figure_panel(
        CLINE / "robustness" / "paired_cline_contrasts.png",
        (65, 2150, 2055, 3510),
        "03  Paired cline differences",
    )

    box = (2105, 2150, 4335, 3510)
    x0, y0, x1, _y1 = box
    draw.rounded_rectangle(box, radius=24, fill="white", outline="#d9e3eb", width=3)
    draw.text(
        (x0 + 40, y0 + 25), "04  Ten-group ancestry readout", fill=INK, font=font(38, bold=True)
    )
    draw.text(
        (x0 + 40, y0 + 85),
        "Local 87,538-site panel • 38 chromosome jackknife blocks",
        fill=MUTED,
        font=font(26),
    )

    with (STUDY / "local_contrasts.csv").open(newline="", encoding="utf-8") as handle:
        rows = {row["group"].split(" (")[0]: row for row in csv.DictReader(handle)}

    labels = [
        ("Great Lakes wolf", "Great Lakes", TEAL),
        ("Eastern/Algonquin wolf", "Eastern / Algonquin", PURPLE),
        ("Red wolf", "Red wolf", RED),
        ("Mexican wolf", "Mexican", GOLD),
        ("Yellowstone wolf", "Yellowstone", "#5b7795"),
        ("Alaskan wolf", "Alaska", "#5b7795"),
    ]
    yy = y0 + 158
    draw.text((x0 + 40, yy), "GROUP", fill=MUTED, font=font(23, bold=True))
    draw.text((x0 + 750, yy), "WOLF-SIDE f4 RATIO", fill=MUTED, font=font(23, bold=True))
    draw.text((x0 + 1530, yy), "COYOTE D  (Z)", fill=MUTED, font=font(23, bold=True))
    yy += 53
    for key, label, color in labels:
        row = rows[key]
        draw.line((x0 + 38, yy - 7, x1 - 38, yy - 7), fill="#e7edf2", width=2)
        draw.ellipse((x0 + 41, yy + 19, x0 + 60, yy + 38), fill=color)
        draw.text((x0 + 76, yy + 8), label, fill=INK, font=font(27, bold=True))
        ratio = float(row["wolf_f4_ratio"]) * 100
        se = float(row["wolf_f4_ratio_se"]) * 100
        draw.text((x0 + 750, yy + 8), f"{ratio:.1f}% ± {se:.1f}", fill=INK, font=font(27))
        d = float(row["coyote_d"])
        z = float(row["coyote_d_z"])
        draw.text(
            (x0 + 1530, yy + 8),
            f"{d:+.3f}  ({z:.2f})",
            fill=color if z > 3 else MUTED,
            font=font(27),
        )
        yy += 75

    draw.rounded_rectangle((x0 + 38, yy + 11, x1 - 38, yy + 116), radius=12, fill=LIGHT)
    draw.text(
        (x0 + 58, yy + 22),
        "Supported coyote-side model complements:",
        fill=INK,
        font=font(24, bold=True),
    )
    draw.text(
        (x0 + 58, yy + 64),
        "Great Lakes 25.4% ± 3.5  •  Algonquin 35.6% ± 2.7  •  Red 64.0% ± 3.4",
        fill=INK,
        font=font(24),
    )
    yy += 144
    draw.text(
        (x0 + 40, yy), "No matched local nuclear genotypes", fill=INK, font=font(28, bold=True)
    )
    yy += 47
    for label in (
        "Historical Plains (C. l. nubilus; C. variabilis Wied 1841)",
        "Pacific Coast",
        "Atlantic Coast",
        "West Arctic",
    ):
        draw.text((x0 + 57, yy), "• " + label, fill=MUTED, font=font(25))
        yy += 43

    yy += 18
    draw.line((x0 + 40, yy, x1 - 40, yy), fill="#d8e2ea", width=2)
    yy += 20
    notes = [
        "f4 ratios are model-dependent. Only the three eastern groups also have",
        "strong local coyote D evidence. Other ratios are not admixture fractions.",
        "Dog: the local panel gives no reliable genome-wide percentage.",
        "Published K-locus evidence exists for Yellowstone and Alaska.",
    ]
    for note in notes:
        draw.text((x0 + 40, yy), note, fill=MUTED, font=font(23))
        yy += 34

    draw.text(
        (80, 3560),
        "Source: results/north_american_wolf_admixture/local_contrasts.csv"
        " • Cline figures from results/wolf_ancestry_cline",
        fill=MUTED,
        font=font(26),
    )
    draw.text(
        (80, 3602),
        "Regional published context and limitations: doi:10.1371/journal.pgen.1007745"
        " • Full methods and accessions in the study README",
        fill=MUTED,
        font=font(25),
    )
    im.save(OUTPUT, format="PNG", optimize=True)
    print(OUTPUT)


if __name__ == "__main__":
    main()
