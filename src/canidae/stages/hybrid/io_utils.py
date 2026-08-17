"""Shared path / JSON helpers for hybrid diagnostic stages."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from canidae.core.errors import StageInputError


def resolve_path(path: Path, root: Path) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    return candidate.resolve()


def require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise StageInputError(f"{label} not found: {path}")
    return path


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_csv_dicts(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_locus_key(key: str) -> tuple[str, int]:
    chrom, pos = key.rsplit(":", 1)
    return chrom, int(pos)


def load_calls_csv(
    path: Path,
    sites: dict[tuple[str, int], tuple[str, str]],
    min_depth: int,
) -> dict[tuple[str, int], tuple[str | None, int, int, int]]:
    """Reapply a depth threshold to saved per-site allele counts."""
    calls: dict[tuple[str, int], tuple[str | None, int, int, int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["chrom"], int(row["position"]))
            if key not in sites:
                continue
            depth = int(row["depth"])
            ref_count = int(row["ref_count"])
            alt_count = int(row["alt_count"])
            if depth < min_depth:
                calls[key] = (None, depth, ref_count, alt_count)
                continue
            fraction = alt_count / depth
            call = "0/0" if fraction <= 0.15 else "1/1" if fraction >= 0.85 else "0/1"
            calls[key] = (call, depth, ref_count, alt_count)
    if set(calls) != set(sites):
        raise StageInputError(f"saved calls do not cover the requested sites: {path}")
    return calls


def load_wgs_panel_json(
    path: Path,
) -> tuple[list[str], dict[tuple[str, int], list[str | None]], list[tuple[str, int]]]:
    payload = load_json(path)
    samples = list(payload["samples"])
    records = {
        parse_locus_key(key): list(values) for key, values in payload["loci"].items()
    }
    keys = sorted(records, key=lambda item: (int(item[0].removeprefix("chr")), item[1]))
    return samples, records, keys
