"""Sample-sheet parsing, validation, and materialization into the domain model.

The sample sheet is the human-editable entry point to a run: one row per sample, with the
metadata that population/geographic analyses need. It is validated strictly (fail early,
before any compute) and normalized into :class:`~canidae.core.model.Individual` objects and
a :class:`~canidae.core.model.Cohort`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from canidae.core.errors import StageInputError
from canidae.core.model import (
    CanidTaxon,
    Cohort,
    GeoLocality,
    Individual,
    Taxon,
)

REQUIRED_COLUMNS = ("sample_id", "taxon", "population")
OPTIONAL_COLUMNS = (
    "subspecies",
    "scientific_name",
    "region",
    "country",
    "latitude",
    "longitude",
    "source_study",
    "accession",
)


def read_sample_sheet(path: Path) -> pd.DataFrame:
    """Read and validate a sample-sheet CSV into a normalized DataFrame.

    Raises :class:`StageInputError` with an actionable message on any schema violation.
    """
    path = Path(path)
    if not path.exists():
        raise StageInputError(f"sample sheet not found: {path}")
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except Exception as exc:  # malformed CSV
        raise StageInputError(f"could not parse sample sheet {path}: {exc}") from exc

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise StageInputError(
            f"sample sheet {path} missing required column(s): {missing}. "
            f"Required: {list(REQUIRED_COLUMNS)}"
        )

    df["sample_id"] = df["sample_id"].str.strip()
    if (df["sample_id"] == "").any():
        raise StageInputError("sample sheet has empty sample_id value(s)")
    dupes = df["sample_id"][df["sample_id"].duplicated()].unique().tolist()
    if dupes:
        raise StageInputError(f"duplicate sample_id(s) in sample sheet: {dupes}")

    df["taxon"] = df["taxon"].str.strip().str.lower()
    _validate_taxa(df["taxon"], path)
    df["population"] = df["population"].str.strip()
    if (df["population"] == "").any():
        raise StageInputError("sample sheet has empty population value(s)")

    _validate_coordinates(df, path)

    # Ensure optional columns exist so downstream code can rely on them.
    for col in OPTIONAL_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    ordered = ["sample_id", "taxon", "population", *OPTIONAL_COLUMNS]
    return df[ordered].reset_index(drop=True)


def _validate_taxa(taxa: pd.Series, path: Path) -> None:
    valid = {t.value for t in CanidTaxon}
    unknown = sorted(set(taxa) - valid)
    if unknown:
        raise StageInputError(
            f"sample sheet {path} has unrecognized taxon label(s): {unknown}. "
            f"Use one of {sorted(valid)} (or 'other')."
        )


def _validate_coordinates(df: pd.DataFrame, path: Path) -> None:
    for axis, lo, hi in (("latitude", -90.0, 90.0), ("longitude", -180.0, 180.0)):
        if axis not in df.columns:
            continue
        for raw in df[axis]:
            if raw == "":
                continue
            try:
                val = float(raw)
            except ValueError:
                raise StageInputError(
                    f"sample sheet {path}: {axis} value '{raw}' is not numeric"
                ) from None
            if not lo <= val <= hi:
                raise StageInputError(
                    f"sample sheet {path}: {axis} {val} out of range [{lo}, {hi}]"
                )


def to_individuals(df: pd.DataFrame, *, dataset_id: str | None = None) -> list[Individual]:
    """Materialize validated rows into :class:`Individual` objects."""
    individuals: list[Individual] = []
    for row in df.itertuples(index=False):
        taxon = Taxon(
            label=CanidTaxon(row.taxon),
            scientific_name=row.scientific_name or None,
            subspecies=row.subspecies or None,
        )
        locality = GeoLocality(
            latitude=float(row.latitude) if row.latitude else None,
            longitude=float(row.longitude) if row.longitude else None,
            region=row.region or None,
            country=row.country or None,
        )
        individuals.append(
            Individual(
                id=row.sample_id,
                taxon=taxon,
                population=row.population,
                locality=locality,
                source_study=row.source_study or None,
                dataset_id=dataset_id,
            )
        )
    return individuals


def to_cohort(df: pd.DataFrame, *, cohort_id: str, dataset_id: str | None = None) -> Cohort:
    return Cohort(id=cohort_id, individuals=to_individuals(df, dataset_id=dataset_id))
