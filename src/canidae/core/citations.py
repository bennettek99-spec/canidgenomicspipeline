"""Source-citation bundles for public presets.

A preset points the report stage at one or more YAML bundles under
``configs/citations/``. Each bundle names the published studies and archive
accessions a recipe draws on, so a rendered report always carries the
attribution and access conditions for the data behind it, whether or not the
run happens to produce a callset artifact.

Fields are deliberately permissive: ``doi`` is optional because this repository
records accessions for every source but not always a DOI, and an empty string is
an honest "not recorded here" rather than a fabricated identifier.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, Field, field_validator

from canidae.core.errors import ConfigError

_ACCESSION_URLS = {
    "bioproject": "https://www.ncbi.nlm.nih.gov/bioproject/{accession}",
    "sra_run": "https://www.ncbi.nlm.nih.gov/sra/{accession}",
}


class SourceCitation(BaseModel):
    """One published data source used by a preset."""

    key: str = Field(..., description="Short stable identifier, unique within a bundle.")
    label: str = Field(..., description="Human-readable name of the dataset or study.")
    study: str = ""
    year: int | None = None
    accession: str = ""
    accession_kind: str = Field(
        default="", description="bioproject | sra_run | other; drives the archive link."
    )
    doi: str = ""
    url: str = ""
    used_for: str = Field(default="", description="What this recipe uses the source for.")
    access: str = Field(default="", description="Licence or access conditions.")

    @field_validator("doi")
    @classmethod
    def _strip_doi_prefix(cls, value: str) -> str:
        """Accept a bare DOI or a doi.org URL; store the bare form."""
        value = value.strip()
        for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
            if value.lower().startswith(prefix):
                return value[len(prefix) :]
        return value

    @property
    def resolved_url(self) -> str:
        """An explicit url, else a DOI link, else the archive link for the accession."""
        if self.url:
            return self.url
        if self.doi:
            return f"https://doi.org/{self.doi}"
        template = _ACCESSION_URLS.get(self.accession_kind)
        if template and self.accession:
            return template.format(accession=self.accession)
        return ""


class CitationBundle(BaseModel):
    """A named group of sources, normally one per upstream dataset."""

    id: str
    label: str
    description: str = ""
    sources: list[SourceCitation] = Field(default_factory=list)


def load_citation_bundle(path: Path) -> CitationBundle:
    """Read and validate a single bundle YAML."""
    try:
        payload = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    except OSError as exc:
        raise ConfigError(f"citation bundle not readable: {path} ({exc})") from exc
    if not isinstance(payload, dict):
        raise ConfigError(f"citation bundle must be a mapping: {path}")
    try:
        return CitationBundle.model_validate(payload)
    except Exception as exc:  # pydantic ValidationError
        raise ConfigError(f"invalid citation bundle {path}: {exc}") from exc


def _search_roots(root: Path) -> list[Path]:
    """Where a relative bundle reference may live, most specific first.

    ``paths.root`` is a run's data root and is often a scratch directory, so the
    shipped ``configs/citations/`` bundles are also looked up from the working
    directory and from the checkout that provides this package.
    """
    checkout = Path(__file__).resolve().parents[3]
    ordered = [Path(root).resolve(), Path.cwd().resolve(), checkout]
    return list(dict.fromkeys(ordered))


def resolve_bundle_path(entry: Path, root: Path) -> Path:
    """Resolve one bundle reference: a bare id, a relative path, or an absolute path."""
    candidate = Path(entry)
    if candidate.is_absolute():
        if candidate.exists():
            return candidate
        raise ConfigError(f"citation bundle not found: {candidate}")

    relative = (
        Path("configs") / "citations" / f"{candidate.name}.yaml"
        if not candidate.suffix
        else candidate
    )
    tried = []
    for base in _search_roots(root):
        resolved = base / relative
        tried.append(resolved)
        if resolved.exists():
            return resolved
    raise ConfigError(f"citation bundle not found: {entry} (looked in {[str(p) for p in tried]})")


def load_citation_bundles(paths: list[Path], root: Path) -> list[CitationBundle]:
    """Load several bundles, resolving each reference via :func:`resolve_bundle_path`.

    A bare name such as ``nhgri_722g_wgs`` resolves to
    ``configs/citations/nhgri_722g_wgs.yaml`` so presets can reference bundles
    without repeating the directory. Repeated bundles are listed once.
    """
    bundles: list[CitationBundle] = []
    seen: set[str] = set()
    for entry in paths:
        bundle = load_citation_bundle(resolve_bundle_path(Path(entry), root))
        if bundle.id in seen:
            continue  # two presets may share a bundle; list it once
        seen.add(bundle.id)
        bundles.append(bundle)
    return bundles


def citation_rows(bundles: list[CitationBundle]) -> list[dict[str, str]]:
    """Flatten bundles into report-table rows."""
    rows: list[dict[str, str]] = []
    for bundle in bundles:
        for source in bundle.sources:
            rows.append(
                {
                    "dataset": bundle.label,
                    "source": source.label,
                    "study": " ".join(
                        part for part in (source.study, str(source.year or "")) if part
                    ),
                    "accession": source.accession or "not recorded",
                    "doi": source.doi or "not recorded",
                    "link": source.resolved_url,
                    "used_for": source.used_for,
                    "access": source.access,
                }
            )
    return rows
