"""Domain model for CANIS.

These types are the shared vocabulary every module speaks. They are deliberately plain
:mod:`dataclasses` (not tool-specific objects) so that a stage depends on *concepts*
(a :class:`Cohort`, a :class:`Callset`, an :class:`Artifact`) rather than on any other
module. Immutable records are frozen; working sets that are assembled incrementally are not.

Nothing here performs I/O or shells out — this module is pure data and stays importable in
any environment.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import StrEnum
from pathlib import Path
from typing import Any

# --------------------------------------------------------------------------------------
# Controlled vocabularies
# --------------------------------------------------------------------------------------


class CanidTaxon(StrEnum):
    """Canonical taxon labels. Free-form subspecies live on :class:`Taxon.subspecies`;
    this enum fixes the coarse label used for grouping, colouring, and QC checks."""

    GRAY_WOLF = "gray_wolf"  # Canis lupus
    DOMESTIC_DOG = "domestic_dog"  # Canis lupus familiaris
    VILLAGE_DOG = "village_dog"  # free-breeding domestic dog
    DINGO = "dingo"  # Canis (lupus) dingo
    COYOTE = "coyote"  # Canis latrans
    RED_WOLF = "red_wolf"  # Canis rufus
    EASTERN_WOLF = "eastern_wolf"  # Canis lycaon / lupus lycaon
    GOLDEN_JACKAL = "golden_jackal"  # Canis aureus
    AFRICAN_GOLDEN_WOLF = "african_golden_wolf"  # Canis lupaster
    ETHIOPIAN_WOLF = "ethiopian_wolf"  # Canis simensis
    BLACK_BACKED_JACKAL = "black_backed_jackal"  # Lupulella mesomelas
    SIDE_STRIPED_JACKAL = "side_striped_jackal"  # Lupulella adusta
    DHOLE = "dhole"  # Cuon alpinus
    AFRICAN_WILD_DOG = "african_wild_dog"  # Lycaon pictus (outgroup)
    OTHER = "other"

    @property
    def is_ingroup_canis(self) -> bool:
        """True for members of genus *Canis* commonly used as ingroup in canid studies."""
        return self in {
            CanidTaxon.GRAY_WOLF,
            CanidTaxon.DOMESTIC_DOG,
            CanidTaxon.VILLAGE_DOG,
            CanidTaxon.DINGO,
            CanidTaxon.COYOTE,
            CanidTaxon.RED_WOLF,
            CanidTaxon.EASTERN_WOLF,
            CanidTaxon.GOLDEN_JACKAL,
            CanidTaxon.AFRICAN_GOLDEN_WOLF,
            CanidTaxon.ETHIOPIAN_WOLF,
        }


class Sex(StrEnum):
    MALE = "male"
    FEMALE = "female"
    UNKNOWN = "unknown"


class SequencingPlatform(StrEnum):
    ILLUMINA = "illumina"
    BGI = "bgi"
    PACBIO = "pacbio"
    NANOPORE = "nanopore"
    UNKNOWN = "unknown"


class LibraryLayout(StrEnum):
    PAIRED = "paired"
    SINGLE = "single"


class ArtifactKind(StrEnum):
    """The typed 'ports' modules connect through. New kinds may be appended freely;
    existing values must never be renumbered/renamed (they persist in run manifests)."""

    RAW_READS = "raw_reads"  # FASTQ (possibly gzipped)
    ALIGNMENT = "alignment"  # CRAM/BAM
    GVCF = "gvcf"  # per-sample gVCF
    JOINT_VCF = "joint_vcf"  # multi-sample joint-genotyped VCF/BCF
    CALLSET = "callset"  # harmonized, analysis-ready variant set (VCF/BCF)
    SAMPLE_SHEET = "sample_sheet"  # validated per-sample metadata table
    PLINK = "plink"  # PLINK1/2 fileset (bed/bim/fam or pgen)
    ZARR = "zarr"  # chunked genotype array (sgkit/xarray)
    GENOTYPES = "genotypes"  # materialized genotype matrix (npz/zarr) for analysis
    QC_TABLE = "qc_table"  # per-sample/site QC metrics
    ANALYSIS_RESULT = "analysis_result"  # any tabular/figure result envelope
    TREE = "tree"  # Newick phylogeny
    REPORT = "report"  # rendered HTML/PDF/MD report
    REFERENCE = "reference"  # reference genome FASTA + indices


class FileFormat(StrEnum):
    FASTQ = "fastq"
    CRAM = "cram"
    BAM = "bam"
    VCF = "vcf"
    BCF = "bcf"
    PLINK_BED = "plink_bed"
    PLINK_PGEN = "plink_pgen"
    ZARR = "zarr"
    FASTA = "fasta"
    NEWICK = "newick"
    CSV = "csv"
    TSV = "tsv"
    PARQUET = "parquet"
    NPZ = "npz"
    JSON = "json"
    PNG = "png"
    SVG = "svg"
    HTML = "html"
    PDF = "pdf"
    MARKDOWN = "markdown"
    OTHER = "other"


# --------------------------------------------------------------------------------------
# Biological / geographic entities
# --------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class GeoLocality:
    """Sampling locality. Coordinates are WGS84 decimal degrees."""

    latitude: float | None = None
    longitude: float | None = None
    region: str | None = None
    country: str | None = None

    def __post_init__(self) -> None:
        if self.latitude is not None and not -90.0 <= self.latitude <= 90.0:
            raise ValueError(f"latitude out of range: {self.latitude}")
        if self.longitude is not None and not -180.0 <= self.longitude <= 180.0:
            raise ValueError(f"longitude out of range: {self.longitude}")

    @property
    def has_coordinates(self) -> bool:
        return self.latitude is not None and self.longitude is not None


@dataclass(frozen=True, slots=True)
class Taxon:
    """A taxonomic assignment: a coarse controlled label plus optional fine detail."""

    label: CanidTaxon
    scientific_name: str | None = None
    subspecies: str | None = None  # e.g. "arctos", "occidentalis", "familiaris"

    def __str__(self) -> str:
        parts = [self.scientific_name or self.label.value]
        if self.subspecies:
            parts.append(f"({self.subspecies})")
        return " ".join(parts)


@dataclass(frozen=True, slots=True)
class Individual:
    """A biological animal. The unit of most population-genetic analyses."""

    id: str
    taxon: Taxon
    population: str | None = None
    sex: Sex = Sex.UNKNOWN
    locality: GeoLocality = field(default_factory=GeoLocality)
    source_study: str | None = None  # DOI or short citation key
    dataset_id: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class Sample:
    """A sequencing library/run belonging to an :class:`Individual`.

    One individual may have several samples (re-sequencing, multiple libraries); the
    ``individual_id`` links them.
    """

    id: str
    individual_id: str
    accession: str | None = None  # e.g. SRR/ERR/DRR
    platform: SequencingPlatform = SequencingPlatform.UNKNOWN
    layout: LibraryLayout = LibraryLayout.PAIRED
    reported_coverage: float | None = None
    observed_coverage: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class GenomicInterval:
    """A half-open genomic interval ``[start, end)`` used as the unit of scatter/gather."""

    contig: str
    start: int = 0
    end: int | None = None  # None => to end of contig

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError("interval start must be >= 0")
        if self.end is not None and self.end <= self.start:
            raise ValueError("interval end must be > start")

    @property
    def length(self) -> int | None:
        return None if self.end is None else self.end - self.start

    def to_region_string(self) -> str:
        """Render as an htslib/samtools region, e.g. ``chr1:1-1000000`` (1-based, inclusive)."""
        if self.end is None:
            return f"{self.contig}:{self.start + 1}-" if self.start else self.contig
        return f"{self.contig}:{self.start + 1}-{self.end}"


@dataclass(frozen=True, slots=True)
class ReferenceGenome:
    """A reference assembly and its coordinate context."""

    id: str  # e.g. "UU_Cfam_GSD_1.0"
    assembly: str
    fasta: Path
    species: str = "Canis lupus familiaris"
    chrom_style: str = "ncbi"  # "ucsc" (chr1) | "ncbi" (1/NC_...) etc.
    autosomes: tuple[str, ...] = ()
    liftover_chains: dict[str, Path] = field(default_factory=dict)  # target_id -> chain

    @property
    def fai(self) -> Path:
        return self.fasta.with_suffix(self.fasta.suffix + ".fai")


@dataclass(frozen=True, slots=True)
class Dataset:
    """An immutable, versioned view of a published data source."""

    id: str
    version: str
    doi: str | None = None
    fetcher: str = "ena"  # name of a registered DatasetFetcher
    description: str = ""
    sample_ids: tuple[str, ...] = ()
    checksums: dict[str, str] = field(default_factory=dict)  # relpath -> sha256/md5

    @property
    def versioned_id(self) -> str:
        return f"{self.id}@{self.version}"


# --------------------------------------------------------------------------------------
# Working sets and artifacts
# --------------------------------------------------------------------------------------


@dataclass
class Cohort:
    """The working set for a run: the individuals/populations to analyse, against a
    reference. Assembled incrementally by acquisition/QC, so it is mutable."""

    id: str
    individuals: list[Individual] = field(default_factory=list)
    reference: ReferenceGenome | None = None
    filters: dict[str, Any] = field(default_factory=dict)

    @property
    def size(self) -> int:
        return len(self.individuals)

    @property
    def populations(self) -> dict[str, list[Individual]]:
        pops: dict[str, list[Individual]] = {}
        for ind in self.individuals:
            pops.setdefault(ind.population or "unassigned", []).append(ind)
        return pops

    @property
    def taxa(self) -> dict[CanidTaxon, int]:
        counts: dict[CanidTaxon, int] = {}
        for ind in self.individuals:
            counts[ind.taxon.label] = counts.get(ind.taxon.label, 0) + 1
        return counts

    def subset(self, individual_ids: set[str], *, new_id: str | None = None) -> Cohort:
        keep = [i for i in self.individuals if i.id in individual_ids]
        return replace(self, id=new_id or f"{self.id}.subset", individuals=keep)


@dataclass(frozen=True, slots=True)
class Artifact:
    """A typed handle to a produced file/dir plus the metadata needed for provenance and
    dependency resolution. This is the *only* thing that crosses between modules."""

    kind: ArtifactKind
    role: str  # semantic role, e.g. "harmonized_callset"
    path: Path
    fmt: FileFormat = FileFormat.OTHER
    checksum: str | None = None  # sha256 of the file (or manifest of a dir)
    produced_by: str | None = None  # stage name
    provenance_id: str | None = None  # links to a ProvenanceRecord
    schema_version: int = 1
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def key(self) -> tuple[ArtifactKind, str]:
        """The (kind, role) pair used to satisfy stage input specs."""
        return (self.kind, self.role)

    def exists(self) -> bool:
        return self.path.exists()


@dataclass(frozen=True, slots=True)
class AnalysisResult:
    """A typed envelope for a numeric/analytic result (PCA, F_ST, D-stats, ...)."""

    name: str
    parameters: dict[str, Any] = field(default_factory=dict)
    table: Path | None = None  # tabular payload (CSV/Parquet)
    figures: tuple[Path, ...] = ()
    summary: dict[str, Any] = field(default_factory=dict)
    provenance_id: str | None = None


@dataclass(frozen=True, slots=True)
class Callset:
    """A harmonized, analysis-ready variant set with its coordinate context."""

    id: str
    reference_id: str
    sample_ids: tuple[str, ...]
    n_sites: int | None = None
    biallelic_snps_only: bool = True
    variants: Artifact | None = None  # BCF/VCF handle
    plink: Artifact | None = None  # PLINK fileset handle
    zarr: Artifact | None = None  # chunked array handle
    filter_lineage: tuple[str, ...] = ()  # ordered list of filters applied

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)
