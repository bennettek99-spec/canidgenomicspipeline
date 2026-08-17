# CANIS — laptop-safe canid genomics

[![CI](https://github.com/bennettek99-spec/canidgenomicspipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/bennettek99-spec/canidgenomicspipeline/actions/workflows/ci.yml)

CANIS is a configuration-driven comparative evolutionary-genomics pipeline for canids. It is designed for a normal laptop first: it processes chromosomes sequentially, puts hard guards around downloads and disk space, resumes only verified unchanged work, and produces a portable HTML analysis package.

> Comparative evolutionary genomics only. CANIS is exploratory software, not a veterinary diagnostic or a standalone conservation-management decision tool.

## What this version does

- Safe resume cache: every stage is fingerprinted from its inputs, stage-relevant configuration, and transitive CANIS code dependencies, so unrelated report/config edits do not invalidate expensive acquisition. Cache records are small JSON metadata; VCFs and genomic arrays are never copied for caching.
- Atomic outputs and provenance: stages write to a temporary workspace, validate outputs, then promote them atomically. Full-content artifact checksums, the resolved YAML configuration, interrupted-stage recovery, and stable promoted paths are retained with the run manifest.
- Enforced QC: sample and site failures are physically removed before genotype loading. `qc_exclusions.csv` records every excluded sample/site and exact reason.
- Correct chromosome handling: D/f-statistic uncertainty uses chromosomes or fixed physical blocks, outgroups must be explicit, D tests include FDR-adjusted values, and local-ancestry HMMs restart at every chromosome and run sequentially.
- Guarded acquisition: `reduced_panel` performs exact byte-range preflight, reuses checksummed ranges after interruption, bounds concurrent memory, reports live ETA, and can enforce GQ/DP hard-call filters. `acquire_reads` retains its streamed, checksummed, pair-preserving safeguards and 9 GB default ceiling.
- Analysis readiness: a dedicated gate applies autosome/optional LD filtering and records small-population, near-duplicate, ascertainment, and hard-call-quality limitations before scientific stages run.
- Scalable local arrays: genotype matrices switch to memory-mapped `.npy` artifacts above a configurable threshold instead of loading one monolithic NPZ into RAM.
- Lightweight VCF harmonization: common-reference verification, SNP normalization, REF/ALT and strand-orientation checks, sample concordance, missingness/batch diagnostics. Liftover is an explicit external preparation step, never an automatic laptop action.
- Honest diversity and demography: fixed-panel diversity/heterozygosity are labelled panel-relative. Full per-callable-site π, θ, Tajima’s D, and Ne require an explicit callable-site denominator.
- Local control: `canidae ui` opens a loopback-only browser interface for dataset selection, preset/resource estimates, optional analyses, start/pause/stop/resume, progress, plain-language errors, and output links.
- Self-contained reporting: `report.html` opens by double-clicking and includes executive summary, QC badges, excluded records, PCA, tree/ancestry/introgression where configured, limitations, source accessions, downloads, methods, and provenance links.

## Quick start (Windows PowerShell)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev,analysis]"

# Validate the install and inspect available stages.
canidae version
canidae stages
pytest -q

# Start the local-only browser interface.
canidae ui
```

For a reproducible simulated cohort:

```powershell
python tests\simulated\make_cohort.py data
canidae run -c configs\examples\popgen_mvp.yaml `
  --set stages.ingest.callset=data\cohort.vcf `
  --set stages.load_genotypes.min_maf=0.0
```

The report is written to `data/store/report/report.html`; double-click it after the run. The adjacent `data/store/` stage folders and `runs/<run-id>/manifest.json` form the companion reproducibility package.

## Hybrid-canid diagnostics (bridge loci)

NYC coydog pedigree check + breed scoring, and eastern coyote three-way
coyote/wolf/dog mixture, are first-class stages (not only `scripts/`):

```powershell
# Requires pre-built bridge panel + call tables under data/ (see docs/ROADMAP.md).
canidae run -c configs\examples\nyc_coydog_validation.yaml `
  -c configs\profiles\laptop.yaml

canidae run -c configs\examples\eastern_coyote_ancestry.yaml `
  -c configs\profiles\laptop.yaml
```

These recipes are **bridge-locus diagnostics** (~250 sites), not whole-genome
ancestry. Stages: `reference_mixture`, `breed_assign`, `multiway_admixture`.
The reports say so, and the `breed_assign` "no single breed supported" flag is
the expected outcome at this marker count.

The equivalent `scripts/` entry points remain as shortcuts for one-off data
preparation, but the YAML recipes above are the supported path: they run through
the same executor, cache, and provenance as every other pipeline.

The estimators are covered by a deterministic synthetic bridge panel with known
mixture fractions, so `pytest` verifies that they recover the truth rather than
merely running:

```powershell
pytest tests\unit\test_hybrid_analysis.py tests\integration\test_hybrid_pipeline.py
```

`tests/golden/` pins both the synthetic results and the study results. The
study-panel snapshots skip unless the prepared panels are present under `data/`;
re-bless either set with `CANIDAE_UPDATE_GOLDEN=1` after reviewing the diff.

## Citing the data behind a run

Each public preset names its data sources through a citation bundle under
`configs/citations/`, and the report renders them under **Sources**:

```yaml
stages:
  report:
    citations: [nhgri_722g_wgs, eastern_coyote_radseq]
```

A bundle lists each study, its accession or DOI, what the recipe uses it for,
and its access conditions. Identifiers this repository does not record are shown
as "not recorded" rather than guessed — fill them in before citing in a
manuscript.

## Laptop reduced-panel workflow

The integrated public Red Wolf / golden jackal example starts from a remote indexed VCF rather than downloading the full source VCF or raw reads:


```powershell
canidae run -c configs\examples\redwolf_jackal_reduced_panel.yaml `
  -c configs\profiles\laptop.yaml
```

It stops after calculating a real byte-range estimate when explicit confirmation is required. Review that estimate and set `stages.reduced_panel.confirm_large_transfer: true` only if it is acceptable. See [the reduced-panel guide](docs/REDUCED_PANEL_PIPELINE.md).

## Pipeline layout

```text
configs/                         safe defaults, laptop profile, reproducible recipes
docs/                            operating guide, architecture, reduced-panel details
src/canidae/
  cli.py                         `canidae run`, `canidae ui`, config/stage commands
  ui.py                          local loopback-only browser controller
  core/                          cache, atomic staging, provenance, resource manager
  stages/
    acquisition/                 ingest and remote indexed reduced-panel extraction
    qc/                          enforced sample/site filtering and exclusion ledger
    processing/                  harmonization and VCF utilities
    popgen/, introgression/,     PCA, FST, diversity, D/f statistics, local ancestry,
    local_ancestry/, hybrid/,    hybrid bridge-locus mixture/breed diagnostics,
    ...                          phylogenetics, demography, reporting

tests/                           unit, simulated, and integration validation
```

The shipped executor is the bounded local backend. Scheduler/container/GPU execution,
automatic liftover, and large raw-read calling are future extension points rather than
current default capabilities; see [the implementation boundary](docs/ARCHITECTURE.md#current-implementation-boundary).

## Data and statistical scope

Use a sample sheet with at least `sample_id`, `taxon`, and `population`. For full-genome demographic estimates, provide a defensible callable-site denominator with `stages.diversity.callable_sites` and/or `stages.demography.callable_sites`. Without it, the report deliberately presents selected-SNP panel-relative values only.

Keep large data outside Git. The repository ignores VCF/BCF/BAM/CRAM/Zarr files, run outputs, staging workspaces, and cache metadata by default.

## Documentation

- [Laptop operation, pause/resume, and cleanup](docs/LAPTOP_OPERATIONS.md)
- [Integrated reduced-panel acquisition](docs/REDUCED_PANEL_PIPELINE.md)
- [Laptop-safe paired/interleaved/single-end FASTQ acquisition](docs/LAPTOP_READ_ACQUISITION.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Legacy public Red Wolf AADR-style extraction notes](docs/REDWOLF_JACKAL_AADR.md)
- [Changelog](CHANGELOG.md)

## License

MIT. See [LICENSE](LICENSE).
