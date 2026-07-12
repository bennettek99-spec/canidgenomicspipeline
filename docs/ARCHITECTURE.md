# Canid Evolutionary Genomics Pipeline — Architecture

> Working repo name: `canid-genomics` · importable package: `canidae` · codename: **CANIS**
> (Canid Ancestry, Networks, Introgression & Selection)

A standalone, publication-quality computational genomics platform for comparative
evolutionary analysis of *Canis* and related genera. Independent codebase; may reuse
*ideas* and utility patterns from the human archaeogenomics work but shares no code by
default. Designed to scale from a handful to **thousands** of whole genomes across
**dozens** of published datasets.

Scope: comparative evolutionary genomics only. **Not** for veterinary diagnostics or
conservation-management decisions.

> **Current operating model (v0.2):** CANIS is laptop-first. It uses a local native
> executor, bounded workers/memory/disk, sequential chromosome analysis, metadata-only safe
> resume, and reduced-panel/same-build VCF entry points. Cluster, automatic liftover, and
> large raw-read workflows are not default execution paths; see [Laptop operations](LAPTOP_OPERATIONS.md).

---

## 1. Design principles

1. **Loose coupling via typed artifacts.** Modules never import or call each other.
   Each stage consumes and produces *artifacts* (typed file handles + metadata)
   registered in a `DataStore`. The only shared contracts are the domain model and the
   artifact schemas. This is what makes the project expandable without restructuring.
2. **Configuration over code.** Every run is fully described by layered YAML configs
   validated by Pydantic. No analysis parameter is hard-coded in a module.
3. **Reproducibility is a first-class output, not an afterthought.** Every stage emits a
   provenance record: tool versions, exact command lines, input content hashes, config
   snapshot, container digest, git commit, RNG seeds, wall-clock, and host. A run is
   re-runnable from its manifest alone.
4. **Wrap, don't reimplement.** The genomics community has battle-tested C/C++ tools
   (bwa-mem2, GATK, bcftools, PLINK2, ADMIXTURE, Dsuite, IQ-TREE). We orchestrate them
   through a uniform runner abstraction and add value in the *glue*, the *domain model*,
   the *harmonization*, and the *analysis/reporting* — not by rewriting aligners.
5. **Scale by design.** Genotype matrices for thousands of samples never fit in memory
   naively. Use chunked, compressed on-disk formats (CRAM, BCF, Zarr) and out-of-core /
   distributed compute (Dask, cluster schedulers) from the start.
6. **Fail loudly and early.** Validate inputs against schemas before launching
   multi-hour jobs. Every external command is checked for exit status, expected outputs,
   and (where cheap) output sanity.

---

## 2. Layered architecture

```
Interfaces        CLI (Typer) · Python API · Notebooks · HTML/PDF reports
Orchestration     Workflow DAG · Snakemake/Nextflow backend · executors (local/SLURM/cloud) · GPU
Pipeline modules  (1) acquisition (2) qc (3) processing (4) popgen
                  (5) phylogenetics (6) introgression (7) comparative (8) geographic (9) reporting
Core foundation   config · provenance/logging · domain model · tool runners · datastore · registry
External          SRA/ENA/Dog10K/Dryad · reference genomes · external bioinformatics binaries
```

### 2.1 Core foundation (`canidae.core`)

The layer every module depends on and that depends on no module.

| Sub-package | Responsibility |
|---|---|
| `core.config` | Pydantic models for global + per-stage config; layered YAML loading (defaults → project → dataset → CLI override); schema validation; config hashing. |
| `core.logging` | Structured logging (JSON + rich console); per-run log file under the run directory; log level from config; correlation IDs per stage. |
| `core.provenance` | `ProvenanceRecord` dataclass; captures tool versions, argv, input/output hashes, config snapshot, container digest, git SHA, seeds, timing, host. Writes a `manifest.json` (and human-readable `manifest.md`) per run. |
| `core.model` | The domain model (§4): `Sample`, `Individual`, `Population`, `Taxon`, `Dataset`, `Cohort`, `ReferenceGenome`, `GenomicInterval`, `Artifact`, `Callset`, `AnalysisResult`. Dataclasses + Pydantic validation. |
| `core.datastore` | Content-addressable workspace. Typed `Artifact` handles for FASTQ/CRAM/BCF/PLINK/Zarr; path layout; caching + skip-if-exists; integrity (md5/sha256, `.done` sentinels). |
| `core.runtime` | `ToolRunner` abstraction: run an external command locally, in a container (Apptainer/Docker), or as a scheduler job (SLURM). Declares CPU/mem/GPU/time; captures stdout/stderr; retries; dry-run. `ToolSpec` records the required binary + version constraint. |
| `core.registry` | Plugin registry (`register_stage`, `register_dataset`, `register_reference`, `register_analysis`) using Python entry points so third parties / future modules self-register. |
| `core.errors` | Typed exception hierarchy (`ConfigError`, `MissingToolError`, `IntegrityError`, `StageInputError`, `ExternalToolError`). |
| `core.parallel` | Thin helpers over `concurrent.futures` and Dask; scatter/gather by `GenomicInterval`; bounded worker pools; progress. |

### 2.2 Stage contract (`canidae.core.stage`)

Every pipeline module implements one or more `Stage`s against a single interface:

```python
class Stage(ABC):
    name: ClassVar[str]
    config_model: ClassVar[type[StageConfig]]

    @abstractmethod
    def required_inputs(self) -> list[ArtifactSpec]: ...
    @abstractmethod
    def produced_outputs(self) -> list[ArtifactSpec]: ...
    def validate_inputs(self, ctx: RunContext) -> None: ...   # schema + existence checks
    @abstractmethod
    def run(self, ctx: RunContext) -> StageResult: ...        # does the work
```

* `ArtifactSpec` = an artifact *kind* (e.g. `CRAM`, `JOINT_VCF`, `PLINK_BED`) plus role
  ("input alignments"). It declares what a stage needs/emits **without** naming files.
* `RunContext` bundles the resolved config, the `DataStore`, a logger, the provenance
  writer, and the resource/executor handle. It is the *only* thing a stage's `run`
  receives — no globals, no cross-module imports.
* `StageResult` returns produced `Artifact`s + a `ProvenanceRecord` + metrics. The
  orchestrator registers the outputs in the `DataStore` so downstream stages can find
  them by kind/role. **This is the interface between modules.**

The upshot: adding a module = writing a `Stage` subclass and a config model, then
registering it. Nothing else in the tree changes.

---

## 3. Pipeline modules

Each is a subpackage under `canidae.stages.*`, exposing one or more `Stage`s.

### (1) `stages.acquisition` — dataset acquisition
* Download public data from **SRA/ENA** (`prefetch`/`fasterq-dump`, or ENA/EBI HTTP/
  Aspera), plus dataset-specific fetchers (Dog10K portal, Dryad, figshare) behind a
  common `DatasetFetcher` interface registered per dataset.
* **Metadata harvest** from ENA/BioSample/BioProject + curated per-study overrides
  (species, subspecies, population label, lat/long, source study DOI, platform,
  reported coverage, sex). Normalized against a controlled vocabulary.
* **Versioning**: each dataset has an immutable version tag + checksum manifest.
* **Integrity**: verify md5/ENA-reported hashes on download; quarantine mismatches.
* **Organization**: deterministic on-disk layout `datasets/<dataset>/<version>/…`.
* Outputs: `RawReads` (FASTQ/CRAM) artifacts + a validated `SampleSheet`.

### (2) `stages.qc` — quality control
* Read QC (fastp/FastQC) and aggregate (MultiQC).
* Coverage & mapping stats (samtools stats/mosdepth) once aligned.
* Missingness per sample/site; heterozygosity outliers; sex check.
* Relatedness / duplicate detection (KING or PLINK `--genome`) to flag replicate
  individuals and cryptic kinship across studies.
* Contamination / species check (e.g. mash/kmer screen; verifyBamID for humans has a
  canid analog via het-based estimators).
* Sample filtering with configurable thresholds → a **QC pass/fail table** and a
  filtered `Cohort`. Outputs QC metrics + a QC report fragment.

### (3) `stages.processing` — genome processing
* **Alignment** to a configurable `ReferenceGenome` (bwa-mem2 for short reads;
  minimap2 for long reads) → sorted, dup-marked **CRAM**.
* Post-processing: duplicate marking (samblaster/Picard), optional BQSR.
* **Variant calling**: pluggable backends — GATK HaplotypeCaller→GenomicsDB→
  GenotypeGVCFs, bcftools, or **DeepVariant** (GPU). Joint genotyping via GLnexus or
  GenomicsDB, scattered by `GenomicInterval`.
* **Variant filtering**: hard filters / VQSR; site- and genotype-level; biallelic SNP
  extraction; LD pruning outputs for downstream popgen.
* **Cross-study harmonization** (critical): reconcile datasets called against different
  references or site sets. Strategy = *lift/normalize to one canonical reference*
  (liftOver/CrossMap where needed) then restrict to a shared, well-covered, filtered
  site panel; track and correct for reference bias and batch effects; document every
  merge decision in provenance. Outputs a harmonized, analysis-ready `Callset`
  (BCF + PLINK2 pgen + optional Zarr).

### (4) `stages.popgen` — population genomics
* PCA (PLINK2 / scikit-allel randomized SVD / PCAngsd for low-coverage GLs).
* ADMIXTURE / fastSTRUCTURE / sNMF; cross-validation for K selection.
* Clustering (K-means / hierarchical / UMAP over PCs).
* Genetic distance matrices; pairwise F_ST (Weir–Cockerham).
* Diversity: π, θ, Tajima's D, individual heterozygosity, ROH — computed with **pixy**
  (correct handling of missing/invariant sites) or scikit-allel.
* Outputs `AnalysisResult`s (tables + figure-ready arrays).

### (5) `stages.phylogenetics` — evolutionary relationships
* Distance trees: Neighbor-Joining from genetic-distance matrices.
* Maximum-likelihood trees (IQ-TREE 2 / RAxML-NG) from SNP alignments; ascertainment-
  bias correction for SNP-only data.
* Bootstrap / ultrafast bootstrap support; consensus trees.
* Allele-frequency trees with migration edges via **TreeMix** (bridges into module 6).
* Outputs Newick trees + support values + rendered figures.

### (6) `stages.introgression` — hybridization & gene flow
* D-statistics (ABBA-BABA) and f3/f4 via **Dsuite** (VCF-native, fast) and/or
  **ADMIXTOOLS 2**; block-jackknife SE.
* f-branch statistics; admixture-graph fitting (qpGraph/`admixtools2`).
* Genome-wide introgression scans (f_d, f_dM, d_f in windows) for localized signals.
* Local-ancestry inference (RFMix / Loter / AncestryHMM) per admixed individual.
* Outputs statistic tables, window scans (BED-like), and local-ancestry karyograms.

### (7) `stages.comparative` — cross-population synthesis
* Assembles population-level comparisons: diversity summaries, shared-drift/outgroup-f3
  matrices, pairwise F_ST heatmaps, cross-population contrasts.
* Pure *aggregation/derivation* layer over modules 4–6; produces the tables that
  reporting turns into figures. No new external tools.

### (8) `stages.geographic` — spatial genetic structure
* Map sample localities (geopandas + cartopy/plotly).
* Regional ancestry summaries (admixture proportions by region/biome).
* Isolation-by-distance (Mantel), spatial PCA, and effective-migration surfaces
  (**FEEMS**). Outputs maps + spatial statistics.

### (9) `stages.reporting` — automated reporting
* Publication-quality static figures (matplotlib, consistent theme, colorblind-safe,
  vector SVG/PDF).
* Interactive visualizations (plotly/bokeh); optional lightweight dashboard.
* Composed **HTML** report (Jinja2), **PDF** (WeasyPrint/Quarto), **Markdown** summary,
  and exportable tables (CSV/Parquet/XLSX).
* Consumes `AnalysisResult`s only — reporting reads the DataStore, never re-computes.

### (10) Future modules (design accommodates without restructuring)
Selection scans (`stages.selection`: iHS/XP-EHH/PBS via selscan, SweepFinder2),
demographic inference (`stages.demography`: PSMC/MSMC2, SMC++, momi2, tsinfer/tsdate),
ancient canid genomes (`stages.ancient`: damage handling, pseudohaploid calls,
genotype-likelihood workflows), structural variation (`stages.sv`: Manta/Delly/
paragraph). Each is just another registered `Stage` set.

---

## 4. Domain model (`canidae.core.model`)

Dataclasses (frozen where immutable) with Pydantic validation:

* `Taxon` — species/subspecies with a controlled vocabulary (gray wolf & subspecies,
  coyote, red wolf, eastern wolf, domestic dog, village dog, dingo, golden jackal, …).
* `Individual` — a biological animal: id, taxon, population label, sex, geo-locality
  (lat/long, region), source study DOI, metadata bag.
* `Sample` — a sequencing library/run for an `Individual`: accession, platform, layout,
  reported/observed coverage, links to `RawReads`.
* `Population` — a named grouping used by analyses (many-to-many with individuals).
* `ReferenceGenome` — id, assembly (CanFam3.1 / UU_Cfam_GSD_1.0 / Dog10K), fasta+index
  artifacts, chrom naming scheme, liftover chains to other assemblies.
* `Dataset` — a published source: id, version, DOI, fetcher, sample list, checksums.
* `Cohort` — the working set for a run: selected individuals/populations + reference +
  filters. The primary input object to most stages.
* `GenomicInterval` — chrom/start/end for scatter/gather.
* `Artifact` — typed handle: kind (enum), path, format, checksum, producing stage,
  provenance id, schema version.
* `Callset` — a harmonized variant set: reference, samples, site count, format handles
  (BCF/pgen/Zarr), filter lineage.
* `AnalysisResult` — a typed result envelope: analysis name, parameters, tabular payload
  (path/DataFrame), figure specs, provenance id.

`SampleSheet` (a validated table, e.g. pandera schema over pandas) is the human-editable
entry point that materializes into `Individual`/`Sample`/`Cohort` objects.

---

## 5. Data flow

```
public sources ──▶ (1) RawReads + SampleSheet
                      │
                      ▼
                  (2) QC metrics ──▶ filtered Cohort
                      │
                      ▼
   reference ──▶ (3) CRAM ──▶ per-sample gVCF ──▶ joint Callset ──▶ harmonized Callset
                                                                        │
                        ┌───────────────────────────────┬──────────────┤
                        ▼                                ▼              ▼
                  (4) PopGen results            (5) Phylo trees   (6) Introgression stats
                        │                                │              │
                        └──────────────┬─────────────────┴──────────────┘
                                       ▼
                        (7) Comparative  +  (8) Geographic  ──▶  AnalysisResults
                                       │
                                       ▼
                             (9) Reporting  ──▶  HTML / PDF / MD / figures / tables
```

Everything crossing an arrow is a registered `Artifact` in the `DataStore`, described in
the run `manifest.json`. Re-running a stage with unchanged inputs+config is a no-op
(cache hit). Many users will **enter at the harmonized `Callset`** by ingesting already-
published VCFs — modules 4–9 run without ever touching raw reads.

---

## 6. Recommended libraries & tools

**Python core:** pydantic, typer (CLI), rich, PyYAML/tomllib, structlog, joblib,
pandas + pandera, numpy/scipy, attrs/dataclasses.
**Genomics in Python:** cyvcf2 & pysam (htslib bindings), scikit-allel, **sgkit**
(xarray/Dask/Zarr — the scalable backbone), msprime/tskit/demes (simulation & validation),
DendroPy / ete3 / Bio.Phylo (trees), pyd4/mosdepth wrappers.
**Geospatial/plots:** geopandas, cartopy, shapely, matplotlib, seaborn, plotly, bokeh,
Jinja2, WeasyPrint / Quarto.
**Scale/parallel:** Dask (+ dask-jobqueue), Zarr, numba, concurrent.futures.
**External binaries (wrapped):** sra-tools, fastp, FastQC, MultiQC, bwa-mem2, minimap2,
samtools/bcftools/htslib, samblaster/Picard, GATK4, GLnexus, DeepVariant, PLINK2,
vcftools, pixy, ADMIXTURE, PCAngsd, ANGSD, Dsuite, ADMIXTOOLS 2 (R), TreeMix, IQ-TREE 2,
RAxML-NG, RFMix/Loter, FEEMS, selscan (future), SMC++/MSMC2 (future).
**Orchestration:** Snakemake (Python-native, recommended default) or Nextflow/nf-core
(sarek for calling) as an execution backend behind `core.runtime`.
**Reproducibility:** conda/mamba + lockfiles, Apptainer/Docker images, per-tool version
pins, optional DVC for large-artifact tracking.
**Dev/test/CI:** pytest + hypothesis, coverage, ruff + black + mypy, pre-commit,
GitHub Actions, mkdocs-material (docs), Zenodo (release DOIs).

---

## 7. Computational bottlenecks & acceleration

| Stage | Bottleneck | Nature | Parallelization | GPU |
|---|---|---|---|---|
| Alignment | bwa-mem2 mapping | CPU, per-sample | Embarrassingly parallel across samples; multithread within | NVIDIA Parabricks `fq2bam` |
| Variant calling | HaplotypeCaller / DeepVariant | CPU (or GPU) | Scatter by `GenomicInterval`, gather | **DeepVariant / Parabricks on GPU — big win** |
| Joint genotyping | GenomicsDB / GLnexus merge | I/O + memory | Scatter by region | — |
| PCA | SVD of N×M genotype matrix | memory/CPU | Randomized SVD; PLINK2 multithread; Dask/sgkit out-of-core | cuML/RAPIDS randomized PCA |
| ADMIXTURE | EM over K | CPU | Multithreaded; parallel across K & seeds | GPU admixture variants (optional) |
| Distance matrix | pairwise O(N²) | CPU/BLAS | Blocked BLAS; Dask | GPU BLAS for large N |
| D/f-statistics | combinatorial over pops×blocks | CPU | Parallel block-jackknife; Dsuite optimized | — |
| Local ancestry | HMM per indiv×chrom | CPU | Parallel across indiv/chrom | — |
| ML phylogeny | tree search + bootstrap | CPU | Bootstrap replicates parallel; IQ-TREE threads | — |
| UMAP/t-SNE/clustering | neighbor graphs | CPU | — | cuML UMAP/KMeans/t-SNE |

**Python parallelism strategy:** the heavy numeric work lives in external C/C++ binaries
run as subprocesses, so the GIL is largely irrelevant — use **process-level** parallelism
(job scheduler / `ProcessPoolExecutor`) for scatter/gather, **Dask** for out-of-core
in-Python arrays, thread-level only for BLAS/htslib-internal threads. Reserve GPU for
DeepVariant (calling) and cuML (PCA/clustering/UMAP) where it pays off.

---

## 8. Repository structure

```
canid-genomics/
├── pyproject.toml                # PEP 621, src-layout, ruff/black/mypy/pytest config
├── README.md  LICENSE  CITATION.cff  CHANGELOG.md
├── environment.yml  conda-lock.yml
├── docs/                         # mkdocs-material: architecture, module guides, tutorials
│   ├── ARCHITECTURE.md  ROADMAP.md  index.md
├── configs/
│   ├── defaults.yaml             # base config
│   ├── datasets/*.yaml           # per-published-dataset descriptors
│   ├── references/*.yaml         # reference-genome descriptors
│   └── profiles/{local,slurm,cloud}.yaml
├── containers/                   # Apptainer/Docker defs per tool group
├── workflows/                    # Snakemake/Nextflow entrypoints (execution backend)
├── src/canidae/
│   ├── __init__.py  version.py  cli.py            # Typer CLI
│   ├── core/  {config, logging, provenance, model, datastore, runtime, registry,
│   │           parallel, errors, stage}.py
│   ├── stages/
│   │   ├── acquisition/  qc/  processing/  popgen/
│   │   ├── phylogenetics/  introgression/  comparative/  geographic/  reporting/
│   │   └── (future) selection/ demography/ ancient/ sv/
│   ├── datasets/          # registered DatasetFetcher plugins
│   ├── viz/               # figure themes + reusable plot functions
│   └── report/            # Jinja2 templates + report assembler
├── tests/
│   ├── unit/  integration/  fixtures/  simulated/   # msprime toy cohorts
│   └── conftest.py
├── data/                  # gitignored: raw/, interim/, processed/, results/, runs/
├── notebooks/             # exploratory, not part of the package
└── .github/workflows/ci.yml
```

---

## 9. Phased roadmap (MVP → publication platform)

**Phase 0 — Foundation & scaffolding.** src-layout package, `core.config/logging/
provenance/model/datastore/runtime/registry`, `Stage` contract, Typer CLI skeleton, conda
env + containers, CI, pytest harness, and **msprime-simulated toy cohorts** with a known
admixture graph as ground-truth fixtures. *Exit:* `canidae run --dry-run` builds and prints
a stage DAG; provenance manifest written.

**Phase 1 — MVP science from published VCFs.** `acquisition` (ENA/Dryad fetchers +
metadata) + `qc` + ingest of an existing harmonized `Callset` (e.g. Dog10K-derived) +
minimal `popgen` (PCA, F_ST, diversity) + a basic HTML report. Deliberately skips building
your own calling first, to reach a real, checkable result fast. *Exit:* one-command run
reproduces a PCA + F_ST + π report on a real canid dataset.

**Phase 2 — Full population genomics + reporting + geographic.** ADMIXTURE + K-selection,
clustering, distance matrices, ROH; polished publication figure theme; `geographic` maps &
IBD. *Exit:* publication-grade popgen report with maps.

**Phase 3 — Evolutionary + introgression.** `phylogenetics` (NJ + IQ-TREE + bootstrap +
TreeMix) and `introgression` (Dsuite/ADMIXTOOLS2 D & f-stats, f-branch, window scans).
Validate against the Phase-0 simulated truth. *Exit:* recovers known simulated gene flow;
runs D-stats on real wolf/coyote/dog data.

**Phase 4 — Raw-read processing & scale.** `processing` (bwa-mem2 → CRAM → DeepVariant/GATK
→ GLnexus) with Snakemake/Nextflow backend; **cross-study harmonization**; scale-out to
thousands via sgkit/Zarr + Dask + SLURM/cloud. *Exit:* end-to-end FASTQ→report on a cluster
for a multi-hundred-genome cohort.

**Phase 5 — Advanced introgression & GPU.** Local-ancestry inference; genome-wide
introgression maps; GPU acceleration (DeepVariant, cuML PCA/UMAP). *Exit:* per-individual
ancestry karyograms; GPU path benchmarked.

**Phase 6 — Future modules & publication packaging.** `selection`, `demography`, `ancient`,
`sv`; full docs site, tutorials, Zenodo DOI, CITATION.cff, and a reproducible
"paper-in-a-repo" demonstrating a real comparative-genomics result. *Exit:* external user
reproduces a figure from a single config + `canidae run`.
