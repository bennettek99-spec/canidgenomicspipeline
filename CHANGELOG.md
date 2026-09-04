# Changelog

## Unreleased — hybrid-canid stages

### Improved

- **LD pruning is now laptop-fast.** `analysis_readiness` used a per-pair
  `np.corrcoef` Python loop (~220 s on a 9.6k-variant simulated cohort); the same
  greedy bp-window prune now standardizes the matrix once and scores each variant
  against its in-window retained set with one vector dot product plus a per-contig
  sliding window (~0.3 s, identical results on complete data, mean-imputed missing
  calls). Full `ingest → report` laptop pipelines drop from minutes to seconds.
- **The report now states its findings.** The executive summary names the retained
  analysis panel, LD-pruning yield, PC1 separation and variance explained, strongest
  and weakest F_ST pairs, selected admixture K with its CV margin, diversity ranking,
  top F_ROH sample, and Mantel isolation-by-distance. PCA/FST sections carry the same
  one-line captions, admixture shows CV errors per K, and ROH renders its per-sample
  table alongside the barplot.

### Fixed

- **Corrected a genotype-scale bug in the hybrid likelihood layer.** `analysis.genotypes.dosage`
  returned an ALT dosage in `[0, 1]` (0 / 0.5 / 1), but every consumer — the binomial
  log-pmf, the Laplace-smoothed allele frequencies, the shrunk panel frequencies, and the
  three-way ML grid — is written for ALT allele *counts* (0 / 1 / 2), the convention used
  everywhere else in CANIS. Two consequences: homozygous-ALT genotypes were scored on the
  heterozygote branch (`dosage > 1.0` was never true), and group allele frequencies were
  compressed into `[0, 0.5]` (a fixed-ALT panel returned 0.5 instead of ~0.95).
  `allele_count()` is now the primary helper and the likelihood layer uses it.
  - `reference_mixture` / NYC dog fractions are **unchanged** — `infer_dog_fraction`
    compares a genotype against group mean dosages on one consistent scale.
  - `breed_assign` leave-one-out top-1 accuracy on the real 722-genome panel improves
    from 0.409 to 0.520. Best-breed calls and every "no single breed supported" verdict
    are unchanged.
  - `multiway_admixture` wolf/dog fractions were inflated (roughly 2x on data with a real
    wolf component). On the eastern-coyote panel the qualitative result is unchanged —
    western controls at zero, wolf/dog still unresolvable at this marker count — but the
    pooled eastern dog fraction moves from 0.12 to 0.03.
  - Renamed for clarity at the call site: `wgs_dosage_matrix` → `wgs_allele_count_matrix`,
    `bridge_panel.dosage_vector` → `allele_count_vector`.

### Added

- Extracted hybrid-canid libraries under `canidae.analysis` and `canidae.io.indexed_vcf`; scripts are thin wrappers.
- Added stages `reference_mixture`, `breed_assign`, and `multiway_admixture` with NYC coydog and eastern coyote YAML presets.
- Report stage accepts hybrid-only recipes (optional PCA/FST/diversity) and renders mixture/breed/multiway sections with bridge-locus limitations.
- Added a deterministic synthetic bridge panel (`tests/simulated/make_hybrid_panel.py`) with known
  mixture fractions, plus unit, stage-DAG, and golden-snapshot coverage for all three hybrid stages.
- Added golden snapshots under `tests/golden/`: a CI-runnable synthetic set, and a
  study-data set that gates promotion locally and skips where `data/` is absent.
- Added offline coverage for the indexed-VCF transport (BGZF decoding, tabix index math,
  transfer budget) and an opt-in live-source test behind a new `network` marker, which is
  deselected by default.
- Added a chunk-iterating genotype API — `iter_genotype_chunks`, `chunked_allele_counts`,
  `chunked_alt_frequency` — that streams the memory-mapped backend a variant block at a time.
- Added optional PLINK2 and Dsuite parity tests that run when those binaries are on PATH
  and skip otherwise.
- Added source citation bundles under `configs/citations/`, wired into every public preset
  and rendered in the report's Sources section. Unrecorded DOIs are reported as
  "not recorded" rather than guessed.
- Added coverage over every shipped example preset and metadata sheet: each recipe
  instantiates its stages, configures no stage it does not run, and each sheet satisfies
  the sample-sheet schema. This corrected `configs/examples/nyc_coydog_samples.csv`, whose
  `hybrid` taxon was not a valid `CanidTaxon` and whose accession columns were named off
  the schema, so their values were being dropped on read.

## 0.3.0 - reliable, analysis-ready laptop pipeline


- Added exact indexed-panel transfer preflight, persistent checksummed HTTP range caching, bounded two-worker fetching, live ETA/progress state, and optional GQ/DP hard-call filters.
- Added an analysis-readiness gate with autosome filtering, optional LD pruning, population-size and near-duplicate audits, and explicit panel/quality limitations.
- Added memory-mapped genotype artifacts for larger callsets while preserving the existing stage contracts.
- Required deliberate outgroups for D/f3 analyses, moved uncertainty defaults to chromosome/physical blocks, added D-statistic FDR values, and added replicated ADMIXTURE fitting diagnostics.
- Made safe-resume fingerprints depend on stage-relevant configuration and transitive code dependencies, while persisted artifact integrity now always uses full-content checksums.
- Added resumable provenance with resolved configuration snapshots, interrupted-stage recovery records, and post-promotion paths that never leak staging locations.
- Replaced optimistic UI transfer guesses with exact preflight, validated optional-analysis inputs, and served outputs through allowlisted local HTTP routes.
- Made type checking a blocking CI gate under an incremental legacy baseline; expanded reliability, readiness, cache, quality-filter, mmap, and integration coverage.

## 0.2.0 — laptop-safe pipeline

- Added metadata-only safe-resume cache fingerprints, transactional stage promotion, atomic artifact manifests, and failure/skip/recovery provenance.
- Added enforced sample/site QC filtering with compact filtered VCFs and exact exclusion reports.
- Added chromosome/fixed-Mb D/f-statistic blocks and chromosome-reset local ancestry HMMs.
- Added integrated indexed-VCF reduced-panel acquisition with 2K/10K/25K presets, transfer safeguards, checksums, and cleanup.
- Added streamed paired, interleaved, and single-end gzip FASTQ acquisition with pair-preserving downsampling, archive-size gates, checksum verification, and automatic alignment-manifest handoff.
- Added panel-relative diversity labels and callable-site gating for whole-genome demographic estimators.
- Added lightweight same-build harmonization diagnostics and allele-orientation reconciliation.
- Added laptop resource controls, disk preflight, cleanup, pause/resume/cancel markers, and `canidae ui`.
- Expanded the self-contained report into a navigable analysis package with source/provenance and limitations.

## 0.1.0

- Initial CANIS population-genomics, introgression, local-ancestry, processing, and reporting foundation.
