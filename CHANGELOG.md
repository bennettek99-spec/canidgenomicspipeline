# Changelog

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
