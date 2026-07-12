# Changelog

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
