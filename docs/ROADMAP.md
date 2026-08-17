# CANIS operating roadmap

CANIS is intentionally laptop-first. The roadmap is not a promise to automatically scale raw-read workflows to a large cluster; it prioritizes correct, auditable reduced-panel and same-reference VCF analyses on ordinary local hardware.

## Delivered laptop foundation

- Metadata-only safe resume keyed to input/config/code fingerprints.
- Atomic outputs, checksum validation, complete manifests, and recovery directions.
- Enforced sample/site filtering with retained/excluded datasets made explicit.
- One-chromosome-at-a-time local ancestry and physical-block introgression uncertainty.
- Indexed remote VCF reduced-panel acquisition with a hard sub-10-GB ceiling and human confirmation for larger estimates.
- Explicit panel-relative diversity and callable-site gating of whole-genome demographic metrics.
- Same-build harmonization checks, orientation diagnostics, and no automatic liftover.
- Resource limits, disk preflight, conservative threads, intermediate cleanup, pause/resume/cancel markers, and browser UI controls.
- A self-contained, navigable HTML report with tables, figures, limitations, provenance, and source accessions.
- Exact indexed-range preflight, persistent verified range reuse, live transfer ETA, and bounded concurrent buffering.
- Analysis-readiness filtering/audits, explicit outgroups, block-aware uncertainty, and multiple-testing output.
- Memory-mapped genotype artifacts with automatic size-based selection.
- Resolved-config snapshots, interrupted-run recovery, and blocking incremental type checks in CI.

## Delivered hybrid-canid diagnostics (v0.4 path)

- Library extraction under `canidae.analysis` and `canidae.io.indexed_vcf` (scripts are thin wrappers).
- Registered stages: `reference_mixture`, `breed_assign`, `multiway_admixture`.
- Presets: `configs/examples/nyc_coydog_validation.yaml`, `configs/examples/eastern_coyote_ancestry.yaml`.
- Report sections for mixture / breed / multiway results with explicit bridge-locus limitations.
- Deterministic synthetic bridge panel with known mixture fractions, backing unit, stage-DAG, and golden-snapshot tests.
- Golden snapshots for both the synthetic fixture (runs in CI) and the real study panels (skipped where `data/` is absent).
- Source citation bundles under `configs/citations/`, referenced by every public preset and rendered in the report.
- Chunk-iterating genotype API over the memory-mapped backend.
- Optional PLINK2 / Dsuite parity tests, and an opt-in `network` marker for live-source transport tests.
- Remaining: optional `bridge_panel` acquisition stage to replace ad-hoc data prep paths.

## Near-term quality work

1. Expand report interactivity with optional self-contained SVG/Canvas views without making Plotly or a web service mandatory.
2. Add input-specific QC modules for coverage/read quality when the user intentionally chooses a raw-read workflow.
3. Extend executable parity coverage to bcftools, and run the existing PLINK2/Dsuite parity tests in a CI job that installs those binaries.
4. Move the remaining whole-matrix stages onto the chunk-iterating API where the algorithm allows it.
5. Fill the `doi:` fields in `configs/citations/` that this repository does not yet record, so every source resolves to a publication.


## Deliberately external preparation steps

- Cross-assembly liftover and reference harmonization should be run, inspected, and recorded before CANIS receives a VCF.
- Full alignment/calling of large raw-read cohorts needs appropriate storage and compute planning; CANIS preserves wrappers but does not launch it as a default laptop action.
- Publication-ready claims need study-specific callable-site masks, sensitivity analyses, independent replication, and reviewed methods beyond a browser report.
