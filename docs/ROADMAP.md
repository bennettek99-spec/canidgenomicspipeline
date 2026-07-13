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

## Near-term quality work

1. Add more tested public preset metadata sheets and source-specific citation bundles.
2. Add a chunk-iterating analysis API on top of the delivered memory-mapped genotype backend for algorithms that do not require whole-matrix access.
3. Expand report interactivity with optional self-contained SVG/Canvas views without making Plotly or a web service mandatory.
4. Add input-specific QC modules for coverage/read quality when the user intentionally chooses a raw-read workflow.
5. Add optional executable parity tests against PLINK2, bcftools, and Dsuite when those binaries are available in CI.

## Deliberately external preparation steps

- Cross-assembly liftover and reference harmonization should be run, inspected, and recorded before CANIS receives a VCF.
- Full alignment/calling of large raw-read cohorts needs appropriate storage and compute planning; CANIS preserves wrappers but does not launch it as a default laptop action.
- Publication-ready claims need study-specific callable-site masks, sensitivity analyses, independent replication, and reviewed methods beyond a browser report.
