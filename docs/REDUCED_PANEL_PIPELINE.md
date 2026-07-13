# Integrated reduced-panel acquisition

`reduced_panel` is a normal CANIS acquisition stage for laptop-sized analysis of a public,
bgzip-compressed VCF with a tabix index. It fetches only the compressed byte ranges that
overlap a fixed SNP panel; it does not download the complete source VCF or raw reads.

The built-in presets are:

| Preset | Target SNPs | Typical use |
| --- | ---: | --- |
| `2k` | 2,000 | Fast connectivity and workflow smoke test |
| `10k` | 10,000 | Default exploratory laptop analysis |
| `25k` | 25,000 | Higher-resolution exploratory analysis |

## Run the integrated Red Wolf example

```powershell
canidae run -c configs/examples/redwolf_jackal_reduced_panel.yaml `
            -c configs/profiles/laptop.yaml
```

The first attempt downloads only the marker panel and small tabix index, then calculates an
exact merged byte-range plan without fetching source VCF ranges. If the estimate exceeds the configured confirmation
threshold, CANIS stops before downloading source VCF ranges and prints the exact estimate.
Review it, then set this one configuration field to proceed:

```yaml
stages:
  reduced_panel:
    confirm_large_transfer: true
```

The default hard ceiling is 9 GB and cannot be configured to 10 GB or more. Every byte-range
response must be HTTP 206; a server that ignores `Range` is rejected so it cannot trigger an
accidental full-source download.

## What the stage does

```text
remote indexed VCF
  -> download marker panel + tabix index
  -> estimate byte ranges and apply confirmation/ceiling policy
  -> validate selected sample names against VCF header
  -> fetch only target ranges through a persistent checksummed cache
  -> retain complete diploid biallelic SNP GT calls, optionally filtered by GQ/DP
  -> autosome/readiness filtering and population/duplicate audits
  -> compact .vcf.gz + normalized sample sheet + manifest
  -> QC and analysis stages
```

The sample sheet is the sample-selection input. Its `sample_id` values must exactly match
remote VCF sample names. Set `selected_samples` to select a deliberate subset from a larger
sheet; otherwise every row is included.

```yaml
stages:
  reduced_panel:
    sample_sheet: configs/my_candidates.csv
    selected_samples: [Wolf25, Wolf26, GoldenJackal01]
    preset: 10k
    max_download_bytes: 9000000000
    confirmation_threshold_bytes: 1000000000
    confirm_large_transfer: false
    range_workers: 2
    min_genotype_quality: 20
    min_genotype_depth: 5
```

## Safety and reproducibility

- The marker-panel publisher checksum is verified before use. An optional tabix-index
  checksum can also be set with `index_checksum` (`md5:<digest>` or `sha256:<digest>`).
- The compact VCF is written through a temporary neighboring file and promoted only once
  complete. Its SHA-256, source URLs, transfer estimate, actual bytes, validation exclusions,
  preset, and selected samples are written to `reduced_panel_<preset>.manifest.json`.
- Temporary marker-panel and index files are removed after every attempt. The retained compact
  VCF, manifest, and normalized sample sheet live under `data/.../store/reduced_panel/`.
- Verified source ranges persist under the configured cache root. An interrupted rerun reuses
  only files whose length and SHA-256 sidecar still match; corrupt or incomplete entries are
  discarded and fetched again. Fetching is parallel but buffered in small batches for laptop RAM.
- Live range counts, downloaded/reused bytes, retained SNPs, and ETA are written to
  `runs/<run-id>/progress.json` for the local UI.
- Extraction requires at least 1,000 retained sites by default (or the selected preset size if
  it is smaller). Increase `min_retained_sites` for a stricter requirement.

This path is deliberately exploratory: a fixed marker panel is suitable for PCA, distances,
and visual comparisons, but it is not a substitute for a callable-site whole-genome analysis.
