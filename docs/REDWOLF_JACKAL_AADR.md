# Red Wolf versus golden jackal: laptop SNP panel

This workflow builds a compact, **AADR-style reduced genotype dataset** from the public
722-canid project and runs a deliberately small comparison on a laptop. It does not download
or store the complete 309 GB source VCF, and it does not download raw FASTQs.

> **Integrated path:** the same indexed extraction is now available as the normal
> `reduced_panel` pipeline entry point. For the 2K/10K/25K presets, transfer estimate and
> confirmation flow, use [`REDUCED_PANEL_PIPELINE.md`](REDUCED_PANEL_PIPELINE.md) and
> `configs/examples/redwolf_jackal_reduced_panel.yaml`. The scripts below remain useful for
> standalone/manual preparation.

“AADR-style” describes the data reduction: a fixed set of informative, genome-wide SNPs and
compact per-sample genotype calls. The output remains `VCF.gz`, rather than EIGENSTRAT,
because that is the native input of this repository.

## Source data

- Project: NCBI BioProject `PRJNA448733`, analysis `SRZ189891`.
- Source VCF: `722g.990.SNP.INDEL.chrAll.vcf.gz`, aligned to CanFam3.1.
- Target positions: an evenly spaced subset of the CanFam3.1 Illumina CanineHD panel.
- Samples:
  - `Wolf25`: Red Wolf, BioSample `SAMN02921317`, source coverage 28.5x.
  - `Wolf26`: Red Wolf, BioSample `SAMN02921318`, source coverage 7.1x.
  - `GoldenJackal01`: golden jackal, BioSample `SAMN03366713`, source coverage 5.5x.

Two Red Wolves are included because the existing PCA and neighbor-joining stages require at
least three individuals. This remains a small exploratory comparison, not a population-scale
study.

## Requirements

The preferred preparation script runs directly on Windows with Python and has no additional
binary dependencies. It reads the tabix index itself and enforces a hard 9 GB transfer ceiling.
The original `bcftools` implementation remains available for Linux/WSL2 users.

The remote VCF is bgzip-compressed, tabix-indexed, and served with HTTP byte ranges. The script
therefore requests only compressed blocks surrounding the selected markers. Network transfer
depends on block density and caching, but is expected to be vastly smaller than 309 GB.

## Prepare the dataset

From the repository root on Windows:

```powershell
.\.venv\Scripts\python.exe scripts\prepare_redwolf_jackal_aadr.py
```

Or under Linux/WSL2 with `bcftools`:

```bash
conda activate canis
bash scripts/prepare_redwolf_jackal_aadr.sh
```

Default output:

```text
data/redwolf_jackal_aadr/
├── CFA31_IlluminaHD.vcf.gz
├── caninehd_autosomal_positions.tsv
├── caninehd_even_25000.regions.tsv
├── redwolf_jackal_caninehd25000.vcf.gz
└── redwolf_jackal_caninehd25000.manifest.json
```

For a smaller network and runtime test, request approximately 10,000 sites:

```bash
TARGET_SITES=10000 bash scripts/prepare_redwolf_jackal_aadr.sh
```

On Windows, the equivalent is:

```powershell
.\.venv\Scripts\python.exe scripts\prepare_redwolf_jackal_aadr.py --target-sites 10000
```

The supplied analysis config points to the default 25,000-site filename. If `TARGET_SITES` is
changed, override the callset path when running CANIS.

## Run the prepared comparison

The preparation script does **not** launch CANIS. After it succeeds, run:

```bash
canidae run \
  -c configs/examples/redwolf_jackal_aadr.yaml \
  -c configs/profiles/laptop.yaml
```

For a 10,000-site preparation:

```bash
canidae run \
  -c configs/examples/redwolf_jackal_aadr.yaml \
  -c configs/profiles/laptop.yaml \
  --set stages.ingest.callset=data/redwolf_jackal_aadr/redwolf_jackal_caninehd10000.vcf.gz
```

The config runs only the stages that make sense for this tiny cohort: call-rate QC, genotype
loading, individual distance, two-component PCA, exploratory FST/diversity summaries, a
three-tip neighbor-joining tree, and the HTML report. It intentionally omits ADMIXTURE,
D-statistics, local ancestry, demography, selection scans, and geography.

Outputs are isolated under:

```text
data/redwolf_jackal_aadr_workspace/store/
runs/redwolf_jackal_aadr/
```

The report will be at:

```text
data/redwolf_jackal_aadr_workspace/store/report/report.html
```

## Interpretation limits

- This panel is designed for a fast visual/genetic-distance comparison.
- The jackal population has one individual and the Red Wolf population has two; FST and
  diversity values are exploratory and should not be treated as population estimates.
- SNP-array ascertainment favors variants selected for the canine IlluminaHD panel.
- The two Red Wolves are modern samples and do not represent historical Red Wolf diversity.
- The workflow cannot answer Red Wolf hybrid-origin or introgression questions; those require
  coyotes, gray wolves, an outgroup, and larger sample sizes.
