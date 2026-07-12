# Running the pipeline on real red wolf genomes

This turns the simulated red-wolf demo into a study on **real public genomes**. Because red
wolves are only deposited as **raw sequencing reads** (no public VCF), you must align + call
first — which needs Linux/bioconda tools (bwa-mem2, samtools, GATK). None of it uses a GPU;
it is CPU-only, but it is a *large* compute job.

## The cohort (all real accessions)

| Taxon | Runs | Study |
|---|---|---|
| red wolf (*Canis rufus*) | SRR7976417, SRR7976418, SRR7976419 | vonHoldt et al. 2016 |
| coyote | SRR12075083, SRR12075084, SRR12075106 | Rougeux et al. 2023 (PRJNA641325) |
| gray wolf (Alberta) | SRR12075090, SRR12075091, SRR12075093 | Rougeux et al. 2023 |
| eastern wolf | SRR12075101, SRR12075104, SRR12075105 | Rougeux et al. 2023 |
| golden jackal (outgroup) | SRR2149876 | Koepfli et al. 2015 |

The red wolves and references come from *different studies* — so this exercises the
cross-study design (`sample_sheet` + `harmonize`) on real data.

## Reality check (please read before you start)

- **Download:** ~40–55 Gbp/sample × 13 ≈ **150–250 GB** of gzipped FASTQ.
- **Working space:** CRAMs + gVCFs + bwa-mem2 index ≈ **1–2 TB** peak. Use a big scratch disk.
- **Compute (CPU only):** bwa-mem2 ≈ 2–6 CPU-h/genome; GATK HaplotypeCaller ≈ several CPU-h/
  genome. On a 12-core laptop expect **~1–3 days** wall-clock for the full cohort.
- **Laptop tip:** set `SUBSAMPLE=0.1` in the reads script for a ~5× proof-of-concept (hours,
  not days, and ~20 GB), or run the whole thing on a cloud VM. Downsampling weakens D-stat Z
  and local-ancestry resolution but still shows the structure.

## No Linux machine? You don't need one

Only the **reads → VCF** step needs Linux tools. The whole analysis half (PCA, ADMIXTURE,
D-stats, local ancestry, report) already runs on Windows, and the joint VCF is small
(~100s of MB) and portable. So you have three routes:

### Route A — WSL2 (free, on your Windows 11 laptop; best for a proof-of-concept)

WSL2 *is* Linux, running inside Windows. In **PowerShell (as admin)**:

```powershell
wsl --install -d Ubuntu     # then reboot if prompted; set a Linux username/password
```

Open the **Ubuntu** terminal and install Miniforge (conda/mamba) + run the bootstrap:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b && ~/miniforge3/bin/conda init bash && exec bash
cd "/mnt/c/Users/benne/OneDrive/Desktop/Claude Work/canid-genomics"   # your repo, seen from WSL
bash scripts/bootstrap.sh          # creates the env + installs the package
```

**Disk reality:** a full 13-genome run needs ~1–2 TB — more than an XPS 13 SSD. On the laptop,
keep it small: `SUBSAMPLE=0.1` **and** trim the cohort to ~2 samples per group (edit
`configs/examples/redwolf_samples.csv` and the accession list in
`scripts/fetch_redwolf_reads.sh`). That fits in ~50–100 GB and runs in a few hours. Put the
`data/` dir on an external SSD if space is tight (`export CANIDAE_DATA=/mnt/e/...` and set
`paths.data_root`).

### Route B — a cloud Linux VM (best for the full-coverage study)

Spin up a Linux VM with a big disk and run the exact same commands — no GPU needed:
- **AWS EC2** `c7i.4xlarge` (16 vCPU) + a 2 TB gp3 EBS volume, or **GCP** `c3-standard-16`.
- Cost is roughly **$5–20** for the full cohort (a day of compute + storage), less on spot/
  preemptible instances.
- `git clone` the repo, run `bash scripts/bootstrap.sh`, then the fetch + run commands below.
- When done, download just `data/store/report/report.html` (+ the result CSVs) and delete the VM.

### Route C — do only the reads→VCF in the cloud, analyze on Windows

Use a free web platform (e.g. **usegalaxy.org**) or a cloud VM to align + joint-call the 13
genomes, download the resulting **joint VCF**, then on your Windows laptop:

```powershell
canidae run -c configs\examples\redwolf_from_vcf.yaml --set stages.ingest.callset=C:\path\to\joint.vcf.gz
```

That config uses `ingest` (registers your VCF *and* the sample sheet), then runs the entire
analysis stack natively on Windows — exactly like the simulated demo did. The VCF's sample
names must match the `sample_id` column in `configs/examples/redwolf_samples.csv`.

## Setup (WSL2 or Linux)

```bash
# 1. tools (bioconda)
mamba create -n canis-tools -c bioconda -c conda-forge \
    sra-tools bwa-mem2 samtools bcftools gatk4 seqtk curl
mamba activate canis-tools

# 2. the package itself (in its venv/conda env)
pip install -e ".[dev,analysis]"

# 3. reference + reads
bash scripts/fetch_reference.sh
SUBSAMPLE=0.1 bash scripts/fetch_redwolf_reads.sh      # drop SUBSAMPLE for full coverage

# 4. run everything -> one HTML report
canidae run -c configs/examples/redwolf_real.yaml -c configs/profiles/laptop.yaml
# -> data/store/report/report.html
```

`executor.resume: true` uses input/config/code fingerprints and validated promoted outputs.
If a long run is interrupted, re-run the same command (or use the same `--run-id`) and only
unchanged completed stages are reused; inspect the manifest for exact recovery directions.

## What you'll get

The same 17-stage report as the simulated demo, but on real animals: PCA, ADMIXTURE, F_ST,
neighbor-joining tree, **D-statistics** (is red wolf admixed between coyote and gray wolf?),
**local-ancestry karyograms** (the chromosomal coyote/wolf mosaic), diversity, and geography.
This is, in miniature, the exact analysis behind the vonHoldt (2016) vs. Hohenlohe (2017) /
Waples (2018) debate — and the pipeline lets you see which signals actually support each side.

## Caveats specific to real data

- **Reference bias:** all samples are mapped to the *dog* assembly (CanFam3.1), which is
  closer to wolves than coyotes; this can bias coyote-side statistics. A wolf/coyote pangenome
  reference reduces it.
- **Batch effects:** cross-study samples differ in coverage, platform, and processing —
  `load_genotypes`' MAF/missingness filters and `harmonize` mitigate but don't erase this.
- **Small n:** 3 per group is enough to *demonstrate* the signals, not to publish; add more of
  the PRJNA641325 runs (up to 10 coyote / 10 gray wolf / 5 eastern wolf) for real power.
- **The hard question** — "hybrid origin vs. ancient distinct lineage with recent
  introgression" — needs demographic modeling of *divergence times* (SMC++/momi2), not just
  these summary statistics.
