# Run the real red-wolf study on your laptop — exact steps

A 9-genome real cohort (2 red wolf + 2 coyote + 2 gray wolf + 2 eastern wolf + 1 golden
jackal outgroup) sized to run on a Dell XPS 13 through WSL2. Everything below is copy-paste.

**Before you start**
- WSL2 + Ubuntu installed (`wsl --install -d Ubuntu` in an admin PowerShell, then reboot).
- **Free disk:** ~80 GB for a first "does it work" run (`SUBSAMPLE=0.1`), ~150–200 GB for good
  results (`SUBSAMPLE=0.3`). Check `Cleanup` at the bottom if space is tight.
- **Time:** ~a few hours. The laptop config uses the **fast `bcftools_call`** path (one
  multi-sample variant call instead of per-sample GATK), so alignment is now the main cost
  (~10–30 min/genome downsampled). Every stage is checkpointed (`resume: true`), so you can
  close the laptop and re-run the same command to continue.

---

### Step 1 — open Ubuntu and install Miniforge (conda)

Open the **Ubuntu** app (Start menu). Then:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
exec bash
```

### Step 2 — go to the project (it lives on your Windows drive)

```bash
cd "/mnt/c/Users/benne/OneDrive/Desktop/Claude Work/canid-genomics"
```

### Step 3 — build the environment (bwa-mem2, GATK, samtools, + the package)

```bash
bash scripts/bootstrap.sh
conda activate canis
```

This takes ~10–20 min the first time. When it finishes, `canidae version` should print a
version number.

### Step 4 — download the dog reference genome and index it

```bash
bash scripts/fetch_reference.sh
```

~2.4 GB download + a ~15–30 min bwa-mem2 index (uses ~25 GB). One-time.

### Step 5 — download the 9 red-wolf-cohort genomes (downsampled)

```bash
SUBSAMPLE=0.1 bash scripts/fetch_redwolf_reads.sh configs/examples/redwolf_laptop_samples.csv
```

`SUBSAMPLE=0.1` keeps ~10% of reads (≈2–3× coverage) → ~15–20 GB, fast, rough results. For
cleaner statistics use `SUBSAMPLE=0.3` (≈8×, ~50 GB) if you have the disk.

### Step 6 — run the whole pipeline → one HTML report

```bash
canidae run -c configs/examples/redwolf_laptop.yaml -c configs/profiles/laptop.yaml
```

It runs 16 stages: align → **bcftools_call** → PCA/ADMIXTURE/F_ST/tree/D-stats/local
ancestry/… → report. If it stops (closed laptop, out of disk), just run the **same command
again** — it resumes from where it left off.

### Step 7 — open the report

The report is written under the project on your Windows drive. From Ubuntu:

```bash
explorer.exe "$(wslpath -w data/store/report/report.html)"
```

or just double-click it in Windows Explorer at:
`C:\Users\benne\OneDrive\Desktop\Claude Work\canid-genomics\data\store\report\report.html`

You'll get the same 14-section report as the simulated demo — PCA, ADMIXTURE, the
neighbor-joining tree, **D-statistics** (is red wolf admixed?), **local-ancestry karyograms**
(the coyote/gray-wolf chromosomal mosaic), diversity, and geography — but on real animals.

---

## Want better results later?

Re-run at higher coverage on the same machine (needs the disk), or move to a cloud VM
(`docs/REAL_DATA_REDWOLF.md`, Route B). To change coverage you must re-fetch + re-align:

```bash
rm -rf data/redwolf data/store            # clear reads + results (keeps the reference)
SUBSAMPLE=0.3 bash scripts/fetch_redwolf_reads.sh configs/examples/redwolf_laptop_samples.csv
canidae run -c configs/examples/redwolf_laptop.yaml -c configs/profiles/laptop.yaml
```

To add more samples for real statistical power, add rows to
`configs/examples/redwolf_laptop_samples.csv` (more coyote/gray-wolf runs from PRJNA641325 —
see `configs/examples/redwolf_samples.csv` for the full list) and re-run.

## Troubleshooting

- **"No space left on device":** you ran out of disk. Use `SUBSAMPLE=0.1`, run the cleanup
  step below between/after stages, or put the whole repo on an external SSD (clone/copy it to
  e.g. `/mnt/e/canid-genomics` and run from there — all paths are relative to the repo).
- **A download stalls:** re-run Step 5 — it skips files already fully downloaded.
- **Calling/alignment feels slow:** whole-genome scale is just heavy; check progress with
  `ls -la data/store/align/` and `ls -la data/store/bcftools_call/`. Fewer samples / lower
  `SUBSAMPLE` = faster.
- **A stage fails:** the run reports which one; fix its input and re-run the same command
  (completed stages are skipped).

## Cleanup (reclaim disk after the run)

```bash
rm -rf data/redwolf/reads          # the FASTQs (biggest; safe once alignment finished)
rm -rf data/store/align            # the CRAMs (safe once calling finished)
# keep data/store/report and the result CSVs — that's your study output
```
