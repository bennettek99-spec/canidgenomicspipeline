#!/usr/bin/env bash
# Download the dog reference genome (CanFam3.1 — the assembly the published red-wolf/canid
# studies used) and build the indices the align/call stages need. Requires: curl, bwa-mem2,
# samtools, gatk (all via bioconda; see docs/REAL_DATA_REDWOLF.md).
set -euo pipefail

REFDIR="${1:-data/references}"
mkdir -p "$REFDIR"
FA="$REFDIR/CanFam3.1.fa"
URL="https://ftp.ensembl.org/pub/release-104/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz"

if [ ! -s "$FA" ]; then
  echo ">> downloading CanFam3.1"
  curl -L --retry 3 -o "$FA.gz" "$URL"
  gunzip -f "$FA.gz"
fi

echo ">> samtools faidx"
samtools faidx "$FA"
echo ">> gatk sequence dictionary"
gatk CreateSequenceDictionary -R "$FA" -O "${FA%.fa}.dict" || true
echo ">> bwa-mem2 index (this takes a while + ~30 GB)"
bwa-mem2 index "$FA"
echo "Reference ready: $FA"
