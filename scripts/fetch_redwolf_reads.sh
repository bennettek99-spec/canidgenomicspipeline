#!/usr/bin/env bash
# Fetch paired-end FASTQs from ENA (direct HTTP — no sra-tools needed) for the accessions in
# a sample sheet's `sample_id` column, and write the fastq manifest the `align` stage consumes.
#
#   bash scripts/fetch_redwolf_reads.sh [SAMPLE_SHEET_CSV]
#   SUBSAMPLE=0.3 bash scripts/fetch_redwolf_reads.sh configs/examples/redwolf_laptop_samples.csv
#
# SUBSAMPLE keeps a fraction of reads (needs seqtk) so a laptop can cope; 0 = full coverage.
# WARNING: whole-genome data is large. Full coverage is ~40-55 Gbp/sample.
set -euo pipefail

SHEET="${1:-configs/examples/redwolf_samples.csv}"
OUT="${OUT:-data/redwolf/reads}"
MANIFEST="${MANIFEST:-data/redwolf/fastq_manifest.csv}"
SUBSAMPLE="${SUBSAMPLE:-0}"
mkdir -p "$OUT" "$(dirname "$MANIFEST")"

# accessions = first column (sample_id) of the sheet, minus the header
mapfile -t ACCESSIONS < <(tail -n +2 "$SHEET" | cut -d, -f1 | tr -d '\r')
echo ">> ${#ACCESSIONS[@]} accessions from $SHEET  (SUBSAMPLE=$SUBSAMPLE)"

echo "sample_id,fastq1,fastq2" > "$MANIFEST"
for acc in "${ACCESSIONS[@]}"; do
  echo ">> $acc"
  urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp&format=tsv" \
         | tail -n +2 | cut -f1)
  i=1
  for u in ${urls//;/ }; do
    dest="$OUT/${acc}_${i}.fastq.gz"
    [ -s "$dest" ] || curl -L --retry 3 -o "$dest" "https://${u}"
    i=$((i+1))
  done
  r1="$OUT/${acc}_1.fastq.gz"; r2="$OUT/${acc}_2.fastq.gz"
  if [ "$SUBSAMPLE" != "0" ]; then
    for f in "$r1" "$r2"; do
      seqtk sample -s100 "$f" "$SUBSAMPLE" | gzip > "${f%.fastq.gz}.sub.fastq.gz"
      mv "${f%.fastq.gz}.sub.fastq.gz" "$f"
    done
  fi
  echo "${acc},${r1},${r2}" >> "$MANIFEST"
done
echo ">> reads in $OUT ; manifest at $MANIFEST"
