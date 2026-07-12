#!/usr/bin/env bash
# Prepare a small, AADR-style canid genotype panel from the indexed 722-genome VCF.
#
# This does NOT download the 309 GB source VCF. bcftools uses the public tabix index and
# HTTP byte-range requests to fetch only selected CanineHD marker positions for:
#   Wolf25          red wolf, 28.5x in the source metadata
#   Wolf26          red wolf,  7.1x in the source metadata
#   GoldenJackal01  golden jackal, 5.5x in the source metadata
#
# Run from the repository root under Linux/WSL2:
#   bash scripts/prepare_redwolf_jackal_aadr.sh
#
# For an even smaller/faster smoke test:
#   TARGET_SITES=10000 bash scripts/prepare_redwolf_jackal_aadr.sh
set -euo pipefail

SOURCE_VCF_URL="${SOURCE_VCF_URL:-https://research.nhgri.nih.gov/dog_genome/downloads/datasets/WGS/722g.990.SNP.INDEL.chrAll.vcf.gz}"
PANEL_URL="${PANEL_URL:-https://download.cncb.ac.cn/dogsd/dog10k/variations/CFA31_IlluminaHD.vcf.gz}"
PANEL_MD5="${PANEL_MD5:-4a87088b17631fb6237210044ac099a6}"
TARGET_SITES="${TARGET_SITES:-25000}"
OUT_DIR="${OUT_DIR:-data/redwolf_jackal_aadr}"

SAMPLES=(Wolf25 Wolf26 GoldenJackal01)
SAMPLE_CSV="$(IFS=,; echo "${SAMPLES[*]}")"
PANEL_VCF="$OUT_DIR/CFA31_IlluminaHD.vcf.gz"
ALL_POSITIONS="$OUT_DIR/caninehd_autosomal_positions.tsv"
TARGET_POSITIONS="$OUT_DIR/caninehd_even_${TARGET_SITES}.regions.tsv"
OUTPUT_VCF="$OUT_DIR/redwolf_jackal_caninehd${TARGET_SITES}.vcf.gz"
MANIFEST="$OUT_DIR/redwolf_jackal_caninehd${TARGET_SITES}.manifest.json"

for tool in bcftools curl md5sum awk python; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required tool '$tool' is not on PATH" >&2
    echo "Install/activate the repository conda environment first." >&2
    exit 1
  fi
done

if ! [[ "$TARGET_SITES" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: TARGET_SITES must be a positive integer; got '$TARGET_SITES'" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

echo ">> checking the remote 722-genome VCF and tabix index"
curl --fail --silent --show-error --head "$SOURCE_VCF_URL" >/dev/null
curl --fail --silent --show-error --head "${SOURCE_VCF_URL}.tbi" >/dev/null

echo ">> downloading the 1.5 MB CanineHD marker panel"
if [[ ! -s "$PANEL_VCF" ]]; then
  curl --fail --location --retry 3 --output "$PANEL_VCF" "$PANEL_URL"
fi
echo "$PANEL_MD5  $PANEL_VCF" | md5sum --check --status || {
  echo "ERROR: CanineHD panel checksum did not match $PANEL_MD5" >&2
  exit 1
}

echo ">> confirming the requested samples are present in the remote VCF header"
REMOTE_SAMPLES="$(bcftools query -l "$SOURCE_VCF_URL")"
for sample in "${SAMPLES[@]}"; do
  if ! grep -Fxq "$sample" <<<"$REMOTE_SAMPLES"; then
    echo "ERROR: sample '$sample' is absent from the 722-genome VCF" >&2
    exit 1
  fi
done

echo ">> building an evenly spaced autosomal target panel"
bcftools query -f '%CHROM\t%POS\t%POS\n' "$PANEL_VCF" \
  | awk '$1 ~ /^[0-9]+$/ && $1 >= 1 && $1 <= 38' \
  > "$ALL_POSITIONS"

N_AVAILABLE="$(wc -l < "$ALL_POSITIONS" | tr -d ' ')"
if [[ "$N_AVAILABLE" -eq 0 ]]; then
  echo "ERROR: no CanFam3.1 autosomal positions were read from the marker panel" >&2
  exit 1
fi
STRIDE=$(( (N_AVAILABLE + TARGET_SITES - 1) / TARGET_SITES ))
awk -v stride="$STRIDE" '(NR - 1) % stride == 0' "$ALL_POSITIONS" \
  > "$TARGET_POSITIONS"
N_TARGETED="$(wc -l < "$TARGET_POSITIONS" | tr -d ' ')"

echo ">> extracting $N_TARGETED indexed sites for $SAMPLE_CSV"
echo "   source VCF is remote; only indexed byte ranges will be transferred"
bcftools view \
  --regions-file "$TARGET_POSITIONS" \
  --samples "$SAMPLE_CSV" \
  --min-alleles 2 \
  --max-alleles 2 \
  --types snps \
  --output-type u \
  "$SOURCE_VCF_URL" \
  | bcftools view \
      --exclude 'GT="mis"' \
      --output-type z \
      --threads 2 \
      --output "$OUTPUT_VCF"

bcftools index --force --tbi "$OUTPUT_VCF"
N_RETAINED="$(bcftools index --nrecords "$OUTPUT_VCF")"
if [[ "$N_RETAINED" -lt 1000 ]]; then
  echo "ERROR: only $N_RETAINED complete biallelic SNPs were retained; expected >=1000" >&2
  exit 1
fi

OUTPUT_SHA256="$(sha256sum "$OUTPUT_VCF" | awk '{print $1}')"
export SOURCE_VCF_URL PANEL_URL PANEL_MD5 TARGET_SITES N_AVAILABLE N_TARGETED
export N_RETAINED OUTPUT_VCF OUTPUT_SHA256 MANIFEST
python - <<'PY'
import json
import os
from datetime import datetime, timezone

manifest = {
    "dataset_id": "redwolf_jackal_722g_caninehd",
    "description": "AADR-style reduced CanFam3.1 SNP panel; not raw genomes.",
    "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "source": {
        "bioproject": "PRJNA448733",
        "analysis_accession": "SRZ189891",
        "vcf_url": os.environ["SOURCE_VCF_URL"],
        "reference": "CanFam3.1",
    },
    "marker_panel": {
        "name": "CFA31_IlluminaHD",
        "url": os.environ["PANEL_URL"],
        "md5": os.environ["PANEL_MD5"],
        "available_autosomal_sites": int(os.environ["N_AVAILABLE"]),
        "requested_max_sites": int(os.environ["TARGET_SITES"]),
        "targeted_sites": int(os.environ["N_TARGETED"]),
    },
    "samples": [
        {"sample_id": "Wolf25", "taxon": "red_wolf", "biosample": "SAMN02921317"},
        {"sample_id": "Wolf26", "taxon": "red_wolf", "biosample": "SAMN02921318"},
        {"sample_id": "GoldenJackal01", "taxon": "golden_jackal", "biosample": "SAMN03366713"},
    ],
    "filters": {
        "chromosomes": "CanFam3.1 autosomes 1-38",
        "variant_type": "biallelic SNP",
        "missing_genotypes": "excluded if missing in any selected sample",
    },
    "output": {
        "path": os.environ["OUTPUT_VCF"],
        "retained_sites": int(os.environ["N_RETAINED"]),
        "sha256": os.environ["OUTPUT_SHA256"],
    },
}
with open(os.environ["MANIFEST"], "w", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2)
    handle.write("\n")
PY

echo
echo ">> dataset ready: $OUTPUT_VCF"
echo ">> retained SNPs: $N_RETAINED"
echo ">> provenance: $MANIFEST"
echo ">> next command (not run by this script):"
echo "   canidae run -c configs/examples/redwolf_jackal_aadr.yaml -c configs/profiles/laptop.yaml"

