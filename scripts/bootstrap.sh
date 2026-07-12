#!/usr/bin/env bash
# One-shot setup for the read-processing pipeline inside a fresh Linux/WSL2 shell.
# Assumes conda/mamba is on PATH (install Miniforge first — see docs/REAL_DATA_REDWOLF.md).
# Run from the repository root:  bash scripts/bootstrap.sh
set -euo pipefail

ENV_NAME="${ENV_NAME:-canis}"
MAMBA="$(command -v mamba || command -v conda)"
if [ -z "$MAMBA" ]; then
  echo "ERROR: install Miniforge first (conda/mamba not found). See docs/REAL_DATA_REDWOLF.md"
  exit 1
fi

echo ">> creating conda env '$ENV_NAME' (bwa-mem2, samtools, GATK, ... + python)"
"$MAMBA" env create -n "$ENV_NAME" -f environment.yml || \
  "$MAMBA" env update -n "$ENV_NAME" -f environment.yml

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo ">> installing the canidae package (analysis extras)"
pip install -e ".[analysis]"

echo
echo "Environment ready. Next:"
echo "  conda activate $ENV_NAME"
echo "  bash scripts/fetch_reference.sh"
echo "  SUBSAMPLE=0.1 bash scripts/fetch_redwolf_reads.sh   # laptop-scale proof-of-concept"
echo "  canidae run -c configs/examples/redwolf_real.yaml -c configs/profiles/laptop.yaml"
