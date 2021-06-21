#!/bin/bash

set -euf -o pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 metatlas_repo_dir base_output_dir experiment_id"
  exit 1
fi

REPO_DIR="$1"
OUT_DIR="$2"
EXP="$3"

function install_kernel () {
  local REPO_DIR="$1"
  local SOURCE="${REPO_DIR}/notebooks/kernels/metatlas-targeted.kernel.json"
  local DEST="${HOME}/.local/share/jupyter/kernels/metatlas-targeted/kernel.json"
  if [[ ! -f "$DEST" ]]; then
    mkdir -p $(dirname "$DEST")
    cp "$SOURCE" "$DEST"
  fi
}

function validate_data_dir () {
  local EXP="$1"
DATA_DIR="/project/projectdirs/metatlas/raw_data/akuftin/${EXP}"

if [ ! -d "${DATA_DIR}" ]; then
  echo "ERROR: could not find data directory ${DATA_DIR}." >&2
  exit 2
fi

function populate_experiment_dir () {
IFS='_' read -ra EXP_ARRAY <<< "$EXP"
NOTEBOOK_BASE="${EXP_ARRAY[3]}_${EXP_ARRAY[4]}"

mkdir -p "${OUT_DIR}/${EXP}"
cp "${REPO_DIR}/notebooks/reference/Workflow_Notebook_VS_Auto_RT_Predict_V2.ipynb" "${OUT_DIR}/${EXP}/${NOTEBOOK_BASE}_RT_Predict.ipynb"
cp "${REPO_DIR}/notebooks/reference/Targeted.ipynb" "${OUT_DIR}/${EXP}/${NOTEBOOK_BASE}.ipynb"
