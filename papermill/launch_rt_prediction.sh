#!/bin/bash

set -euf -o pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage $: experiment_name analysis_number project_directory"
fi

EXP="$1"
ANALYSIS_NUM="$2"
PROJECT_DIR="$3"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
EXP_DIR="${PROJECT_DIR}/$EXP"

IFS='_' read -ra TOKENS <<< "$EXP"
PROPOSAL="${TOKENS[0]}"

mkdir -p "$EXP_DIR"

export IN_FILE="${REPO_DIR}/notebooks/reference/RT_Prediction.ipynb"
export OUT_FILE="$/503256_RT_Prediction_papermill_12.ipynb"
export PARAMETERS="-p experiment $EXP -p metatlas_repo_path $REPO_DIR -p project_directory $PROJECT_DIR -p max_cpus 32 -p analysis_number $ANALYSIS_NUM"

sbatch -J "${PROPOSAL}_RT_Pred" "${REPO_DIR}/papermill/slurm_template.sh"
