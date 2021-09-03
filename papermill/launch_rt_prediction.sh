#!/bin/bash
set -euf -o pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage $0: experiment_name analysis_number project_directory"
    exit 0
fi

EXP="$1"
ANALYSIS_NUM="$2"
PROJECT_DIR="$3"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
EXP_DIR="${PROJECT_DIR}/$EXP"
ANALYSIS_DIR="${EXP_DIR}/${USER}${ANALYSIS_NUM}"
KERNEL_SOURCE="${SCRIPT_DIR}/notebooks/kernels/metatlas-targeted.kernel.json"
KERNEL_DESTINATION="${HOME}/.local/share/jupyter/kernels/metatlas-targeted/kernel.json"

IFS='_' read -ra TOKENS <<< "$EXP"
PROPOSAL="${TOKENS[3]}"

export IN_FILE="${REPO_DIR}/notebooks/reference/RT_Prediction.ipynb"
export OUT_FILE="${ANALYSIS_DIR}/${PROPOSAL}_RT_Prediction_papermill.ipynb"
export PARAMETERS="-p experiment $EXP -p metatlas_repo_path $REPO_DIR -p project_directory $PROJECT_DIR -p max_cpus 32 -p analysis_number $ANALYSIS_NUM"

mkdir -p "${HOME}/.local/share/jupyter/kernels/metatlas-targeted"
cp "$KERNEL_SOURCE" "$KERNEL_DESTINATION"

mkdir -p "$ANALYSIS_DIR"
sbatch -J "${PROPOSAL}_RT_Pred" "${REPO_DIR}/papermill/slurm_template.sh"
