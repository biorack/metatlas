#!/bin/bash
set -euf -o pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage $0: experiment_name analysis_number project_directory"
    exit 0
fi

if [[ $(pwd) == /global/common/software/* ]]; then
   >&2 echo ""
   >&2 echo "ERROR: You cannot submit a SLURM job from a directory that will be"
   >&2 echo "read only from a job execution node, such as any directory under"
   >&2 echo "/global/common/software"
   >&2 echo "Please change to a different directory and try again."
   >&2 echo "No SLURM jobs have been submitted."
   >&2 echo ""
   exit 1
fi

EXP="$1"
ANALYSIS_NUM="$2"
PROJECT_DIR="$3"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
EXP_DIR="${PROJECT_DIR}/$EXP"
ANALYSIS_DIR="${EXP_DIR}/${USER}${ANALYSIS_NUM}"

IFS='_' read -ra TOKENS <<< "$EXP"
PROPOSAL="${TOKENS[3]}"

export IN_FILE="/src/notebooks/reference/RT_Prediction.ipynb"
export OUT_FILE="${ANALYSIS_DIR}/${PROPOSAL}_RT_Prediction_papermill.ipynb"
export PARAMETERS="-p experiment $EXP -p project_directory $PROJECT_DIR -p max_cpus 32 -p analysis_number $ANALYSIS_NUM"

KERNEL_PATH="${HOME}/.local/share/jupyter/kernels/metatlas-targeted/kernel.json"
cp "${SCRIPT_DIR}/kernel.json" "$KERNEL_PATH"

mkdir -p "$ANALYSIS_DIR"
sbatch -J "${PROPOSAL}_RT_Pred" "${SCRIPT_DIR}/slurm_template.sh"
