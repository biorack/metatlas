#!/bin/bash
#SBATCH --image=docker:doejgi/metatlas_shifter:latest
#SBATCH --nodes=1
#SBATCH --licenses=cfs
#SBATCH --mail-type=ALL

set -euo pipefail

# PAPERMILL_EXECUTION is set to True so we can determine if notebooks
# are being run in papermill or interactively
shifter_flags="--module=none --clearenv --env=PAPERMILL_EXECUTION=True"

log_dir="/global/cfs/projectdirs/m2650/jupyter_logs/slurm"

# make output notebook accessible for troubleshooting purposes
# want this to run even if we exit on a papermill error
# shellcheck disable=SC2027,SC2086,SC2046
trap \
  "{ cp "$OUT_FILE" "${log_dir}/${SLURM_JOB_ID}_$(basename "$OUT_FILE")" ; }" \
  EXIT

log="${log_dir}/${SLURM_JOB_ID}.log"

output () {
  printf "%s\n" "$1" | tee --append "$log"
}

output "start time: $(date)"
output "user: $USER"
output "input file: $IN_FILE"
output "output file: $OUT_FILE"
output "parameters: $PARAMETERS"
output "yaml parameters: $(echo "$YAML_BASE64" | base64 --decode)"

# this creates the cache black uses and prevents some error messages
# doesn't need --entrypoint and is faster to leave it off
# shellcheck disable=SC2086
shifter $shifter_flags /bin/bash -c \
  'set -euo pipefail && \
   black --quiet --check /metatlas_image_version && \
   papermill \
     /src/notebooks/reference/RT_Alignment.ipynb \
     - \
     -p model_only True \
     --prepare-only \
     -k papermill > /dev/null'

# shellcheck disable=SC2086
shifter --entrypoint $shifter_flags \
  /usr/local/bin/papermill \
  -k "papermill" \
  "$IN_FILE" \
  "$OUT_FILE" \
  --parameters_base64 "$YAML_BASE64" \
  $PARAMETERS 2>&1 | tee --append "$log"
