#!/bin/bash
#SBATCH --image=docker:doejgi/metatlas_shifter:latest
#SBATCH --nodes=1
#SBATCH --licenses=cfs
#SBATCH --mail-type=ALL

set -euo pipefail

shifter_flags="--module=none --clearenv \
   --env=SOURCE_CODE_VERSION_ID=${SOURCE_CODE_VERSION_ID:-2.x}"

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
# shellcheck disable=SC2086,SC2016
shifter --entrypoint $shifter_flags /bin/bash -c \
  'black --quiet --check /metatlas_image_version && \
   papermill \
     "${METATLAS_WORKING_SOURCE_DIR}/notebooks/reference/RT_Prediction.ipynb" \
     - \
     -p model_only True \
     --prepare-only \
     -k papermill > /dev/null'

# IN_FILE contains an environmental variable that should be evaluated
# in the context of the new shifter shell.
# shellcheck disable=SC2086,SC2016,SC2116
shifter --entrypoint $shifter_flags \
        "--env=PARAMETERS=$PARAMETERS" \
	"--env=YAML_BASE64=$YAML_BASE64" \
	"--env=log=$log" \
	/bin/bash -c \
          '/usr/local/bin/papermill \
             -k "papermill" \
	     "'"$IN_FILE"'" \
             "$OUT_FILE" \
             --parameters_base64 "$YAML_BASE64" \
             $PARAMETERS 2>&1 | tee --append "$log"'
