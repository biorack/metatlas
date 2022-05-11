#!/bin/bash
#SBATCH --image=docker:doejgi/metatlas_shifter:latest
#SBATCH --constraint=haswell
#SBATCH --nodes=1
#SBATCH --licenses=cfs
#SBATCH --mail-type=ALL
#SBATCH --time=36:00:00

set -euo pipefail

log_dir="/global/cfs/projectdirs/m2650/jupyter_logs/slurm"

# make output notebook accessible for troubleshooting purposes
# want this to run even if we exit on a papermill error
trap \
  "{ cp "$OUT_FILE" "${log_dir}/${SLURM_JOB_ID}_$(basename "$OUT_FILE")" ; }" \
  EXIT

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export HDF5_USE_FILE_LOCKING=FALSE

log="${log_dir}/${SLURM_JOB_ID}.log"

output () {
  printf "%s\n" "$1" | tee --append "$log"
}

output "start time: $(date)"
output "user: $USER"
output "input file: $IN_FILE"
output "output file: $OUT_FILE"
output "parameters: $PARAMETERS"

# this creates the cache black uses and prevents some error messages
# doesn't need --entrypoint and is faster to leave it off
shifter /bin/bash -c 'black --quiet --check /metatlas_image_version &&
        papermill \
         /src/notebooks/reference/RT_Prediction.ipynb \
         - \
         -p model_only True \
	 --prepare-only \
	 -k papermill > /dev/null'

shifter --entrypoint /usr/local/bin/papermill -k "papermill" "$IN_FILE" "$OUT_FILE" $PARAMETERS 2>&1 | tee --append "$log"
