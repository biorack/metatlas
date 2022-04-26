#!/bin/bash
#SBATCH --image=docker:wjhjgi/metatlas_shifter:latest
#SBATCH --constraint=haswell
#SBATCH --nodes=1
#SBATCH --licenses=cfs
#SBATCH --mail-type=ALL
#SBATCH --time=36:00:00

set -euo pipefail

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export HDF5_USE_FILE_LOCKING=FALSE

log_dir="/global/cfs/projectdirs/m2650/jupyter_logs/slurm"
log="${log_dir}/${SLURM_JOB_ID}.log"

output () {
  printf "%s\n" "$1" | tee --append "$log"
}

output "start time: $(date)"
output "user: $USER"
output "input file: $IN_FILE"
output "output file: $OUT_FILE"
output "parameters: $PARAMETERS"

shifter --entrypoint /usr/local/bin/papermill -k "papermill" "$IN_FILE" "$OUT_FILE" $PARAMETERS 2>&1 | tee --append "$log"

# make output notebook accessible for troubleshooting purposes
nb_file_name="$(basename "${OUT_FILE}")"
cp "$OUT_FILE" "${log_dir}/${SLURM_JOB_ID}_${nb_file_name}"
