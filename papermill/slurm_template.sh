#!/bin/bash
#SBATCH --image=docker:wjhjgi/metatlas_shifter:latest
#SBATCH --constraint=haswell
#SBATCH --nodes=1
#SBATCH --licenses=cfs
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export HDF5_USE_FILE_LOCKING=FALSE

LOG="/global/cfs/projectdirs/m2650/jupyter_logs/${USER}.log"

set -o pipefail

date
echo "input file: $IN_FILE"
echo "output file: $OUT_FILE"
echo "parameters: $PARAMETERS"

(shifter --entrypoint /usr/local/bin/papermill -k "papermill" "$IN_FILE" "$OUT_FILE" $PARAMETERS) 2>&1 | tee --append "$LOG"
