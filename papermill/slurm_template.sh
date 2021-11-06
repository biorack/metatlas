#!/bin/bash
#SBATCH --image=docker:wjhjgi/metatlas_shifter:latest
#SBATCH --constraint=haswell
#SBATCH --account=gtrnd
#SBATCH --qos=genepool
#SBATCH --licenses=cfs
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

date
echo "input file: $IN_FILE"
echo "output file: $OUT_FILE"

shifter --entrypoint /usr/local/bin/papermill -k "papermill" "$IN_FILE" "$OUT_FILE" $PARAMETERS
