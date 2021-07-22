#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH --account=gtrnd
#SBATCH --qos=genepool
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

date
echo "input file: $IN_FILE"
echo "output file: $OUT_FILE"
eval "$(conda shell.bash hook)"
conda activate /global/common/software/m2650/metatlas-targeted-2021-07-16

srun -n 1 -c 64 --cpu_bind=cores papermill "$IN_FILE" "$OUT_FILE" $PARAMETERS

