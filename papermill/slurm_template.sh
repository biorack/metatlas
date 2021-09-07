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

CONDA_DIR="$(dirname "$(dirname "$(grep 'metatlas-targeted' ../notebooks/kernels/metatlas-targeted.kernel.json | cut -d\" -f 2)")")"
date
echo "input file: $IN_FILE"
echo "output file: $OUT_FILE"
module load python/3.8-anaconda-2020.11
eval "$(conda shell.bash hook)"
conda activate "$CONDA_DIR"

srun -n 1 -c 64 --cpu_bind=cores papermill "$IN_FILE" "$OUT_FILE" $PARAMETERS
