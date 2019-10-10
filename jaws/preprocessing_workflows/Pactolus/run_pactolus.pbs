#!/bin/bash -l
#PBS -A m1541
#PBS -N gnps_mxd5
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l mppwidth=120
#PBS -l mppnppn=24
#PBS -q regular


export INSTALL_DIR=/scratch1/scratchdirs/curt/  # Location where the packages/ dir is located
export INCHI_FILE=/scratch1/scratchdirs/curt/gnps_molecule_inchis.txt
export ISOTOPE_FILE=/scratch1/scratchdirs/curt/packages/meta-iq/pactolus/data/max_abundance_isotopes.csv
export OUTPUT_DIR=/scratch1/scratchdirs/curt/gen_tree_result_enzo_abc
export ERROR_LOG=/scratch1/scratchdirs/curt/gen_tree_error_log.txt
export MAX_DEPTH=5
export TASK_PER_NODE=24  # -n of aprun is the total number of MPI tasks, always equal to mppnppn
export NUM_TASKS=120      # -N of aprun is the number of tasks per node. Usually we use the maximum of 24. 
# -n and -N of aprun must be set according to mppwidth and mppnppn. If mppnppn is 24 then TASK_PER_NODE=24
# and NUM_TASKS=mppwidth. Otherwise, if less than 24 tasks are used per node, the mppwidth must be increased
# accordingly.  NUM_TASKS = mppwdith / 24 * mppnppn


cd $PBS_O_WORKDIR
module load python/2.7-anaconda

export PYTHONPATH=$INSTALL_DIR/packages/bastet:$PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR/anaconda/lib/python2.7/site-packages
export LD_LIBRARY_PATH=/usr/common/usg/python/2.7-anaconda/lib:$LD_LIBRARY_PATH

aprun -n $NUM_TASKS -N $TASK_PER_NODE python $INSTALL_DIR/packages/meta-iq/pactolus/pactolus/generate_frag_dag.py --inchi_file $INCHI_FILE --isotope_file $ISOTOPE_FILE --output_dir $OUTPUT_DIR --error_log $ERROR_LOG --max_depth $MAX_DEPTH



