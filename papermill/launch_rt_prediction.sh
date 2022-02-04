#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

max_cpus=8

usage() { 
>&2 echo "Usage:
$(basename "$0") experiment_name analysis_number project_directory [-p notebook_parameter=value]

   where:
      experiment_name:   experiment identifier 
      analysis_number:   integer, use 0 for a new analysis
                         and increment if reworking one
      project_directory: output directory will be created within this directory
      -p:                optional notebook parameters, can have more than one 

for more information see:
https://github.com/biorack/metatlas/blob/main/docs/Targeted_Analysis.md
"
exit 128
}

validate_extra_parameters() {
  array_name="$1[@]"
  parameters=("${!array_name}")
  for i in "${parameters[@]}"
  do
    if [ "${i:0:1}" == "=" ]; then
      usage
    fi
    IFS='=' read -r -a arr <<< "$i"
    if [ ${#arr[@]} -ne 2 ]; then
      usage
    fi
  done
}

install_jupyter_kernel() {
  kernel_path="${HOME}/.local/share/jupyter/kernels/metatlas-targeted"
  mkdir -p "$kernel_path"
  cp "${script_dir}/../docker/shifter.kernel.json" "${kernel_path}/kernel.json"
}

declare -a positional_parameters=()
declare -a extra_parameters=()
while [ $OPTIND -le "$#" ]
do
  if getopts p: option; then
    case $option in
      p) extra_parameters+=("$OPTARG");;
      \?) usage;;
    esac
  else
    positional_parameters+=("${!OPTIND}")
    ((OPTIND++))
  fi
done

if [  ${#positional_parameters[@]} -ne 3 ]; then
  usage
fi

if [  ${#extra_parameters[@]} -ne 0 ]; then
  validate_extra_parameters extra_parameters
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

exp="${positional_parameters[0]}"
analysis_num="${positional_parameters[1]}"
project_dir="${positional_parameters[2]}"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
exp_dir="${project_dir}/$exp"
analysis_dir="${exp_dir}/${USER}${analysis_num}"

IFS='_' read -ra TOKENS <<< "$exp"
proposal="${TOKENS[3]:-}"
exp_check_len="${TOKENS[8]:-}"

if [[ $exp_check_len == "" ]]; then
  >&2 echo "ERROR: experiment_name parameter is invalid. Must have 9 fields when split on '_'."
  usage
fi

if id -nG "$USER" | grep -qw "gtrnd"; then
  flags="--account=gtrnd --qos=genepool_shared"
elif id -nG "$USER" | grep -qw "m2650"; then
  # could also use '--qos=flex' for lower cost and lower priority
  flags="--account=m2650 --qos=shared"
else
  echo "WARNING: ${USER} is not a member of gtrnd or m2650. Attempting to use ${USER}'s default account."
  flags="--qos=shared"
fi
flags+=" --cpus-per-task=${max_cpus}"

IN_FILE="/src/notebooks/reference/RT_Prediction.ipynb"
OUT_FILE="${analysis_dir}/${proposal}_RT_Prediction_papermill.ipynb"
PARAMETERS="-p experiment $exp -p project_directory $project_dir -p max_cpus $max_cpus -p analysis_number $analysis_num"
if [  ${#extra_parameters[@]} -ne 0 ]; then
  for i in "${extra_parameters[@]}"
  do
    IFS='=' read -r -a arr <<< "$i"
    PARAMETERS+=" -p ${arr[0]} ${arr[1]}"
  done
fi
export IN_FILE
export OUT_FILE
export PARAMETERS

install_jupyter_kernel

mkdir -p "$analysis_dir"
sbatch $flags -J "${proposal}_RT_Pred" "${script_dir}/slurm_template.sh"
