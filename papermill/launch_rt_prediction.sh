#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# fix threads at 8 even when requesting more CPUs,
# as the larger jobs need more memory per thread.
threads_to_use=8

usage() {
  >&2 echo "Usage:
  $(basename "$0") experiment_name analysis_number project_directory [-p notebook_parameter=value] [-y yaml_string]

     where:
        experiment_name:   experiment identifier
        analysis_number:   integer, use 0 for a new analysis
                           and increment if reworking one
        project_directory: output directory will be created within this directory
        -p:                optional notebook parameters, can have more than one
        -y:                optional notebook parameters in YAML string

  for more information see:
  https://github.com/biorack/metatlas/blob/main/docs/Targeted_Analysis.md
  "
  exit 128
}

validate_extra_parameters() {
  local array_name="$1[@]"
  local parameters=("${!array_name}")
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
  local kernel_path="${HOME}/.local/share/jupyter/kernels/metatlas-targeted"
  mkdir -p "$kernel_path"
  cp "${script_dir}/../docker/shifter.kernel.json" "${kernel_path}/kernel.json"
}

is_group_member() {
  local group_name=$1
  id -nG "$USER" | grep -qw "${group_name}"
}

is_C18_experiment() {
  local experiment_name="$1"
  [[ $experiment_name == *"_C18_"* ]]
}

get_num_cpus() {
  local experiment_name="$1"
  if is_C18_experiment "$experiment_name"; then
    echo "64"
  else
    echo "8"
  fi
}

get_slurm_account() {
  if is_group_member "m342"; then
    echo "m342"
  elif is_group_member "m2650"; then
    # could also use '--qos=flex' for lower cost and lower priority
    echo "m2650"
  else
    >&2 echo "WARNING: ${USER} is not a member of m342 or m2650. Attempting to use ${USER}'s default account."
    echo ""
  fi
}

get_slurm_queue() {
  local experiment_name="$1"
  if is_C18_experiment "$experiment_name"; then
    if is_group_member "m342"; then
      echo "genepool"
    else
      # could also use 'flex' for lower cost and lower priority
      echo "regular"
    fi
  else
    if is_group_member "m342"; then
      echo "genepool_shared"
    else
      # could also use 'flex' for lower cost and lower priority
      echo "shared"
    fi
  fi
}

declare -a positional_parameters=()
declare -a extra_parameters=()
while [ $OPTIND -le "$#" ]
do
  if getopts p:y: option; then
    case $option in
      p) extra_parameters+=("$OPTARG");;
      y) PARAMETERS="-y $OPTARG";;
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

IFS='_' read -ra TOKENS <<< "$exp"
proposal="${TOKENS[3]:-}"
exp_check_len="${TOKENS[8]:-}"

if [[ $exp_check_len == "" ]]; then
  >&2 echo "ERROR: experiment_name parameter is invalid. Must have 9 fields when split on '_'."
  usage
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
exp_dir="${project_dir}/$exp"
analysis_dir="${exp_dir}/${USER}${analysis_num}"

account="$(get_slurm_account)"
cpus_requested="$(get_num_cpus "$exp")"
queue="$(get_slurm_queue "$exp")"
IFS=$' ' flags="${account:+--account=$account} --qos=${queue} --cpus-per-task=${cpus_requested}"

IN_FILE="/src/notebooks/reference/RT_Prediction.ipynb"
OUT_FILE="${analysis_dir}/${proposal}_RT_Prediction_papermill.ipynb"
PARAMETERS+=" -p experiment $exp \
	      -p project_directory $project_dir \
	      -p max_cpus $threads_to_use \
	      -p analysis_number $analysis_num"
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
