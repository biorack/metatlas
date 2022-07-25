#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# fix threads at 8 even when requesting more CPUs,
# as the larger jobs need more memory per thread.
threads_to_use=8

rclone='/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone'

trap "exit 1" TERM
export TOP_PID=$$

die() {
  kill -s TERM $TOP_PID
}

usage() {
  >&2 echo "Usage:
  $(basename "$0") workflow_name experiment_name [rt_predict_number] [project_directory] [-p notebook_parameter=value] [-y yaml_string]

     where:
        workflow_name:     name associated with a workflow definition in the configuration file
        experiment_name:   experiment identifier
        rt_predict_number: integer, use 0 the first time generating an RT correction for an experiment
	                   and increment if re-generating an RT correction (default: 0)
	project_directory: output directory will be created within this directory (default: $HOME/metabolomics_data)
        -p:                optional notebook parameters, can use multiple times
        -y:                optional notebook parameters in YAML or JSON string

  for more information see:
  https://github.com/biorack/metatlas/blob/main/docs/Targeted_Analysis.md
  "
  die
}

validate_extra_parameters() {
  local array_name="$1[@]"
  local parameters=("${!array_name}")
  for i in "${parameters[@]}"
  do
    if [ "${i:0:1}" == "=" ]; then
      >&2 echo "ERROR: invalid usage of -p in '-p ${i}'"
      die
    fi
    IFS='=' read -r -a arr <<< "$i"
    if [ ${#arr[@]} -ne 2 ]; then
      >&2 echo "ERROR: invalid usage of -p in '-p ${i}'"
      die
    fi
  done
}

is_valid_yaml() {
  echo "$1" | shifter --module=none --clearenv --image=doejgi/metatlas_shifter:latest \
	  /src/metatlas/scripts/yaml_validation.py
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

is_perlmutter() {
  [ "$NERSC_HOST" = "perlmutter" ]
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
  if is_group_member "gtrnd" && ! is_perlmutter; then
    echo "gtrnd"
  elif is_group_member "m2650"; then
    echo "m2650"
  else
    >&2 echo "WARNING: ${USER} is not a member of gtrnd or m2650. Attempting to use ${USER}'s default account."
    echo ""
  fi
}

get_slurm_time() {
  if is_perlmutter; then
    echo "12:00:00"
  else
    echo "36:00:00"
  fi
}

get_slurm_constraint() {
  if is_perlmutter; then
    echo "cpu"
  else
    echo "haswell"
  fi
}

get_slurm_queue() {
  local experiment_name="$1"
  if is_perlmutter; then
    echo "regular"
  else  # cori
    if is_C18_experiment "$experiment_name"; then
      if is_group_member "gtrnd"; then
        echo "genepool"
      else
        # could also use 'flex' for lower cost and lower priority
        echo "regular"
      fi
    else
      if is_group_member "gtrnd"; then
        echo "genepool_shared"
      else
        # could also use 'flex' for lower cost and lower priority
        echo "shared"
      fi
    fi
  fi
}

get_rclone_remote() {
  # this assumes you only have one google account that is set up within rclone...
  remote="$("$rclone" listremotes --long 2> /dev/null | grep "drive$" | head -1 | cut -d' ' -f1)"
  if [ -z "$remote" ]; then
    >&2 echo "ERROR: rclone has not been configured to access Google Drive. Follow instructions at:"
    >&2 echo "       https://github.com/biorack/metatlas/blob/main/docs/Targeted_Analysis.md#rclone-configuration"
    die
  fi
  echo "$remote"
}

check_yaml_is_valid() {
  if ! is_valid_yaml "$1"; then
    >&2 echo "ERROR: invalid YAML or JSON for -y value."
    die
  fi
}

check_gdrive_authorization() {
  if !  "$rclone" lsf "$(get_rclone_remote)" > /dev/null 2>&1; then
    >&2 echo "ERROR: rclone authoriation to Google Drive has expired. Please run:"
    >&2 echo "       ${rclone} config reconnect $(get_rclone_remote)"
    die
  fi
}

check_analysis_dir_does_not_exist() {
  if [ -d "$1" ]; then
    >&2 echo "ERROR: Output directory already exists. Not overwriting:"
    >&2 echo "       ${1}"
    >&2 echo "       Consider incrementing rt_predict_number."
    die
  fi
}

check_project_dir_does_exist() {
  if ! [ -d "$1" ]; then
    >&2 echo "ERROR: project_directory '${1}' does not exist."
    >&2 echo "       Please run 'mkdir \"${1}\"' first."
    die
  fi
}

check_exp_id_has_atleast_9_fields() {
  # inputs: the 9th field (1-indexed) of the experiment_name split on '_'
  if [[ $1 == "" ]]; then
    >&2 echo "ERROR: experiment_name parameter is invalid. Must have 9 fields when split on '_'."
    die
  fi
}

check_not_in_commom_software_filesystem() {
  if [[ $(pwd) == /global/common/software/* ]]; then
    >&2 echo ""
    >&2 echo "ERROR: You cannot submit a SLURM job from a directory that will be"
    >&2 echo "read only from a job execution node, such as any directory under"
    >&2 echo "/global/common/software"
    >&2 echo "Please change to a different directory and try again."
    >&2 echo "No SLURM jobs have been submitted."
    >&2 echo ""
    die
  fi
}

check_rt_predict_number_is_non_neg_int() {
  re='^[0-9]+$'
  if ! [[ $1 =~ $re ]] ; then
    >&2 echo ""
    >&2 echo "ERROR: rt_predict_number must be a non-negative integer"
    >&2 echo ""
    die
  fi
}

get_script_dir() {
  # adapted from https://stackoverflow.com/questions/59895/
  SOURCE=${BASH_SOURCE[0]}
  while [ -L "$SOURCE" ]; do
    DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
    SOURCE=$(readlink "$SOURCE")
    [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE
  done
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  echo "$DIR"
}

YAML_BASE64="$(echo "{}" | base64 --wrap=0)"
declare -a positional_parameters=()
declare -a extra_parameters=()
while [ $OPTIND -le "$#" ]
do
  if getopts p:y: option; then
    case $option in
      p) extra_parameters+=("$OPTARG");;
      y) YAML_BASE64="$(echo "${OPTARG}" | base64 --wrap=0)";;
      \?) usage;;
    esac
  else
    positional_parameters+=("${!OPTIND}")
    ((OPTIND++))
  fi
done

if [  ${#positional_parameters[@]} -lt 2 ]; then
  >&2 echo "ERROR: one of workflow_name or experiment_name was not supplied."
  >&2 echo ""
  usage
fi

if [  ${#positional_parameters[@]} -gt 4 ]; then
  >&2 echo "ERROR: too many positional parameters supplied."
  >&2 echo ""
  usage
fi

if [  ${#extra_parameters[@]} -ne 0 ]; then
  validate_extra_parameters extra_parameters  # pass extra_parameters by name
fi

workflow_name="${positional_parameters[0]}"
exp="${positional_parameters[1]}"
rt_predict_num="${positional_parameters[2]:-0}"
project_dir="${positional_parameters[3]:-$HOME/metabolomics_data}"

script_dir="$(get_script_dir)"
exp_dir="${project_dir}/$exp"
analysis_dir="${exp_dir}/${USER}_${rt_predict_num}_0"

IFS='_' read -ra TOKENS <<< "$exp"
proposal="${TOKENS[3]:-}"
exp_check_len="${TOKENS[8]:-}"

check_exp_id_has_atleast_9_fields "$exp_check_len"
check_analysis_dir_does_not_exist "$analysis_dir"
check_project_dir_does_exist "$project_dir"
check_rt_predict_number_is_non_neg_int "$rt_predict_num"
check_yaml_is_valid "$(echo "${YAML_BASE64:-}" | base64 --decode)"
check_gdrive_authorization
check_not_in_commom_software_filesystem

account="$(get_slurm_account)"
cpus_requested="$(get_num_cpus "$exp")"
queue="$(get_slurm_queue "$exp")"
constraint="$(get_slurm_constraint)"
time="$(get_slurm_time)"
IFS=$' ' flags="${account:+--account=$account} --qos=${queue} --cpus-per-task=${cpus_requested} --constraint=${constraint} --time=${time}"

IN_FILE="/src/notebooks/reference/RT_Alignment.ipynb"
OUT_FILE="${analysis_dir}/${proposal}_${workflow_name}_RT_Alignment_SLURM.ipynb"

PARAMETERS+=" -p experiment $exp \
	      -p workflow_name ${workflow_name} \
	      -p project_directory $project_dir \
	      -p max_cpus $threads_to_use \
	      -p rt_predict_number $rt_predict_num \
	      -p config_file_name /global/cfs/cdirs/m2650/targeted_analysis/metatlas_config.yaml"
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
export YAML_BASE64

install_jupyter_kernel

mkdir -p "$analysis_dir"
# shellcheck disable=SC2086
sbatch $flags -J "${proposal}_RT_Align" "${script_dir}/slurm_template.sh"
