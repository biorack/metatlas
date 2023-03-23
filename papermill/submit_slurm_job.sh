#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

rclone='/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone'
raw_data='/global/cfs/cdirs/metatlas/raw_data'

# fix threads with the notebook even when requesting more CPUs,
# as the larger jobs need more memory per thread.
threads_to_use=8

# system definitions
cori_cpus=64 # hyperthreads per Cori Haswell node
cori_mem=128 # GB per Cori Haswell node
perlmutter_cpus=256 # 64 cores per chip, 2 chips per node for CPU nodes, 2 threads per core
perlmutter_mem=512 # GB per node

# default memory request for SLURM job
memory=48 # GB

trap "exit 1" TERM
export TOP_PID=$$

die() {
  kill -s TERM $TOP_PID
}

usage() {
  >&2 echo "Usage:
  $(basename "$0") workflow_name experiment_name [rt_alignment_number] [project_directory] [-m memory] [-p notebook_parameter=value] [-y yaml_string]

     where:
        workflow_name:     name associated with a workflow definition in the configuration file
        experiment_name:   experiment identifier
        rt_alignment_number:  integer, use 0 the first time generating an RT alignment for an experiment
	                   and increment if re-generating an RT alignment (default: 0)
	project_directory: output directory will be created within this directory (default: $HOME/metabolomics_data)
	-m:                memory in GB to request for the SLURM job (default: ${memory})
        -p:                optional notebook parameters, can use multiple times
        -y:                optional notebook parameters in YAML or JSON string

  for more information see:
  https://github.com/biorack/metatlas/blob/main/docs/Targeted_Analysis.md
  "
  die
}

contains_parameter_name() {
  local array_name="$1[@]"
  local parameter_name="$2"
  local parameters=("${!array_name}")
  for i in "${parameters[@]}"
  do
    IFS='=' read -r -a arr <<< "$i"
    if [ ${#arr[0]} = "${parameter_name}" ]; then
      return 0
    fi
  done
  return 1
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

ceiling_divide() {
  echo "($1 + $2 - 1)/$2" | bc
}

total_cpus() {
  if is_perlmutter; then
    echo "$perlmutter_cpus"
  else
    echo "$cori_cpus"
  fi
}

total_mem() {
  # all values are in GB
  if is_perlmutter; then
    echo "$perlmutter_mem"
  else
    echo "$cori_mem"
  fi
}

get_num_cpus() {
  local needed_cpus mem
  mem="$1"
  needed_cpus="$(ceiling_divide "$(( "$(total_cpus)" * mem ))" "$(total_mem)")"
  if [ "$needed_cpus" -gt "$(( "$(total_cpus)" / 2 ))" ]; then
    # need more than half a node, can't use shared queue, so get full node
    total_cpus
  else
    echo "$needed_cpus"
  fi
}

needs_full_node() {
  # takes memory in GB as an input
  [ "$(get_num_cpus "$1")" -eq "$(total_cpus)" ]
}

get_slurm_account() {
  if is_group_member "gtrnd" && ! is_perlmutter; then
    echo "gtrnd"
  elif is_group_member "m342"; then
    echo "m342"
  elif is_group_member "m2650"; then
    echo "m2650"
  else
    >&2 echo "WARNING: ${USER} is not a member of gtrnd, m342, or m2650."
    >&2 echo "WARNING: Attempting to use ${USER}'s default account."
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
  local mem="$1"
  if ! is_perlmutter; then  # cori
    if is_group_member "gtrnd"; then
      if needs_full_node "$mem"; then
        echo "genepool"
      else
        echo "genepool_shared"
      fi
      return
    fi
  fi
  if needs_full_node "$mem"; then
    echo "regular"
  else
    echo "shared"
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

check_memory_limit() {
  if [ "$1" -gt "$(total_mem)" ]; then
    >&2 echo "ERROR: maximum memory for this system is $(total_mem) GB."
    die
  fi
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

check_alignment_dir_does_not_exist() {
  if [ -d "$1" ]; then
    >&2 echo "ERROR: Output directory already exists. Not overwriting:"
    >&2 echo "       ${1}"
    >&2 echo "       Consider incrementing rt_alignment_number."
    die
  fi
}

check_experiment_dir_does_exist() {
  if ! find "${raw_data}" -type d -name "${1}" -print -quit | grep -q '.'; then
    >&2 echo "ERROR: directory for experiment '${1}' was not found."
    die
  fi
}

check_experiment_contains_h5_file() {
  local experiment_dir
  experiment_dir="$(find "${raw_data}" -type d -name "${1}" -print -quit)"
  if ! find "${experiment_dir}" -maxdepth 1 -name '*.h5' -print -quit | grep -q '.'; then
    >&2 echo "ERROR: directory for experiment '${experiment_dir}' contains no .h5 files."
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

check_rt_alignment_number_is_non_neg_int() {
  re='^[0-9]+$'
  if ! [[ $1 =~ $re ]] ; then
    >&2 echo ""
    >&2 echo "ERROR: rt_alignment_number must be a non-negative integer"
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
  if getopts m:p:y: option; then
    case $option in
      m) memory="$OPTARG";;
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
rt_alignment_number="${positional_parameters[2]:-0}"
project_dir="${positional_parameters[3]:-$HOME/metabolomics_data}"

IFS='_' read -ra TOKENS <<< "$exp"
proposal="${TOKENS[3]:-}"
exp_token="${TOKENS[4]:-}"
sample_set="${TOKENS[5]:-}"
exp_check_len="${TOKENS[8]:-}"

script_dir="$(get_script_dir)"
short_id="${proposal}_${exp_token}_${sample_set}"
exp_dir="${project_dir}/${exp}"
alignment_dir="${exp_dir}/${workflow_name}/${rt_alignment_number}"
alignment_dir="${exp_dir}/${USER}_${workflow_name}_${rt_alignment_number}_0"

check_exp_id_has_atleast_9_fields "$exp_check_len"
check_experiment_dir_does_exist "$exp"
check_experiment_contains_h5_file "$exp"
check_alignment_dir_does_not_exist "$alignment_dir"
check_rt_alignment_number_is_non_neg_int "$rt_alignment_number"
check_yaml_is_valid "$(echo "${YAML_BASE64:-}" | base64 --decode)"
check_gdrive_authorization
check_not_in_commom_software_filesystem
check_memory_limit "$memory"

account="$(get_slurm_account)"
cpus_requested="$(get_num_cpus "$memory")"
queue="$(get_slurm_queue "$memory")"
constraint="$(get_slurm_constraint)"
time="$(get_slurm_time)"
IFS=$' ' flags="${account:+--account=$account} --qos=${queue} --cpus-per-task=${cpus_requested} --constraint=${constraint} --time=${time}"

IN_FILE="/src/notebooks/reference/RT-Alignment.ipynb"
OUT_FILE="${exp_dir}/${short_id}_${workflow_name}_RT-Alignment_SLURM.ipynb"

PARAMETERS+=" -p experiment ${exp} \
	      -p workflow_name ${workflow_name} \
	      -p project_directory ${project_dir} \
	      -p max_cpus ${threads_to_use}
	      -p rt_alignment_number ${rt_alignment_number}"
if [ ${#extra_parameters[@]} -eq 0 ] || \
   ! contains_parameter_name extra_parameters 'config_file_name'; then
  PARAMETERS+=" -p config_file_name /global/cfs/cdirs/m2650/targeted_analysis/metatlas_config.yaml"
fi
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

mkdir -p "$alignment_dir"
# shellcheck disable=SC2086
sbatch $flags -J "${short_id}_${workflow_name}" "${script_dir}/slurm_template.sh"
