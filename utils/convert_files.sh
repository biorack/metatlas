#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [ "$#" -ne 1 ] || { [ "$1" != "jgi" ] && [ "$1" != "egsb" ]; }; then
    >&2 echo "Usage $0: search_directory"
    >&2 echo "    where search_directory is one of 'jgi' or 'egsb'"
    exit 128
fi

function finish {
  module unload parallel
  # kill child processes as well
  pids="$(jobs -p)"
  if [ "$pids" != "" ]; then
    kill -TERM "$pids"
  fi
}
trap finish EXIT SIGHUP SIGINT SIGTERM

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

module load parallel

base_dir="/global/cfs/cdirs/metatlas/raw_data/$1"

log="/global/cfs/cdirs/m2650/file_converter_logs/${1}.log"

converter="${script_dir}/raw_to_h5.sh"

find "$base_dir" -mindepth 2 -maxdepth 2 -type f \
     \( -name '*.raw' -o -name '*.h5' -o -name '*.failed' \) | \
  sed 's%.h5$%.raw%' | \
  sed 's%.failed$%.raw%' | \
  sort | \
  uniq -u | \
  sed "s%^\(.*\)$%${converter} '\1'%" | \
  parallel  >> "$log" 2>&1
