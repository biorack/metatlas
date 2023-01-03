#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

raw_image_base_name='quay.io/biocontainers/thermorawfileparser'
# this is version 1.4.1 of ThermoRawFileParser:
hash='sha256:2aa0ea0aa2c5f43e5244ed7eaad4a14e6ed6319b91a9d917850ac2b5ff73741c'
raw_image="${raw_image_base_name}@${hash}"

if [ "$#" -ne 1 ]; then
    >&2 echo "Usage $0: raw_ms_file"
    exit 128
fi

raw_file="$(realpath "$1")"
mzml_file="${raw_file%.raw}.mzML"
h5_file="${raw_file%.raw}.h5"
progress_file="${raw_file%.raw}.progress"
failure_file="${raw_file%.raw}.failed"

function on_error {
  rm -rf "${mzml_file}" "${h5_file}"
  mv "${progress_file}" "${failure_file}"
  printf "INFO: Reached end of on_error function.\n" | \
    process_output "${failure_file}"
}
trap on_error ERR

function on_term {
  printf "INFO: recieved TERM signal, cleaning up.\n" | \
       process_output "${progress_file}"
  rm -rf "${mzml_file}" "${h5_file}" "${progress_file}"
}
trap on_term SIGTERM

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

function process_output {
  "${script_dir}/ts.py" | tee -a "${1:-${progress_file}}"
}

rm -rf "${progress_file}"
touch "${progress_file}"
chgrp metatlas "${progress_file}"
chmod 640 "${progress_file}"

if [ -f "$h5_file" ]; then
   printf "INFO: .h5 file found for %s - skipping conversion.\n" "$raw_file" \
      process_output
   exit 0
fi

if [ -f "$failure_file" ]; then
   printf "INFO: .failed file found for %s - skipping conversion.\n" "$raw_file" \
     process_output
   exit 0
fi

function metatlas {
  shifter "--image=doejgi/metatlas_shifter:latest" \
          "--env=HDF5_USE_FILE_LOCKING=FALSE" \
          "--env=PYTHONPATH=/src" \
	  "--clearenv" \
	  "--module=none" \
	  "$@" 2>&1 | \
  process_output
}

shifterimg pull "${raw_image}"

metatlas "${script_dir}/validate_file_name.py" "$raw_file"

shifter "--clearenv" \
        "--module=none" \
	"--image=${raw_image}" \
	ThermoRawFileParser.sh "-i=${raw_file}" "-o=$(dirname "$raw_file")" -f=1 2>&1 | \
  sed 's%^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \([A-Z]+\)%\1, ThermoRawFileParser:%' | \
  process_output

if [ ! -f "${mzml_file}" ] || [ ! -s "${mzml_file}" ]; then
   printf "ERROR: raw to mzml conversion failed for %s, but return code was still 0.\n" "$raw_file" \
     process_output
   exit 1
fi

metatlas "${script_dir}/mzml_to_h5.py" "${mzml_file}"

printf "INFO: raw to h5 conversion sucessfully completed.\n" | \
       process_output "${progress_file}"
rm -rf "${progress_file}"
