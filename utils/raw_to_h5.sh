#!/bin/bash
set -Eeuo pipefail
IFS=$'\n\t'

raw_image="quay.io/biocontainers/thermorawfileparser@sha256:3b930ef774b3d4e0d559f38903da2390f9b24b96a016a1761805b88ae78c2b40"

if [ "$#" -ne 1 ]; then
    >&2 echo "Usage $0: raw_ms_file"
    exit 128
fi

raw_file="$(realpath "$1")"
mzml_file="${raw_file%.raw}.mzML"
h5_file="${raw_file%.raw}.h5"

mzml_staged="${mzml_file}.staged"
h5_staged="${h5_file}.staged"

progress_file="${raw_file%.raw}.progress"
failure_file="${raw_file%.raw}.failed"
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

function process_output {
  tee -a "${1:-${progress_file}}"
}

touch "${progress_file}"
chgrp metatlas "${progress_file}" 2>/dev/null || true
chmod 640 "${progress_file}"

function on_error {
  # Replace .progress with .failed if failure occurs
  [ -f "${progress_file}" ] && mv "${progress_file}" "${failure_file}"
  printf "INFO: Reached end of on_error function.\n" | process_output "${failure_file}"
}
trap on_error ERR

function on_term {
  printf "INFO: received TERM signal.\n" | process_output "${progress_file}"
}
trap on_term SIGTERM

function metatlas {
  shifter "--image=ghcr.io/biorack/metatlas/metatlas_shifter:latest" \
          "--env=HDF5_USE_FILE_LOCKING=FALSE" \
          "--env=PYTHONPATH=/src" \
	  "--clearenv" \
	  "--module=none" \
	  "$@" 2>&1 | process_output
}

metatlas "${script_dir}/validate_file_name.py" "$raw_file"

# Convert RAW to mzML
shifter "--clearenv" \
        "--module=none" \
	"--image=${raw_image}" \
	ThermoRawFileParser.sh "-i=${raw_file}" "-o=$(dirname "$raw_file")" -f=1 2>&1 | \
  sed 's%^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \([A-Z]+\)%\1, ThermoRawFileParser:%' | \
  process_output

if [ -f "${mzml_file}" ]; then
    mv "${mzml_file}" "${mzml_staged}"
fi

if [ ! -f "${mzml_staged}" ] || [ ! -s "${mzml_staged}" ]; then
   printf "ERROR: raw to mzml conversion failed for %s\n" "$raw_file" | process_output
   exit 1
fi

# Convert mzML to H5
metatlas "${script_dir}/mzml_to_h5.py" "${mzml_staged}"

# Python script creates .h5 by replacing .mzML with .h5
if [ -f "${h5_staged}" ]; then
    mv "${mzml_staged}" "${mzml_file}"
    mv "${h5_staged}" "${h5_file}"
else
    printf "ERROR: h5 file not generated\n" | process_output
    exit 1
fi

printf "INFO: raw to h5 conversion successfully completed.\n" | process_output "${progress_file}"

# Success: Delete the .progress marker now that final files are written
rm -f "${progress_file}"