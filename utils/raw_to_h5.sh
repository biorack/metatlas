#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

raw_image_base_name='quay.io/biocontainers/thermorawfileparser'
# this is version 1.4.0 of ThermoRawFileParser:
hash='sha256:76c17e3124b723f271bc3d5cf0555650288676bb1a827bd1bae9bb684444a404'
raw_image="${raw_image_base_name}@${hash}"

if [ "$#" -ne 1 ]; then
    >&2 echo "Usage $0: raw_ms_file"
    exit 128
fi

raw_file="$(realpath "$1")"
failure_file="${raw_file%.raw}.failed"
mzml_file="${raw_file%.raw}.mzML"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

shifter "--env=PYTHONPATH=/src" "--image=doejgi/metatlas_shifter:latest" \
	"${SCRIPT_DIR}/validate_file_name.py" "$raw_file" 2>&1 | \
	"${SCRIPT_DIR}/ts.py" | \
	tee "${failure_file}"

# ThermoRawFileParser.sh should return non-zero exit code on error, but it doesn't
# https://github.com/compomics/ThermoRawFileParser/issues/140
# But the mzML to h5 conversion will fail nicely if the mzML files does not exist
shifter "--image=${raw_image}" ThermoRawFileParser.sh \
	"-i=${raw_file}" "-o=$(dirname "$raw_file")" -f=1 2>&1 | \
	sed 's%^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} INFO%INFO ThermoRawFileParser:%' | \
	"${SCRIPT_DIR}/ts.py" | \
	tee -a "${failure_file}"

shifter "--env=PYTHONPATH=/src" "--image=doejgi/metatlas_shifter:latest" \
	"${SCRIPT_DIR}/mzml_to_h5.py" "${mzml_file}" | \
	"${SCRIPT_DIR}/ts.py" | \
	tee -a "${failure_file}"

# if we made it here, none of the commands failed and the conversion is complete
rm "${failure_file}"
