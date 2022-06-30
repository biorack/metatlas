#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

raw_image='quay.io/biocontainers/thermorawfileparser:1.3.4--ha8f3691_0'

if [ "$#" -ne 1 ]; then
    >&2 echo "Usage $0: raw_ms_file"
    exit 128
fi

raw_file="$(realpath "$1")"

validation="\
from metatlas.tools.validate_filenames import validate_file_name
assert validate_file_name('${raw_file}')"

shifter "--env=PYTHONPATH=/src" "--image=doejgi/metatlas_shifter:latest" \
        python -c "$validation"

# ThermoRawFileParser.sh should return non-zero exit code on error, but it doesn't
# https://github.com/compomics/ThermoRawFileParser/issues/140
# But the mzML to h5 conversion will fail nicely if the mzML files does not exist
shifter "--image=${raw_image}" ThermoRawFileParser.sh \
	"-i=${raw_file}" "-o=$(dirname "$raw_file")" -f=1 2>&1 | \
	ts "%Y-%m-%d %H:%M:%.S"

mzml_file="${raw_file%.raw}.mzML"

mzml_to_h5="\
from metatlas.io.file_converter import mzml_to_h5_and_add_to_db
mzml_to_h5_and_add_to_db('${mzml_file}')"

shifter "--env=PYTHONPATH=/src" "--image=doejgi/metatlas_shifter:latest" \
        python -c "$mzml_to_h5"
