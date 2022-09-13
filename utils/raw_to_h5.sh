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
mzml_file="${raw_file%.raw}.mzML"
h5_file="${raw_file%.raw}.h5"
progress_file="${raw_file%.raw}.progress"
failure_file="${raw_file%.raw}.failed"

if [ -f "$h5_file" ]; then
    >&2 echo "INFO: .h5 file found for ${raw_file} - skipping conversion."
   exit 0
fi

if [ -f "$failure_file" ]; then
    >&2 echo "INFO: .failed file found for ${raw_file} - skipping conversion."
   exit 0
fi

function finish {
  if [ "$done" = "true" ]; then
    rm "${progress_file}"
  else
    mv "${progress_file}" "${failure_file}"
  fi
}
trap finish EXIT

done="false"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

time_stamper="${script_dir}/ts.py"

shifter "--env=PYTHONPATH=/src" "--image=doejgi/metatlas_shifter:latest" \
	"${script_dir}/validate_file_name.py" "$raw_file" 2>&1 | \
	"${time_stamper}" | \
	tee "${progress_file}"

# ThermoRawFileParser.sh should return non-zero exit code on error, but it doesn't
# https://github.com/compomics/ThermoRawFileParser/issues/140
# But the mzML to h5 conversion will fail nicely if the mzML files does not exist
shifter "--image=${raw_image}" ThermoRawFileParser.sh \
	"-i=${raw_file}" "-o=$(dirname "$raw_file")" -f=1 2>&1 | \
	sed 's%^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \([A-Z]+\)%\1, ThermoRawFileParser:%' | \
	"${time_stamper}" | \
	tee -a "${progress_file}"

shifter "--env=PYTHONPATH=/src" "--image=doejgi/metatlas_shifter:latest" \
	"${script_dir}/mzml_to_h5.py" "${mzml_file}" | \
	"${time_stamper}" | \
	tee -a "${progress_file}"

# if we made it here, none of the commands failed and the conversion is complete
done="true"
