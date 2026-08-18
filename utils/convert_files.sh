#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [ "$#" -ne 1 ] || { [ "$1" != "jgi" ] && [ "$1" != "egsb" ]; }; then
    >&2 echo "Usage $0: search_directory (jgi|egsb)"
    exit 128
fi

function finish {
  pids="$(jobs -p)"
  [ -n "$pids" ] && kill -TERM "$pids"
}
trap finish EXIT SIGHUP SIGINT SIGTERM

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

base_dir="/global/cfs/cdirs/metatlas/raw_data/$1"
log="/global/cfs/cdirs/m2650/file_converter_logs/${1}.log"
converter="${script_dir}/raw_to_h5.sh"

raw_image_hash='sha256:3b930ef774b3d4e0d559f38903da2390f9b24b96a016a1761805b88ae78c2b40'
shifterimg pull "quay.io/biocontainers/thermorawfileparser@${raw_image_hash}"

# Get a list of all .raw files
all_raws=$(mktemp)
find "$base_dir" -mindepth 2 -maxdepth 2 -type f -name "*.raw" | sort > "$all_raws"

# Get a list of files that already have a result (.h5 or .failed)
done_raws=$(mktemp)
find "$base_dir" -mindepth 2 -maxdepth 2 -type f \( -name '*.h5' -o -name '*.failed' \) | \
  sed -E 's/\.(h5|failed)$/.raw/' | sort > "$done_raws"

# Subtract "Done" from "All" to get only the files that actually need processing
files_to_process=$(comm -23 "$all_raws" "$done_raws")

# Cleanup temporary lists
rm "$all_raws" "$done_raws"

if [ -z "$files_to_process" ]; then
    echo "No files need conversion. Exiting."
    exit 0
fi

# Filter for actual files (safety check) and run with lowest priority
echo "$files_to_process" | while read -r file; do
    if [ -f "$file" ]; then
        echo "$file"
    fi
done | \
  parallel -j 4 --line-buffer --progress "nice -n 19 ionice -c 3 ${converter} '{}'" 2> >(tee -a "$log") | tee -a "$log"