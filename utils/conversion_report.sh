#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [ "$#" -ne 1 ] || { [ "$1" != "jgi" ] && [ "$1" != "egsb" ]; }; then
    >&2 echo "Usage $0: search_directory"
    >&2 echo "    where search_directory is one of 'jgi' or 'egsb'"
    exit 128
fi

base_dir="/global/cfs/cdirs/metatlas/raw_data/$1"
days="7"

converted="$(find "$base_dir" -mindepth 2 -maxdepth 2 -type f \
	          -mtime "-${days}" -name '*.h5' | \
             sed 's%.h5$%.raw%' | \
	     sort)"

failed="$(find "$base_dir" -mindepth 2 -maxdepth 2 -type f -mtime "-${days}" \
	       -name '*.failed' | \
          sed 's%.failed$%.raw%' | \
          sort)"

num_failed="$(printf '%s\n' "$failed" | wc -l)"

num_converted="$(find "$base_dir" -mindepth 2 -maxdepth 2 -type f \
	              -mtime "-${days}" -name '*.h5' | \
	         wc -l)"

converted_exp="$(printf '%s\n' "$converted" | \
	         xargs -l dirname | \
	         xargs -l basename | \
	         uniq -c)"

failed_exp="$(printf '%s\n' "$failed" | \
	      xargs -l dirname | \
	      xargs -l basename | \
	      uniq -c)"

formated_failed="$(printf '%s\n' "$failed" | sed "s%${base_dir}/%    %")"

{ printf 'File conversion report for %s data.\n' "$1"; \
  printf 'In the past %s days...\n' "$days"; \
  printf 'successfull conversions: %s\n' "$num_converted"; \
  printf 'failed conversions: %s\n\n' "$num_failed"; \
  printf 'successfull conversions per experiment:\n%s\n\n' "$converted_exp"; \
  printf 'failed conversions per experiment:\n%s\n\n' "$failed_exp"; \
  printf 'Failed files:\n%s\n' "$formated_failed"; }
