#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ] || { [ "$1" != "jgi" ] && [ "$1" != "egsb" ]; }; then
    >&2 echo "Usage $0 search_directory [days]"
    >&2 echo "    search_directory is one of 'jgi' or 'egsb'"
    >&2 echo "    days is how many days back to report on (default 7)"
    exit 128
fi

base_dir="/global/cfs/cdirs/metatlas/raw_data/$1"
days="${2:-7}"

function get_filenames {
  # shellcheck disable=SC2048,SC2086
  find "$base_dir" -mindepth 2 -maxdepth 2 -type f \
	          -mtime "-${days}" -name "$@"
}

function parent_dir {
  xargs -d '\n'  --no-run-if-empty -l dirname \
  |  xargs -d '\n' --no-run-if-empty -l basename
}

function format_runs {
  sed 's%.failed$%.raw%' | sed "s%${base_dir}/%    %"
}
num_converted="$(get_filenames '*.h5' -printf '.' | wc -c)"
num_failed="$(get_filenames '*.failed' -printf '.' | wc -c)"
num_in_progress="$(get_filenames '*.progress' -printf '.' | wc -c)"
num_to_convert="$(find "$base_dir" -mindepth 2 -maxdepth 2 -type f \
     \( -name '*.raw' -o -name '*.h5' -o -name '*.failed' -o -name '*.progress' \) | \
  sed -E 's%.(h5|failed|progress)$%.raw%' | \
  sort | \
  uniq -u | \
  wc -l)"

readarray -t failed < <(get_filenames '*failed' | sort)
if [ "${#failed[@]}" -eq 0 ]; then
  failed_without_error=()
else
  readarray -t failed_without_error < <(grep -L -i 'ERROR' "${failed[@]}" | format_runs || true)
fi
num_failed_without_error="${#failed_without_error[@]}"

printf 'File conversion report for %s data. Ran on %s\n\n' "$1" "$(date)"

printf 'Conversions not yet attempted: %s\n' "$num_to_convert"
printf 'Conversions currently in progress: %s\n' "$num_in_progress"
printf '\n'
printf 'In the past %s days...\n' "$days"
printf '    successfull conversions: %s\n' "$num_converted"
printf '    failed conversions: %s\n' "$num_failed"
if [ "$num_failed_without_error" -gt 0 ]; then
  printf '    failed conversions without error messages: %s\n' "$num_failed_without_error"
fi
printf '\n'

[ "$num_converted" = "0" ] && [ "$num_failed" = "0" ] && \
  [ "$num_failed_without_error" = "0" ] && exit 0

converted="$(get_filenames '*.h5' | sed 's%.h5$%.raw%'| sort)"
converted_exp="$(printf '%s\n' "$converted" | parent_dir | uniq -c)"
if [ "$num_converted" -ne "0" ]; then
  printf 'Successful conversions per experiment:\n%s\n\n' "$converted_exp"
fi

[ "$num_failed" = "0" ] && [ "$num_failed_without_error" = "0" ] && exit 0

failed_exp="$(printf '%s\n' "${failed[@]}" | parent_dir | uniq -c)"
formated_failed="$(printf '%s\n' "${failed[@]}" | format_runs)"

# shellcheck disable=SC2086
errors="$(grep 'ERROR' "${failed[@]}" | sed 's%^%    %' | \
	  sed "s%${base_dir}/%%" | \
	  sed -E 's%\.failed.*ERROR(,|)%.raw:\n        %' | \
	  sed -E "s%([-: ]*|)${base_dir}/[^[:space:]].*%%" )" || true

printf 'Failed conversions per experiment:\n%s\n\n' "$failed_exp"
printf 'Failed files:\n%s\n\n' "$formated_failed"

if [ "$num_failed_without_error" -gt 0 ]; then
  printf 'Failed files without error messages:\n'
  printf '    %s\n' "${failed_without_error[@]}"
fi

if [ -n "$errors" ]; then
  printf 'Extracted ERROR lines from .failed files:\n%s\n' "$errors"
fi
