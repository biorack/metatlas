#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [ "$#" -ne 2 ] || { [ "$1" != "jgi" ] && [ "$1" != "egsb" ]; }; then
    >&2 printf "Usage: %s search_directory days\n" "$0"
    >&2 printf "    search_directory is one of 'jgi' or 'egsb'\n"
    >&2 printf "    days is how many days back to report on\n"
    exit 128
fi

printf "INFO: generating conversion report for %s covering %s days...\n" "$1" "$2" | ts

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
addresses="$(cat "/global/cfs/cdirs/metatlas/raw_data/email_${1}_reports")"

# shellcheck disable=SC2086
"${SCRIPT_DIR}/conversion_report.sh" "$1" "$2" | \
  mailx -s "${1} file conversion report"  -R "wjholtz@lbl.gov" ${addresses}

printf "INFO: report generation and emailing complete\n" | ts
