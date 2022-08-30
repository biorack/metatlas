#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if [ "$#" -ne 2 ] || { [ "$1" != "jgi" ] && [ "$1" != "egsb" ]; }; then
    >&2 echo "Usage $0 search_directory days"
    >&2 echo "    search_directory is one of 'jgi' or 'egsb'"
    >&2 echo "    days is how many days back to report on"
    exit 128
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

addresses="$(cat "/global/cfs/cdirs/metatlas/raw_data/email_${1}_reports")"

# shellcheck disable=SC2086
"${SCRIPT_DIR}/conversion_report.sh" "$1" "$2" | \
  mailx -s "${1} file conversion report"  -R "wjholtz@lbl.gov" ${addresses}
