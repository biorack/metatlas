#!/bin/bash
set -euo pipefail

# Check if number of untargeted run cycles is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <cycles>"
    echo "Where <cycles> is the number of days to look back"
    exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
addresses="$(cat "/global/cfs/cdirs/metatlas/raw_data/email_untargeted_errors-warnings")"

# shellcheck disable=SC2086
"${SCRIPT_DIR}/untargeted_warning-error_report.sh" "$1" | \
  mailx -s "untargeted pipeline warnings and errors report" -R "bkieft@lbl.gov" ${addresses}

printf "INFO: report generation and emailing complete\n" | ts
