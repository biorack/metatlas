#!/usr/bin/env bash
# ----------------------------------------------------------------------
# aws_sync_report.sh
#
#  Email a weekly report of the “aws s3 sync” log for the requested
#  data set (jgi|egsb).  The script reads the log file produced by the
#  sync script, extracts entries that fall within the requested number
#  of days, and sends a short summary + the full list of upload lines.
#
#  Usage:
#      aws_sync_report.sh <jgi|egsb> <days>
#
#  Example:
#      aws_sync_report.sh jgi 7
# ----------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# --------------------------- argument parsing -------------------------
if [[ $# -ne 2 ]] || { [[ $1 != jgi ]] && [[ $1 != egsb ]]; }; then
    >&2 printf "Usage: %s <jgi|egsb> <days>\n" "$(basename "$0")"
    exit 128
fi

DATASET=$1               # jgi or egsb
DAYS_BACK=$2            # how many days back to look

# ---------------------------- constants -------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
LOGFILE="/global/cfs/cdirs/m2650/aws_s3_sync_log/${DATASET}.log"
ADDRESS_FILE="/global/cfs/cdirs/metatlas/raw_data/email_${DATASET}_reports"
FROM_ADDR="bkieft@lbl.gov"

# ----------------------------------------------------------------------
# Helper: turn a log line timestamp (e.g. "Dec 12 12:27:18") into epoch secs.
# The log always starts with:  <Mon> <Day> <HH:MM:SS>
# ----------------------------------------------------------------------
epoch_from_log_line() {
    # $1 = line
    local trimmed_line="${1#"${1%%[![:space:]]*}"}"
    
    # Extract the LAST valid timestamp pattern (handles repeated timestamps)
    local timestamp
    timestamp=$(echo "$trimmed_line" | grep -oE '[A-Za-z]{3} [0-9]{1,2} [0-9]{2}:[0-9]{2}:[0-9]{2}' | tail -n 1)
    
    # If no valid timestamp found, return empty
    if [[ -z $timestamp ]]; then
        return 1
    fi
    
    # Try to convert to epoch, return empty on failure
    local epoch
    if ! epoch=$(date -d "$timestamp" +%s 2>/dev/null); then
        return 1
    fi
    
    echo "$epoch"
}

# ---------------------- compute cutoff epoch -------------------------
NOW_EPOCH=$(date +%s)
CUTOFF_EPOCH=$(( NOW_EPOCH - DAYS_BACK*86400 ))

# --------------------------- scan log --------------------------------
declare -a upload_lines=()
declare -a error_lines=()

while IFS= read -r line; do
    # Skip empty lines
    [[ -z $line ]] && continue

    # Try to extract epoch from line
    if ! line_epoch=$(epoch_from_log_line "$line"); then
        printf "Failed to parse date from line: %s\n" "$line"
        continue
    fi

    if (( line_epoch < CUTOFF_EPOCH )); then
        continue
    fi    

    if [[ $line == *"upload:"* ]]; then
        upload_lines+=("$line")
    elif [[ $line == *"error:"* ]]; then
        error_lines+=("$line")
    fi
done < <(tail -n 20000 "$LOGFILE" | grep -E "\.raw|error")

# ----------------------- build the report ----------------------------
NUM_UPLOADS=${#upload_lines[@]}
NUM_ERRORS=${#error_lines[@]}
NUM_TOTAL=$(( NUM_UPLOADS + NUM_ERRORS ))

summary_line="uploaded=$NUM_UPLOADS errored=$NUM_ERRORS total=$NUM_TOTAL"
report_body=$(
    printf "AWS S3 sync report for %s (last %s days)\n" "$DATASET" "$DAYS_BACK"
    printf "Uploads   : %s\n" "$NUM_UPLOADS"
    printf "Errors    : %s\n" "$NUM_ERRORS"
    printf "Total syncs examined : %s\n\n" "$NUM_TOTAL"
)

# -------------------------- email handling --------------------------
subject="${DATASET} sync report: ${NUM_UPLOADS} uploaded, ${NUM_ERRORS} errors"

# Load recipient list (space‑separated, one per line), skip commented lines
addresses=$(grep -v '^[[:space:]]*#' "$ADDRESS_FILE" | grep -v '^[[:space:]]*$' | tr '\n' ' ')

# Send the e‑mail.  Using `mailx -s "$subject" -R "$FROM_ADDR" $addresses`
printf '%s\n' "$report_body" |
    mailx -s "$subject" -R "$FROM_ADDR" $addresses

# --------------------------- final log -------------------------------
printf "INFO: %s report for %s (last %s days) - %s uploaded, %s errors - emailed to %s\n" \
    "$DATASET" "$(date '+%Y-%m-%d')" "$DAYS_BACK" "$NUM_UPLOADS" "$NUM_ERRORS" "$addresses" | ts