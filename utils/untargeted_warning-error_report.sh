#!/bin/bash

# Debug: Check if the log file exists
if [ ! -f /global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline_log.txt ]; then
    echo "Log file does not exist."
    exit 1
fi

# Check if number of untargeted run cycles is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <cycles>"
    echo "Where <cycles> is the number of days to look back"
    exit 1
fi

cycles=$1

printf 'Error and warning report for untargeted pipeline. Ran on %s\n\n' "$(date)"

printf 'Checked the last %s untargeted pipeline runs.\n\n' "$cycles"

errors=$(awk -v cycles="$cycles" '
/Successfully started/ { positions[++count] = NR }
{ lines[NR] = $0 }
END {
    if (count >= cycles) {
        start = positions[count - (cycles - 1)]
        for (i = start + 1; i <= NR; i++) {
            print lines[i]
        }
    }
}
' /global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline_log.txt | grep -E "WARNING|CRITICAL|ERROR")

if [ -z "$errors" ]; then
    printf "There were no warnings or errors!\n\n"
else
    printf 'Errors or warnings that occurred during last untargeted pipeline cycle:\n\n'
    printf '%s\n\n' "$errors"
    printf 'Check log here: /global/cfs/cdirs/m2650/untargeted_logs/untargeted_pipeline_log.txt\n\n'
fi
