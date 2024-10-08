#!/bin/bash

# Check if a timeback in days is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <timeback> <department>"
    echo "Where <timeback> is the number of days to look back and <department> is the department to filter on ('jgi' or 'eb')"
    exit 1
fi

timeback=$1
department=$2

printf '%s completion report for untargeted pipeline.\n\t- Run on %s\n\n' "${department^^}" "$(date +"%Y-%m-%d")"

# Get the status of all projects
cmd="/global/homes/m/msdata/.conda/envs/python-untargeted-v1/bin/python3.11 /global/common/software/m2650/metatlas-repo/metatlas/untargeted/check_untargeted_status.py --print_recent 100"
temp_status_table=$(mktemp)
$cmd > "$temp_status_table"

# Only get lines with project statuses
temp_status_table_subset=$(mktemp)
grep -E '^[0-9]+' "$temp_status_table" > "$temp_status_table_subset"

# Sort the combined file by the date in the first column in descending order
temp_sorted_file=$(mktemp)
sort -k1,1 "$temp_status_table_subset" > "$temp_sorted_file"

# Extract the first column (timestamps) from the sorted file
temp_timestamps=$(mktemp)
cut -f1 "$temp_sorted_file" > "$temp_timestamps"

# Extract the timestamp from the last line
current_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
current_timestamp=$(date -Iseconds -d "$current_timestamp")

# Define target timestamp
target_timestamp=$(date -d "$current_timestamp - $timeback days" +"%Y-%m-%d %H:%M:%S")

# Initialize variables to track the closest timestamp
closest_timestamp=""
closest_diff=9999999999

# Iterate over each timestamp in the temp file
while IFS= read -r timestamp; do
    # Calculate the difference in seconds between the current timestamp and the target date
    current_diff=$(($(date -d "$timestamp" +%s) - $(date -d "$target_timestamp" +%s)))
    # Take the absolute value of the difference
    current_diff=${current_diff#-}
    # Update the closest timestamp if the current difference is smaller
    if [ "$current_diff" -lt "$closest_diff" ]; then
        closest_diff=$current_diff
        closest_timestamp=$timestamp
    fi
done < "$temp_timestamps"

# Extract the projects that were completed since the closest timestamp
temp_projects_subset=$(mktemp)
awk -v ts="$closest_timestamp" '
    $0 ~ ts { found = 1; next }
    found { print }
' "$temp_sorted_file" > "$temp_projects_subset"

# Get the relevant department for emails
temp_projects_subset_filtered=$(mktemp)
cat "$temp_projects_subset" | awk -F'_' -v dept="$department" 'tolower($2) == tolower(dept)' > "$temp_projects_subset_filtered"

# Count the number of errors and successes in project status table
errors=$(cat "$temp_projects_subset_filtered" | grep "09 error" | wc -l)
errored_projects=$(cat "$temp_projects_subset_filtered" | grep "09 error" | cut -f2 | sed 's/^/\t- /1')
runnings=$(cat "$temp_projects_subset_filtered" | grep -E "04 running|01 initiation" | wc -l)
running_projects=$(cat "$temp_projects_subset_filtered" | grep -E "04 running|01 initiation" | cut -f2 | sed 's/^/\t- /1')
successes=$(cat "$temp_projects_subset_filtered" | grep -P "07 complete\t07 complete\t07 complete\t07 complete\t08 uploaded" | wc -l)
successful_projects=$(cat "$temp_projects_subset_filtered" | grep -P "07 complete\t07 complete\t07 complete\t07 complete\t08 uploaded" "$temp_projects_subset_filtered" | cut -f2,8 | rev | sed 's/ /\t/1' | rev | awk -F'\t' '{print $1" (Uploaded on "$2")"}' | sed 's/^/\t- /1')

# Print the email summary report
printf 'Status of %s untargeted pipline runs in the last %s days\n' "${department^^}" "$timeback"
printf '\nProjects unsuccessfully ended with error status: %s\n' "$errors"
if [ "$errors" -gt 0 ]; then
    printf '%s\n' "$errored_projects"
fi
printf '\nProjects currently running: %s\n' "$runnings"
if [ "$runnings" -gt 0 ]; then
    printf '%s\n' "$running_projects"
fi
printf '\nProjects successfully completed and uploaded: %s\n' "$successes"
if [ "$successes" -gt 0 ]; then
    printf '%s\n' "$successful_projects"
fi

# Clean up temporary files
rm "$temp_status_table" "$temp_status_table_subset" "$temp_sorted_file" "$temp_timestamps" "$temp_projects_subset" "$temp_projects_subset_filtered"
