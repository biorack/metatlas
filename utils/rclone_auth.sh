#!/bin/bash

set -euf -o pipefail

config_file="$HOME/.config/rclone/rclone.conf"

cmd=rclone
if ! which "$cmd" > /dev/null 2>&1; then
  cmd="/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone"
fi

# authorize connections to Google Drive and add a folder mapping
"$cmd" config create --all rclone_test drive \
        config_change_team_drive false \
        config_is_local true \
        config_fs_advanced false \
        client_id "" \
        client_secret "" \
        root_folder_id "15F2RhHs4hXxPNoF_UUNPi5n16I0M7oCu" \
        scope "" \
        service_account_file ""

function add_google_drive_folder {
  # args:
  #   name:       RClone identifier to assign for this folder
  #   folder_id:  unique identifier from Google Drive folder URL

  # add an additional folder mapping on Google Drive without re-doing the authorization process
  # This is done in two steps because we are reading and writing to the same file
  append="$(printf '\n[%s]\n' "$1" && \
            grep -A8 "^\[rclone_test\]" "$config_file" | \
            tail -n +2 | \
            sed "s%^root_folder_id = .*%root_folder_id = ${2}%" )"
  echo "$append" >> "$config_file"
}

add_google_drive_folder 'JGI_Metabolomics_Projects' '0B-ZDcHbPi-aqZzE5V3hOZFc0dms'
add_google_drive_folder 'JGI_Analysis_uploads' '19Ofs5AHB3O8-NYApJUwj4YvH8TbKCGJW'
add_google_drive_folder 'EGSB_Data_outputs' '178Qs5KXuyPayTbz6YaZkplIvwqeVZbKa'
