#!/bin/bash

set -euf -o pipefail

CMD=rclone
if ! which "$CMD" > /dev/null 2>&1; then
  CMD="/global/cfs/cdirs/m342/USA/shared-envs/rclone/bin/rclone"
fi

"$CMD" config create --all metabolomics drive \
	config_change_team_drive false \
	config_is_local false \
	config_fs_advanced false \
	client_id "" \
	client_secret "" \
	root_folder_id "0B-ZDcHbPi-aqZzE5V3hOZFc0dms" \
	scope "" \
	service_account_file ""
