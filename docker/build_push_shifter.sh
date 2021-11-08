#!/bin/bash

BRANCH="$(git branch  --show-current)"

if [ "$BRANCH" != "main" ]; then
  >&2 echo "ERROR: You must be on the main branch when building the shifter image."
  exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
"${SCRIPT_DIR}/most_recent_commits.py" > "${SCRIPT_DIR}/.version_info"
"${SCRIPT_DIR}/build_x.sh" --registry docker.io --project wjhjgi --id shifter --tag latest --user wjhjgi
