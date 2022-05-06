#!/bin/bash -l

# user has to be whitelisted to use cori21 for cron
# pasteur and bpb are on the list
# crontab -e from cori21
# CRONTAB LOOKS LIKE THIS:
# */15 * * * * /global/common/software/m2650/metatlas-repo/metatlas/scripts/run_file_converter.sh >/dev/null
# 48 2 * * * /usr/sbin/logrotate --state /global/cfs/projectdirs/m2650/file_converter_logs/logrotate.status /global/cfs/projectdirs/m2650/file_converter_logs/logrotate.conf
# COPY FROM FIRST *
# */15 will run it every 15 minutes
# will only email if stderr has content since stdout goes to null

set -euo pipefail

export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p "$HOME/tmp"
PIDFILE="$HOME/tmp/file_converter.pid"

if [ -e "${PIDFILE}" ] && (ps -u "$(whoami)" -opid= |
                           grep -P "^\s*$(cat ${PIDFILE})$" &> /dev/null); then
  echo "Already running."
  exit 99
fi

LOGFILE="/global/cfs/cdirs/m2650/file_converter_logs/file_converter.log"
RAW_DATA_PATH=/global/cfs/cdirs/metatlas/raw_data

# pass in PYTHONPATH and don't use "--entrypoint" so that it doesn't copy the git repo to /tmp/metatlas.*
BIN_PATH="shifter -e PYTHONPATH=/src --image=wjhjgi/metatlas_shifter:latest /usr/local/bin/python"

$BIN_PATH -m metatlas.io.file_converter "$RAW_DATA_PATH" &>> "${LOGFILE}" &

echo $! > "${PIDFILE}"
chmod 644 "${PIDFILE}"
