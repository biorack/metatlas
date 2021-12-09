#!/bin/bash -l

# user has to be whitelisted to use cori21 for cron
# pasteur and bpb are on the list
# crontab -e from cori21
# CRONTAB LOOKS LIKE THIS:
# */15 * * * * /global/u2/b/bpb/repos/metatlas/scripts/run_file_converter.sh >/dev/null
# COPY FROM FIRST *
# */15 will run it every 15 minutes
# will only email if stderr has content since stdout goes to null

mkdir -p "$HOME/tmp"
PIDFILE="$HOME/tmp/file_converter.pid"

if [ -e "${PIDFILE}" ] && (ps -u $(whoami) -opid= |
                           grep -P "^\s*$(cat ${PIDFILE})$" &> /dev/null); then
  echo "Already running."
  exit 99
fi
#export PYTHONPATH="/global/homes/b/bpb/repos/metatlas:${PYTHONPATH}"
#export PYTHONPATH="/global/cscratch1/sd/bpb/repos/metatlas_backup_20201020/metatlas:${PYTHONPATH}"
LOGFILE="/global/cfs/cdirs/m2650/file_converter_logs/file_converter.log"
MET_PATH=/global/cfs/cdirs/metatlas/raw_data
#MET_PATH=/global/homes/b/bpb

BIN_PATH="shifter --entrypoint --image=wjhjgi/metatlas_shifter:latest /usr/local/bin/python"
#BIN_PATH=/global/common/software/m2650/python3-metatlas-cori
$BIN_PATH -m metatlas.io.file_converter $MET_PATH &> ${LOGFILE} &

echo $! > "${PIDFILE}"
chmod 644 "${PIDFILE}"

