#!/bin/bash -l

# CRONTAB LOOKS LIKE THIS:
# */15 * * * * /global/u2/b/bpb/repos/metatlas/scripts/run_rawdata_rsync.sh >/dev/null
# COPY FROM FIRST *
# */15 will run it every 15 minutes
# will only email if stderr has content since stdout goes to null

mkdir -p "$HOME/tmp"
PIDFILE="$HOME/tmp/rawdata_rsync.pid"

if [ -e "${PIDFILE}" ] && (ps -u $(whoami) -opid= |
                           grep -P "^\s*$(cat ${PIDFILE})$" &> /dev/null); then
  echo "Already running."
  exit 99
fi


rsync -zavh /project/projectdirs/metatlas/raw_data/ /global/cscratch1/sd/bpb/raw_data

echo $! > "${PIDFILE}"
chmod 644 "${PIDFILE}"
