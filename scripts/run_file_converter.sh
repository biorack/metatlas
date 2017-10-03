#!/bin/bash -l

mkdir -p "$HOME/tmp"
PIDFILE="$HOME/tmp/file_converter.pid"

if [ -e "${PIDFILE}" ] && (ps -u $(whoami) -opid= |
                           grep -P "^\s*$(cat ${PIDFILE})$" &> /dev/null); then
  echo "Already running."
  exit 99
fi

LOGFILE="/global/project/projectdirs/metatlas/file_converter.log"
MET_PATH=/global/common/software/m2650/python-cori/bin/python
$MET_PATH -m metatlas.file_converter "$MET_PATH/raw_data/" &> ${LOGFILE} &

echo $! > "${PIDFILE}"
chmod 644 "${PIDFILE}"
