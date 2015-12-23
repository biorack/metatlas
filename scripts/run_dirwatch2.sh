#!/bin/bash -l

mkdir -p "$HOME/tmp"
PIDFILE="$HOME/tmp/dirwatch2.pid"

if [ -e "${PIDFILE}" ] && (ps -u $(whoami) -opid= |
                           grep -P "^\s*$(cat ${PIDFILE})$" &> /dev/null); then
  echo "Already running."
  exit 99
fi

LOGFILE="/global/project/projectdirs/metatlas/dirwatch2.log"
MET_PATH=/global/project/projectdirs/metatlas
"$MET_PATH/anaconda/bin/python" -m metatlas.directory_watcher "$MET_PATH/raw_data/lpsilva/20150504_LPSilva_Actino_HILIC_POS_51isolates" &> ${LOGFILE} &

echo $! > "${PIDFILE}"
chmod 644 "${PIDFILE}"
