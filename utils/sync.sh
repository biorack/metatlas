#!/bin/bash

# this script is run by user msdata via a cronjob on dtn04.nersc.gov

# # crontab entries:
# */5 * * * * $HOME/bin/sync.sh jgi
# */5 * * * * $HOME/bin/sync.sh egsb >> $HOME/logs/egsb.log 2>> $HOME/logs/egsb.error.log

netapp_mount="/net/storage.jgi.lbl.gov/vol_metabolomics"

rsyncd_user="msdata"
dest_host="dtn04.nersc.gov"
dest_port="39653"

password_file="/home/svc-metabolomics__lbl.gov/msdata_pw"

if [ "$#" -ne 1 ]; then
    >&2 echo "Usage %0: vol_name"
    >&2 echo "    where vol_name is one of 'jgi' or 'egsb'"
    exit 128
fi

vol="$1"

src="${netapp_mount}/${vol}/"
dest="rsync://${rsyncd_user}@${dest_host}:${dest_port}/${vol}/"
log_file="${HOME}/logs/${vol}.log"
log_error_file="${HOME}/logs/${vol}.error.log"
lock_file="/tmp/msdata.${vol}.lockfile"

/usr/bin/flock -n "$lock_file" /bin/rsync -a -no-p --no-g --chmod=ugo=rwX "--password-file=$password_file" "--log-file=${log_file}" "$src" "$dest" 2>> "$log_error_file"
