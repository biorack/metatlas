#!/bin/bash

# this script is run by user svc-metabolomics@lbl.gov on metabo-dtn.jgi.lbl.gov

# # crontab entries:
# */5 * * * * $HOME/bin/sync.sh jgi
# */5 * * * * $HOME/bin/sync.sh egsb

netapp_mount="/net/storage.jgi.lbl.gov/vol_metabolomics"

rsyncd_user="msdata"
dest_host="dtn04.nersc.gov"
dest_port="39653"

password_file="/home/svc-metabolomics__lbl.gov/msdata_pw"

vol="${1:-}"

if [ "$#" -ne 1 ] || ( [ "$vol" != "jgi" ] && [ "$vol" != "egsb" ] ); then
    >&2 echo "Usage $0: vol_name"
    >&2 echo "    where vol_name is one of 'jgi' or 'egsb'"
    exit 128
fi

src="${netapp_mount}/${vol}/"
dest="rsync://${rsyncd_user}@${dest_host}:${dest_port}/${vol}/"
log_file="${HOME}/logs/${vol}.log"
log_error_file="${HOME}/logs/${vol}.error.log"
lock_file="/tmp/msdata.${vol}.lockfile"

flock -n "$lock_file" stdbuf -oL \
   rsync \
      --archive \
      --no-perms \
      --no-group \
      --chmod=ugo=rwX \
      --prune-empty-dirs \
      --exclude=\$RECYCLE.BIN/ \
      --exclude=filechecks/ \
      --exclude=Robo_copy_test/ \
      --include="*/" \
      --include="*.raw" \
      --exclude="*" \
      "--password-file=$password_file" \
      "--log-file=${log_file}" \
      "$src" "$dest" 2>&1 | ts >> "$log_error_file"
