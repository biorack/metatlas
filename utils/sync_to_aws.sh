#!/bin/bash

usage() {
  >&2 echo "ERROR: Usage: $0 jgi|egsb"
  exit 128
}

[ "$#" -ne "1" ] && usage
! { [ "$1" = "jgi" ] || [ "$1" = "egsb" ]; } && usage

nersc_raw_data='/global/cfs/cdirs/metatlas/raw_data'
s3_raw_data='s3://northen-archive/raw_data'

log="/global/cfs/cdirs/m2650/aws_s3_sync_log/${1}.log"

aws s3 sync --storage-class GLACIER \
            --color off \
            --no-progress \
            --exclude '*' \
	    --include '*.raw' \
	    --include '*.cmbx' \
            --include '*/other/*.sld' \
            --include "*/other/*.csv' \
	    --exclude '*/other/*' \
            "${nersc_raw_data}/${1}" \
	    "${s3_raw_data}/${1}" \
	    2>&1 \
| ts \
| tee -a "${log}"
