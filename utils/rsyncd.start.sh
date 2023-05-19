#!/bin/bash

# add this to msdata's crontab on dtn04.nersc.gov:

# */5 * * * * /bin/flock -w 0 /tmp/msdata.rsyncd.lock /global/homes/m/msdata/rsyncd.start.sh


set -euo pipefail

umask 027
machine="$(uname -n)"
sg metatlas "rsync --daemon --config /global/common/software/m2650/metatlas-repo/utils/rsyncd.${machine}.conf --port=39653"
