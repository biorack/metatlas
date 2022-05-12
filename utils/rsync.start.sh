#!/bin/bash

set -euo pipefail

umask 022
machine="$(uname -n)"
sg msdata "rsync --daemon --config /global/common/software/m2650/metatlas-repo/utils/rsyncd.${machine}.nersc.gov.conf --port=39653 &"
