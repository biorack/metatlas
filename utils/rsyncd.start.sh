#!/bin/bash

set -euo pipefail

umask 026
machine="$(uname -n)"
sg metatlas "rsync --daemon --config /global/common/software/m2650/metatlas-repo/utils/rsyncd.${machine}.conf --port=39653 &"
