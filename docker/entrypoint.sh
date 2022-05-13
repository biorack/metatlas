#!/bin/bash

WORK_DIR="$(mktemp -d /tmp/metatlas.XXXXXXXXXX)"

if [[ ! "$WORK_DIR" || ! -d "$WORK_DIR" ]]; then
  echo "Could not create temp dir"
  exit 1
fi

export PYTHONPATH="${WORK_DIR}:${PYTHONPATH}"

function cleanup {
  rm -rf "$WORK_DIR"
}

trap cleanup EXIT

cp -a "$1/." "$WORK_DIR"

shift

export HDF5_USE_FILE_LOCKING="FALSE"
export OMP_NUM_THREADS=1
export OMP_PLACES="threads"
export OMP_PROC_BIND="spread"

# don't do the usual 'exec "$@"' here as that will break the trap
$@
