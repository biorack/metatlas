#!/bin/bash

METATLAS_WORKING_SOURCE_DIR="$(mktemp -d /tmp/metatlas.XXXXXXXXXX)"
PYTHONPATH="${METATLAS_WORKING_SOURCE_DIR}:${PYTHONPATH}"

if [[ ! "$METATLAS_WORKING_SOURCE_DIR" || ! -d "$METATLAS_WORKING_SOURCE_DIR" ]]; then
  echo "Could not create temp dir"
  exit 1
fi

function cleanup {
  rm -rf "$METATLAS_WORKING_SOURCE_DIR"
}

trap cleanup EXIT

export PYTHONPATH
export METATLAS_WORKING_SOURCE_DIR

cp -a "$1/." "$METATLAS_WORKING_SOURCE_DIR"

shift

export HDF5_USE_FILE_LOCKING="FALSE"
export OMP_NUM_THREADS=1
export OMP_PLACES="threads"
export OMP_PROC_BIND="spread"

if [ -n "$SOURCE_CODE_VERSION_ID" ]; then
   git checkout -C "$METATLAS_WORKING_SOURCE_DIR" "$SOURCE_CODE_VERSION_ID"
fi

# don't do the usual 'exec "$@"' here as that will break the trap
# shellcheck disable=SC2068
$@
