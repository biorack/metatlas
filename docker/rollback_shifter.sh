#!/bin/bash

#### Rollback will not be complete until "shifterimg pull ..." runs
#### on NERSC. There is a cron job that does this every 5 minutes
#### on cori. Even once the "shifterimg pull ..." runs, users
#### with open notebooks will still need to reload their kernel
#### in order to get back to the older docker image.

set -euf -o pipefail

SHORT="metatlas_shifter"
LONG="wjhjgi/$SHORT"
TAG=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -t|--tag) TAG="$2"; shift ;;
    -h|--help)
        echo -e "$0 [options]"
        echo ""
        echo "   -h, --help          show this command reference"
        echo "   -t, --tag string    tag value to revert to" 
        exit 0
        ;;
    *)echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [ -z "$TAG" ]; then
  >&2 echo "ERROR: You must supply a tag value."
  exit 1
fi

FULL_ID="${LONG}:${TAG}"
LATEST_ID="${LONG}:latest"

docker pull "$FULL_ID"
docker tag "$FULL_ID" "$LATEST_ID"
docker push "$LATEST_ID"
