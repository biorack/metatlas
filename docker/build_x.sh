#!/bin/bash

# if image registry requires authorization, login before running this script

set -euf -o pipefail

PROJECT="doejgi"
REGISTRY="docker.io"
DOCKER="docker"
TAG=""
ID=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -d|--docker) DOCKER="$2"; shift ;;
    -i|--id) ID="$2"; shift ;;
    -p|--project) PROJECT="$2"; shift ;;
    -r|--registry) REGISTRY="$2"; shift ;;
    -t|--tag) TAG="$2"; shift ;;
    -h|--help)
        echo -e "$0 [options]"
        echo ""
        echo "   -h, --help              show this command reference"
	echo "   -d, --docker            name of docker command (default ${DOCKER})"
	echo "   -i, --id                image id"
        echo "   -p, --project string    project name within the registry (default ${PROJECT})"
        echo "   -r, --registry string   FQDN of container registry to push to"
        echo "                           use 'NONE' to not push (default ${REGISTRY})"
        echo "   -t, --tag string        image tag"
        exit 0
        ;;
    *)echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [[ "$TAG" == "" ]]; then
  >&2 echo "ERROR: no tag value given"
  exit 1
fi

if [[ "$ID" == "" ]]; then
  >&2 echo "ERROR: no id value given"
  exit 2
fi

DOCKERFILE="Dockerfile.$ID"
IMAGE_NAME="metatlas_$ID"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
"${SCRIPT_DIR}/build.sh" --image "$IMAGE_NAME" --tag "$TAG" --docker "$DOCKER" --project "$PROJECT" --registry "$REGISTRY" --user "$USER" --dockerfile "${SCRIPT_DIR}/$DOCKERFILE"
