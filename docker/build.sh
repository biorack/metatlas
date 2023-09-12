#!/bin/bash

# if image registry requires authorization, login before running this script

# catch some common errors, terminate if a command returns non-zero exit code
set -euf -o pipefail

PROJECT="doejgi"
REGISTRY="docker.io"
DOCKER="docker"
DOCKERFILE=""

IMAGE_NAME=""
TAG=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -d|--docker) DOCKER="$2"; shift ;;
    -f|--dockerfile) DOCKERFILE="$2"; shift ;;
    -i|--image) IMAGE_NAME="$2"; shift ;;
    -r|--registry) REGISTRY="$2"; shift ;;
    -p|--project) PROJECT="$2"; shift ;;
    -t|--tag) TAG="$2"; shift ;;
    -h|--help)
        echo -e "$0 [options]"
        echo ""
        echo "   -h, --help              show this command reference"
	echo "   -d, --docker            name of docker command (default ${DOCKER})"
	echo "   -f, --dockerfile        name of Dockerfile"
        echo "   -i, --image string      name of image to build"
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

if [[ "$IMAGE_NAME" == "" ]]; then
  >&2 echo "ERROR: no Dockerfile value given"
  exit 10
fi

if [[ "$DOCKERFILE" == "" ]]; then
  >&2 echo "ERROR: no Dockerfile value given"
  exit 10
fi

if [[ "$TAG" == "" ]]; then
  >&2 echo "ERROR: no tag value given"
  exit 9
fi

SHORT_TAG="${IMAGE_NAME}:${TAG}"
LONG_TAG="${REGISTRY}/${PROJECT}/${SHORT_TAG}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [[ ! -r "${DOCKERFILE}" ]]; then
  >&2 echo "ERROR: Could not find readable Dockerfile at ${DOCKERFILE}."
  exit 1
fi

${DOCKER} image build --tag "${SHORT_TAG}" --file "${DOCKERFILE}" "$SCRIPT_DIR"

if [[ "$REGISTRY" != "NONE" ]]; then
  if [[ $(uname -s) == "Darwin" ]]; then
    # no readlink on macOS...
    if [[ $(basename "$(which "${DOCKER}")") == 'podman' ]]; then
      PUSH_FLAGS="--format=docker"
    fi
  else
    if [[ $(basename "$(readlink -f "$(which "${DOCKER}")")") == 'podman' ]]; then
      PUSH_FLAGS="--format=docker"
    fi
  fi
  ${DOCKER} image tag "${SHORT_TAG}" "${LONG_TAG}"
  # shellcheck disable=SC2086
  ${DOCKER} image push ${PUSH_FLAGS:-} "${LONG_TAG}"
  TAG="${LONG_TAG}"
else
  TAG="${SHORT_TAG}"
fi
