#!/bin/bash

set -euf -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_DIR=$(dirname "$SCRIPT_DIR")
OUT_DIR="${SCRIPT_DIR}/out"
IMAGE='ghcr.io/biorack/metatlas/metatlas_ci:1.3.0'
PORT=8888

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -i|--image) IMAGE="$2"; shift ;;
    -o|--out_dir) OUT_DIR="$2"; shift ;;
    -h|--help)
        echo -e "$0 [options]"
        echo ""
        echo "   -h, --help              show this command reference"
        echo "   -i, --image string      name of docker image to run (default ${IMAGE})"
        echo "   -o, --out_dir string    directory to write output files (default ${OUT_DIR})"
        exit 0
        ;;
    *)echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

docker run \
   --rm \
   --name "jupyter" \
   -p "${PORT}:${PORT}" \
   -v "${REPO_DIR}:/src" \
   -v "${OUT_DIR}:/out" \
   "$IMAGE"
