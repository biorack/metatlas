#!/usr/bin/env bash
set -euo pipefail

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

shifter --clearenv --module=none --entrypoint \
	--image=ghcr.io/biorack/metatlas/metatlas_shifter:latest \
	python "${script_dir}/validate_config.py" "$@"
