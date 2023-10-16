#!/bin/bash

shifter -e PYTHONPATH=/src --image=ghcr.io/biorack/metatlas/metatlas_shifter:latest python /src/metatlas/scripts/copy_ms_to_new_names.py $@
