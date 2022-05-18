#!/bin/bash

shifter -e PYTHONPATH=/src --image=doejgi/metatlas_shifter:latest python /src/metatlas/scripts/copy_ms_to_new_names.py
