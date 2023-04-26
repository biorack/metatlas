#!/bin/bash
#docker run -p 8888:8888 -v /Users/bpb:/home/micromamba -t analysis:2 jupyter notebook --ip=0.0.0.0 --no-browser --allow-root --port=8888 --notebook-dir=/home/micromamba
docker run \
  -p 8888:8888 \
  -v /Users/bpb:/home/micromamba \
  -v /Users/bpb/Downloads/labkey_user.txt:/global/cfs/cdirs/metatlas/labkey_user.txt \
  -v "/Users/bpb/Downloads/ipython to sheets demo-9140f8697062.json":"/global/cfs/cdirs/metatlas/projects/google_sheets_auth/ipython to sheets demo-9140f8697062.json" \
  ghcr.io/biorack/metatlas/analysis:latest \
  jupyter notebook \
  --ip=0.0.0.0 \
  --no-browser \
  --allow-root \
  --port=8888 \
  --notebook-dir=/home/micromamba
