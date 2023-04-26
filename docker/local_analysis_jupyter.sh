#!/bin/bash
docker run -p 8888:8888 -v /Users/bpb:/home/micromamba -t analysis:2 jupyter notebook --ip=0.0.0.0 --no-browser --allow-root --port=8888 --notebook-dir=/home/micromamba

