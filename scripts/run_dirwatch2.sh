#!/bin/bash -l
echo "Starting"
touch dirwatch2.tmp
echo "Created temp file"
MET_PATH=/global/project/projectdirs/metatlas
"$MET_PATH/anaconda/bin/python" -m metatlas.directory_watcher "$MET_PATH/raw_data/"
echo "Done"
