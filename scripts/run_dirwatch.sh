#!/bin/bash -l

module load python/2.7.3
if [ "$NERSC_HOST" == "carver" ] || [ "$NERSC_HOST" == "edison" ]; then
    source "/global/project/projectdirs/metatlas/scidb_load_tools/venv_carver/bin/activate"
else
    source "/global/project/projectdirs/metatlas/scidb_load_tools/venv_portal_auth/bin/activate"
fi
python "/global/project/projectdirs/metatlas/scidb_load_tools/directory_watcher.py" "/global/project/projectdirs/metatlas/original_data/"
