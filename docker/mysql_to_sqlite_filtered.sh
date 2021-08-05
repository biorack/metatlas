#!/bin/bash

set -euf -o pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage $0 script_file mysql_password"
    exit 1
fi

# cat /global/cfs/cdirs/metatlas/mysql_user.txt

TIMESTAMP="$(date "+%Y-%m-%d-%H-%M")"
mkdir -p dumps
docker run -it --rm -v "$(pwd)/dumps:/dumps" mysql:5.7 \
	/usr/bin/mysqldump -h nerscdb04.nersc.gov -u meta_atlas_admin -p --all-databases "--result-file=/dumps/meta_atlas_all_dbs-${TIMESTAMP}.sql"
MYSQL_ID="$(docker run -it --rm -e "MYSQL_ROOT_PASSWORD=$2"-v "$(pwd)/dumps:/docker-entrypoint-initdb.d" -v "$(pwd):/script" mysql:5.7)"
docker exec -it "$MYSQL_ID" /script/internal_extract.sh "$1" "$2"
docker stop "$MYSQL_ID"
