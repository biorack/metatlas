#!/bin/bash

set -euf -o pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage $0 script_file mysql_password"
    exit 1
fi

# cat /global/cfs/cdirs/metatlas/mysql_user.txt

TIMESTAMP="$(date "+%Y-%m-%d-%H-%M")"
DUMP_FILE="meta_atlas_all_dbs-${TIMESTAMP}.sql"
ssh cori.nersc.gov shifter --image=docker:mysql/mysql-server:5.7.14 mysqldump \
	-h nerscdb04.nersc.gov -u meta_atlas_admin --all-databases \
	--set-gtid-purged=OFF "--result-file=$DUMP_FILE"
mkdir -p dumps
rm -rf dump/*
scp "dtn01.nersc.gov:$DUMP_FILE" "dumps/$DUMP_FILE"
MYSQL_ID="$(docker run --rm -e "MYSQL_ROOT_PASSWORD=$2" -v "$(pwd)/dumps:/docker-entrypoint-initdb.d" -v "$(pwd):/script" mysql:5.7)"
docker exec -it "$MYSQL_ID" /script/internal_extract.sh "$1" "$2"
docker stop "$MYSQL_ID"
SQLITE=$(ls -lt1 "dumps/*.sqlite3" | head -1)
scp -C "dumps/$SQLITE" dtn01.nersc.gov:/global/cfs/cdirs/m2650/www/metatlas/test_data/
