#!/bin/bash
set -euf -o pipefail
set -o xtrace

# dumps production DB on nerscdb04.nersc.gov
# loads the dump into a local, dockerized mysql instance
# transforms (filters) the database using script_file
# exports the database to a sqlite3 file
# uploads the sqlite3 file to a web accessible directory


if [ "$#" -ne 2 ]; then
    echo "Usage $0 script_file mysql_password"
    exit 1
fi

# cat /global/cfs/cdirs/metatlas/mysql_user.txt

ssh_env() {
  if which nersc; then
    # shellcheck disable=SC2068
    nersc $@
  else
    # shellcheck disable=SC2068
    $@
  fi
}

rm -rf dumps
mkdir -p dumps

TIMESTAMP="$(date "+%Y-%m-%d-%H-%M")"
DUMP_FILE="meta_atlas_all_dbs-${TIMESTAMP}.sql"
# shellcheck disable=SC2029 # $DUMP_FILE should expand locally...
ssh_env ssh cori.nersc.gov "shifter --image=docker:mysql/mysql-server:5.7.14 mysqldump \
	-h nerscdb04.nersc.gov -u meta_atlas_admin --all-databases \
	--set-gtid-purged=OFF \"--result-file=$DUMP_FILE\""
ssh_env scp "dtn01.nersc.gov:$DUMP_FILE" "dumps/$DUMP_FILE"
# shellcheck disable=SC2029 # $DUMP_FILE should expand locally...
ssh_env ssh cori.nersc.gov "rm '$DUMP_FILE'"
MYSQL_ID="$(docker run -d --rm -e "MYSQL_ROOT_PASSWORD=$2" -v "$(pwd)/dumps:/docker-entrypoint-initdb.d" -v "$(pwd):/script" mysql:5.7-debian)"
docker exec -it "$MYSQL_ID" /script/internal_extract.sh "$1" "$2"
docker stop "$MYSQL_ID"
SQLITE=$(find ~+/dumps -name '*.sqlite3' | head -1)
ssh_env scp -C "$SQLITE" dtn01.nersc.gov:/global/cfs/cdirs/m2650/www/metatlas/test_data/
