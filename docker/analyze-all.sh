#!/bin/bash
# shellcheck disable=SC2086

if [[ $# == 0 ]]; then
	echo "Usage $0 password"
	exit 1;
fi

DB=meta_atlas
OPTS="--password=$1 -D $DB -Bs"

for table in $(mysql $OPTS -e "show tables"); do
	mysql $OPTS -e "analyze table $table";
done

mysql $OPTS -e "SELECT table_name, table_rows FROM information_schema.TABLES WHERE table_schema='$DB' ORDER BY table_rows DESC"
