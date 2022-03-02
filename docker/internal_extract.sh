#!/bin/bash

set -euf -o pipefail
set -o xtrace

if [ "$#" -ne 2 ]; then
    echo "Usage $0 script_file mysql_password"
    exit 1
fi

# install mysql-to-sqlite3
apt-get update
apt-get install -y python3 python3-pip
pip3 install mysql-to-sqlite3==1.4.10 python-slugify==5.0.2

while ! grep "mysqld: Shutdown complete" /var/log/mysql/error.log; do
    echo 'still waiting for mysql server to finish loading data...'
    sleep 5
done

sleep 30  # shouldn't be needed, but not working without it

# wait for mysql to be ready for connections
while ! mysqladmin ping "--password=$2" --silent; do
    echo 'still waiting for mysql server to be ready...'
    sleep 5
done

# reduce database contents using SQL script
# this is a hack - we rerun the filter step until it works
# this is because we don't have a good check for when the database is
# truely ready. But this does eventually work.
while ! mysql "--password=$2" meta_atlas < "/script/$1"; do
    echo "Database filtering failed. Trying again..."
    sleep 15
done

TIMESTAMP="$(date "+%Y-%m-%d-%H-%M")"

# save reduced database to sqlite3 format
mysql2sqlite -f "/docker-entrypoint-initdb.d/meta_atlas-${TIMESTAMP}.sqlite3" -d meta_atlas -u root --mysql-password "$2"
