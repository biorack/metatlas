#!/bin/bash

set -euf -o pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage $0 script_file mysql_password"
    exit 1
fi

# install mysql-to-sqlite3
apt-get update
apt-get install -y python3 python3-pip
pip3 install mysql-to-sqlite3

# wait for mysql to be ready for connections
while ! mysqladmin ping "--password=$2" --silent; do
    sleep 1
done

# reduce database contents using SQL script
mysql "--password=$2" meta_atlas < "/script/$1"

TIMESTAMP="$(date "+%Y-%m-%d-%H-%M")"

# save reduced database to sqlite3 format
mysql2sqlite -f "/docker-entrypoint-initdb.d/meta_atlas-${TIMESTAMP}.sqlite3" -d meta_atlas -u root --mysql-password "$2"
