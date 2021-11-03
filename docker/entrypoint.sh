#!/bin/bash

cp -R "$1" "$2"

shift
shift

exec "$@"
