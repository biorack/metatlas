FROM mysql:5.7.34

ARG DUMP_FILE=meta_atlas.dump.sql
ENV MYSQL_ROOT_PASSWORD=mypw

COPY $DUMP_FILE /docker-entrypoint-initdb.d/dump.sql
