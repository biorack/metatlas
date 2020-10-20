#!/usr/bin/env python
import argparse
import sys
import sqlite3
from sqlite3 import Error
import sqlite_functions

# parse arguments
sqlite_db=''
parser = argparse.ArgumentParser(description='This file contains functions to create, update, and read tables from an sqlite3 database.')
parser.add_argument("-d","--db",help='name of the sqlite database that will be used', type=str, required=True )
args = parser.parse_args()
    
    
# create a database connection
conn = sqlite_functions.sql_connect(args.db)
    
sql_update='UPDATE metadata SET filesize = 1234132,error = "some error" WHERE limskey = 123'
conn.update(sql_update)
conn.commit()

sql="""SELECT * FROM metadata"""
conn.select(sql)
    
