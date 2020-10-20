#!/usr/bin/env python
import argparse
import sys
import sqlite3
from sqlite3 import Error


class sql_connect:


    def __init__(self, db):
        # create a database connection
        self.conn = self.create_connection(db)

        sql_create_metadata_table = """ CREATE TABLE IF NOT EXISTS metadata (
                                        id integer PRIMARY KEY,
                                        limskey integer NOT NULL,
                                        task text NOT NULL,
                                        input text NOT NULL,
                                        output text NOT NULL,
                                        filesize integer,
                                        error text,
                                        passed text
                                    ); """
        if self.conn is not None:
                # create projects table
                self.create_table(self.conn, sql_create_metadata_table)
        else:
          print("Error! cannot create the database connection.")

    def create_connection(self,db_file):
        """ create a database connection to the SQLite database
            specified by db_file
        :param db_file: database file
        :return: Connection object or None
        """
        conn = None
        try:
            conn = sqlite3.connect(db_file)
            return conn
        except Error as e:
            print(e)
    
        return conn

    
    def create_table(self, conn, create_table_sql):
        """ create a table from the create_table_sql statement
        :param conn: Connection object
        :param create_table_sql: a CREATE TABLE statement
        :return:
        """
        try:
            c = conn.cursor()
            c.execute(create_table_sql)
        except Error as e:
            print(e)
    
    def insert_bulk(self, sql_insert, records):
        """ insert data into table using a list of lists (i.e. executemany). This should be an efficient bulk insert.
        :param records: a list of list
        :return:
        """
        try:
            c = self.conn.cursor()
            c.executemany(sql_insert, records)
            self.conn.commit()

        except Error as e:
            print(e)

    def insert_data(self, sql_insert):
        """ insert data into a table 
        :param conn: Connection object
        :param sql_insert: a INSERT INTO statement
        :return:
        """
        try:
            c = self.conn.cursor()
            c.execute(sql_insert)
            #print(sql_insert)
        except Error as e:
            print(e)
    
    def commit(self):
        self.conn.commit()

    def select(self, sql_select):
        """ select data from a table 
        :param conn: Connection object
        :param sql_select: a SELECT FROM statement
        :return: results of the select statement 
        """
        try:
            c = self.conn.cursor()
            c.execute(sql_select)
    
            rows = c.fetchall()
    
            for row in rows:
                print(row)
    
        except Error as e:
            print(e)


    def update(self, sql_update):
        """ update data from a table 
        :param sql_update: a UPDATE table SET x=y WHERE z=m
        :return: 
        """
        try:
            c = self.conn.cursor()
            c.execute(sql_update)
    
        except Error as e:
            print(e)

def main():

    # parse arguments
    sqlite_db=''
    parser = argparse.ArgumentParser(description='This file contains functions to create, update, and read tables from an sqlite3 database.')
    parser.add_argument("-d","--db",help='name of the sqlite database that will be used', type=str, required=True )
    args = parser.parse_args()


    sql_create_metadata_table = """ CREATE TABLE IF NOT EXISTS metadata (
                                        id integer PRIMARY KEY,
                                        limskey integer NOT NULL,
                                        task text NOT NULL,
                                        input text NOT NULL,
                                        output text NOT NULL,
                                        filesize text,
                                        error text,
                                        passed text
                                    ); """


    # create a database connection
    conn = sql_connect(args.db)

    # insert
    #INSERT INTO langs(name) VALUES(?)
    k=123
    i='/some/path/in'
    o='/some/path/out'
    t='raw_to_mzml'
    sql_insert="INSERT INTO metadata(limskey,task,input,output) VALUES(%s,\"%s\",\"%s\",\"%s\")" % (k,t,i,o)
    conn.insert_data(sql_insert)

    records = [(1,2,3,'y'),(4,5,6,'y'),(7,8,9,'y')]
    sql_bulk='INSERT INTO metadata(limskey,task,input,output) VALUES(?,?,?,?)'
    conn.insert_bulk(sql_bulk, records)

    sql_update='UPDATE metadata SET filesize = 1234132,error = "some error" WHERE limskey = 123'
    conn.update(sql_update)

    sql_select="""SELECT * FROM metadata"""
    conn.select(sql_select)

    conn.commit()

if __name__ == '__main__':
    main()


