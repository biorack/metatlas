from __future__ import absolute_import
from __future__ import print_function

import inspect
import logging
import sys
import os
import getpass
import six
import uuid
from collections import defaultdict
import functools

import dataset
import pandas as pd

import socket
import os.path
import yaml
from six.moves import input
from sqlalchemy import create_engine

try:
    from traitlets import (
        HasTraits, CUnicode, List, CInt, Instance, Enum,
        CFloat, CBool)
except ImportError:
    from IPython.utils.traitlets import (
        HasTraits, CUnicode, List, CInt, Instance, Enum,
        CFloat, CBool)

logger = logging.getLogger(__name__)

# Whether we are running from NERSC
ON_NERSC = 'METATLAS_LOCAL' not in os.environ
logger.info('NERSC=%s', ON_NERSC)
os.environ["METATLAS_LOCAL"] = "1"

# Observable List from
# http://stackoverflow.com/a/13259435

def callback_method(func):
    def notify(self, *args, **kwargs):
        if not hasattr(self, '_callbacks'):
            return func(self, *args, **kwargs)
        for _, callback in self._callbacks:
            callback()
        return func(self, *args, **kwargs)
    return notify


class NotifyList(list):
    extend = callback_method(list.extend)
    append = callback_method(list.append)
    remove = callback_method(list.remove)
    pop = callback_method(list.pop)
    __delitem__ = callback_method(list.__delitem__)
    __setitem__ = callback_method(list.__setitem__)
    __iadd__ = callback_method(list.__iadd__)
    __imul__ = callback_method(list.__imul__)

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.__class__(list.__getitem__(self, item))
        else:
            return list.__getitem__(self, item)

    def __init__(self, *args):
        list.__init__(self, *args)
        self._callbacks = []
        self._callback_cntr = 0

    def register_callback(self, cb):
        self._callbacks.append((self._callback_cntr, cb))
        self._callback_cntr += 1
        return self._callback_cntr - 1

    def unregister_callback(self, cbid):
        for idx, (i, cb) in enumerate(self._callbacks):
            if i == cbid:
                self._callbacks.pop(idx)
                return cb
        else:
            return None

def set_docstring(cls):
    """Set the docstring for a MetatlasObject object"""
    doc = cls.__doc__
    if not doc:
        doc = cls.__name__ + ' object.'
    doc += '\n\nParameters\n----------\n'
    for (tname, trait) in sorted(cls.class_traits().items()):
        if tname.startswith('_'):
            continue
        descr = trait.__class__.__name__.lower()
        if descr.startswith('c'):
            descr = descr[1:]
        elif descr == 'enum':
            descr = '{' + ', '.join(trait.values) + '}'
        doc += '%s: %s\n' % (tname, descr)
        help_text = trait.help#get_metadata('help')
        if not help_text:
            help_text = '%s value.' % tname
        help_text = help_text.strip()
        if help_text.endswith('.'):
            help_text = help_text[:-1]
        if trait.metadata.get('readonly', False):
            help_text += ' (read only)'
        help_text += '.'
        doc += '    %s\n' % help_text
    cls.__doc__ = doc
    return cls

class Workspace(object):
    instance = None

    def __init__(self):
            self.path = self.get_database_path(with_password=True)
            mysql_kwargs = {
                "pool_recycle": 3600,
                "connect_args": {
                    "connect_timeout": 120,
                    "init_command": "SET SESSION innodb_lock_wait_timeout=10000"  # This sets the production mysql timeout
                }
            }
            sqlite_kwargs = {"connect_args": {"timeout": 7200}}  # This sets the development sqlite timeout
            self.engine_kwargs = sqlite_kwargs if self.path.startswith("sqlite") else mysql_kwargs
            self.tablename_lut = {}
            self.subclass_lut = {}

            # Delayed import to avoid circular import
            from .metatlas_objects import MetatlasObject, MetList, MetInstance, MetFloat, Stub
            for klass in _get_subclasses(MetatlasObject):
                name = klass.__name__.lower()
                self.subclass_lut[name] = klass
                if name.endswith('s'):
                    self.subclass_lut[name + 'es'] = klass
                    self.tablename_lut[klass] = name + 'es'
                else:
                    self.subclass_lut[name + 's'] = klass
                    self.tablename_lut[klass] = name + 's'
            # handle circular references
            self.seen = {}
            Workspace.instance = self

    @classmethod
    def get_instance(cls):
        """Returns a existing instance of Workspace or creates a new one"""
        if Workspace.instance is None:
            return Workspace()
        return Workspace.instance

    def get_database_path(self, with_password: bool = True) -> str:
        """Returns database access string, use with_password=False to hide password"""
        config_dir = os.path.dirname(sys.modules[self.__class__.__module__].__file__)
        if ON_NERSC:
            with open(os.path.join(config_dir, 'nersc_config', 'nersc.yml'), encoding="utf-8") as fid:
                nersc_info = yaml.safe_load(fid)
            if with_password:
                with open(nersc_info['db_passwd_file'], encoding="utf-8") as fid:
                    password = fid.read().strip()
            else:
                password = '***********'
            return f"mysql+pymysql://meta_atlas_admin:{password}@nerscdb04.nersc.gov/{nersc_info['db_name']}"
        local_config_file = os.path.join(config_dir, 'local_config', 'local.yml')
        if os.path.isfile(local_config_file):
            with open(local_config_file, encoding="utf-8") as fid:
                local_info = yaml.safe_load(fid)
            hostname = 'localhost' if 'db_hostname' not in local_info else local_info['db_hostname']
            login = ''
            if 'db_username' in local_info:
                if 'db_password' in local_info:
                    password = local_info['db_password'] if with_password else '***********'
                    login = f"{local_info['db_username']}:{password}@"
                else:
                    login = f"{local_info['db_username']}@"
            return f"mysql+pymysql://{login}{hostname}/{local_info['db_name']}"
        filename = f"{getpass.getuser()}_workspace.db"
        if os.path.exists(filename):
            os.chmod(filename, 0o775)
        return f"sqlite:///{filename}"

    def get_connection(self):
        """
        Get a re-useable connection to the database.
        Each activity that queries the database needs to have this function preceeding it.
        """
        logger.debug("Getting connection to database at: %s with engine kwargs: %s", self.path, self.engine_kwargs)
        return dataset.connect(self.path, engine_kwargs=self.engine_kwargs)

    def convert_to_double(self, table, entry):
        """Convert a table column to double type."""
        db = self.get_connection()
        db.begin()
        try:
            db.query('alter table `%s` modify `%s` double' % (table, entry))
            db.commit()
        except Exception as err:
            rollback_and_log(db, err)
        finally:
            close_db_connection(db)

    def save_objects(self, objects, _override=False):
        """Save objects to the database"""
        logger.debug('Entering Workspace.save_objects')
        if not isinstance(objects, (list, set)):
            objects = [objects]
        self._seen = dict()
        self._link_updates = defaultdict(list)
        self._updates = defaultdict(list)
        self._inserts = defaultdict(list)
        for obj in objects:
            logger.debug("Getting saved data for %s", obj.unique_id)
            self._get_save_data(obj, _override)

        logger.debug('Connecting to database...')
        db = self.get_connection()
        logger.debug("DB: %s", db)
        if db is None:
            logger.error("Failed to establish a database connection.")
            print("Failed to establish a database connection.")
            return
        db.begin()
        logger.debug("Database type: %s", type(db))
        logger.debug("Database path: %s", self.path)

        if self.path.startswith("sqlite"):
            pass
        else:
            # Confirm that session initiation set the correct lock wait time
            lock_wait_time = db.query("SHOW SESSION VARIABLES LIKE '%INNODB_LOCK_WAIT_TIMEOUT%';")
            for row in lock_wait_time:
                logger.debug("Database lock wait time from engine kwargs: %s", row)

            # Set the isolation level to READ-COMMITTED based on
            # https://stackoverflow.com/questions/5836623/getting-lock-wait-timeout-exceeded-try-restarting-transaction-even-though-im
            db.query("SET SESSION transaction_isolation = 'READ-COMMITTED';")

        try:
            for (table_name, updates) in self._link_updates.items():
                logger.debug("Linking updates for object %s in table %s", updates, table_name)
                if table_name not in db:
                    logger.debug("Table %s not in database", table_name)
                    continue
                sql = f"UPDATE `{table_name}` SET source_id = %s WHERE source_id = %s"
                self._execute_bulk(db, sql, [(prev_uid, uid) for (uid, prev_uid) in updates])
                logger.debug("Update query completed")

            for (table_name, updates) in self._updates.items():
                logger.debug("Updating object %s in table %s", updates, table_name)
                if '_' not in table_name and table_name not in db:
                    try:
                        table = db.get_table(table_name, primary_id='unique_id', primary_type=db.types.string(32))
                        # Ensure the table has the required columns
                        for col in updates[0].keys():
                            if col != 'unique_id':
                                table.create_column(col, db.types.text)  # Adjust data types as necessary
                    except Exception as e:
                        logger.error("Failed to create table %s: %s", table_name, e)
                    if 'sqlite' not in self.path:
                        self.fix_table(table_name)
                sql = f"UPDATE `{table_name}` SET unique_id = ? WHERE unique_id = ?"
                self._execute_bulk(db, sql, [(prev_uid, uid) for (uid, prev_uid) in updates])
                logger.debug("Update query completed")

            for (table_name, inserts) in self._inserts.items():
                logger.debug("Table name: %s", table_name)
                logger.debug("Inserts: %s", inserts)
                if table_name not in db:
                    logger.debug("Table %s not in database, creating it.", table_name)
                    # Ensure the table exists with the appropriate columns
                    if inserts:
                        table = db.get_table(table_name, primary_id='unique_id', primary_type=db.types.string(32))
                        logger.debug("Created table %s", table)
                        # Ensure the table has the required columns
                        for col in inserts[0].keys():
                            if col != 'unique_id':
                                table.create_column(col, db.types.text)  # Adjust data types as necessary
                                logger.debug("Added column %s to table %s", col, table_name)
                if inserts:
                    logger.debug("Inserting data into table %s", table_name)
                    columns = ', '.join(inserts[0].keys())
                    placeholders = ', '.join(['?'] * len(inserts[0]))
                    sql = f"INSERT INTO `{table_name}` ({columns}) VALUES ({placeholders})"
                    logger.debug("Running SQL insert query: %s", sql)
                    self._execute_bulk(db, sql, [tuple(row.values()) for row in inserts])
                    logger.debug("Insert query completed")

            logger.debug("Committing changes to database")
            db.commit()
            logger.debug('Exiting Workspace.save_objects')
        except Exception as err:
            rollback_and_log(db, err)
        finally:
            logger.debug('Closing database connection')
            close_db_connection(db)

    def _execute_bulk(self, db, sql, params):
        """Execute a bulk operation using the underlying database connection."""
        logger.debug("Executing bulk operation to db %s with SQL: %s", self.path, sql)
        engine = create_engine(self.path, **self.engine_kwargs)
        with engine.connect() as connection:
            with connection.begin() as transaction:
                try:
                    # Determine the placeholder based on the database type
                    if self.path.startswith("sqlite"):
                        sql = sql.replace('%s', '?')
                    connection.execute(sql, params)
                    transaction.commit()
                except Exception as err:
                    transaction.rollback()
                    raise err
            
    def create_link_tables(self, klass):
        """
        Create a link table in the database of the given trait klass
        """
        name = self.tablename_lut[klass]
        db = self.get_connection()
        db.begin()
        try:
            for (tname, trait) in klass.class_traits().items():
                if isinstance(trait, MetList):
                    table_name = '_'.join([name, tname])
                    if table_name not in db:
                        db.create_table(table_name)
                        link = dict(source_id=uuid.uuid4().hex,
                                    head_id=uuid.uuid4().hex,
                                    target_id=uuid.uuid4().hex,
                                    target_table=uuid.uuid4().hex)
                        db[table_name].insert(link)
            db.commit()
        except Exception as err:
            rollback_and_log(db, err)
        finally:
            close_db_connection(db)

    def _get_save_data(self, obj, override=False):
        """Get the data that will be used to save an object to the database"""
        if obj.unique_id in self._seen:
            return
        if isinstance(obj, Stub):
            return
        name = self.tablename_lut[obj.__class__]
        self._seen[obj.unique_id] = True
        changed, prev_uid = obj._update(override)
        state = dict()
        for (tname, trait) in obj.traits().items():
            if tname.startswith('_'):
                continue
            if isinstance(trait, MetList):
                logger.debug("Getting saved data for List trait")
                # handle a list of objects by using a Link table
                # create the link table if necessary
                table_name = '_'.join([name, tname])
                logger.debug("Table name: %s", table_name)
                if changed and prev_uid:
                    logger.debug("Adding link updates")
                    self._link_updates[table_name].append((obj.unique_id,
                                                           obj.prev_uid))
                else:
                    logger.debug("Skipping link updates")
                value = getattr(obj, tname)
                #logger.debug("List type value to store: %s", value)
                # do not store this entry in our own table
                if not value:
                    continue
                # create an entry in the table for each item
                # store the item in its own table
                for subvalue in value:
                    logger.debug("Adding saved List data for subvalue to Workspace.self")
                    self._get_save_data(subvalue, override)
                    link = dict(source_id=obj.unique_id,
                                head_id=obj.head_id,
                                target_id=subvalue.unique_id,
                                target_table=subvalue.__class__.__name__.lower() + 's')
                    if changed:
                        self._inserts[table_name].append(link)
            elif isinstance(trait, MetInstance):
                logger.debug("Getting saved data for MetInstance trait")
                value = getattr(obj, tname)
                #logger.debug("MetInstance type value to store: %s", value)
                # handle a sub-object
                # if it is not assigned, use and empty unique_id
                if value is None:
                    state[tname] = ''
                # otherwise, store the uid and allow the object to store
                # itself
                else:
                    logger.debug("Getting saved data")
                    state[tname] = value.unique_id
                    logger.debug("Adding saved MetInstance data for subvalue to Workspace.self")
                    self._get_save_data(value, override)
            elif changed:
                logger.debug("Storing raw value of changed data")
                value = getattr(obj, tname)
                # store the raw value in this table
                state[tname] = value
        if prev_uid and changed:
            self._updates[name].append((obj.unique_id, obj.prev_uid))
        else:
            state['prev_uid'] = ''
        if changed:
            self._inserts[name].append(state)

        logger.debug("Exiting Workspace._get_save_data")

    def fix_table(self, table_name):
        """Fix a table by converting floating point values to doubles"""
        klass = self.subclass_lut.get(table_name, None)
        if not klass:
            return
        table_name = self.tablename_lut[klass]
        for (tname, trait) in klass.class_traits().items():
            if isinstance(trait, MetFloat):
                self.convert_to_double(table_name, tname)

    def retrieve(self, object_type, **kwargs):
        """Retrieve an object from the database."""
        object_type = object_type.lower()
        logger.debug("Retrieving %s object type", object_type)
        logger.debug("Running retrieve with kwargs: %s", kwargs)
        klass = self.subclass_lut.get(object_type, None)
        items = []
        db = self.get_connection()
        db.begin()
        try:
            if object_type not in db:
                if not klass:
                    raise ValueError('Unknown object type: %s' % object_type)
                object_type = self.tablename_lut[klass]
            if '_' not in object_type:
                if kwargs.get('username', '') in ['*', 'all']:
                    kwargs.pop('username')
                else:
                    kwargs.setdefault('username', getpass.getuser())
            query = 'select * from `%s` where (' % object_type
            clauses = []
            for (key, value) in kwargs.items():
                if isinstance(value, list):
                    if len(value) == 0:
                        return []
                    clauses.append('%s in ("%s")' % (key, '", "'.join(value)))
                elif not isinstance(value, str):
                    clauses.append("%s = %s" % (key, value))
                elif '%%' in value:
                    clauses.append('%s = "%s"' % (key, value.replace('%%', '%')))
                elif '%' in value:
                    clauses.append('%s like "%s"' % (key, value.replace('*', '%')))
                else:
                    clauses.append('%s = "%s"' % (key, value))
            query += ' and '.join(clauses) + ')'
            if not clauses:
                query = query.replace(' where ()', '')
            try:
                logger.debug("Getting items from db with query: %s", query)
                items = list(db.query(query))
            except Exception as err:
                if 'Unknown column' in str(err):
                    keys = [k for k in klass.class_traits().keys()
                            if not k.startswith('_')]
                    raise ValueError('Invalid column name, valid columns: %s' % keys) from err
                raise err
            items = [klass(**i) for i in items]
            uids = [i.unique_id for i in items]
            if not items:
                return []
            for (tname, trait) in items[0].traits().items():
                if isinstance(trait, MetList):
                    table_name = '_'.join([object_type, tname])
                    if table_name not in db:
                        for i in items:
                            setattr(i, tname, [])
                        continue
                    querystr = 'select * from `%s` where source_id in ("' % table_name
                    querystr += '" , "'.join(uids)
                    result = db.query(querystr + '")')
                    sublist = defaultdict(list)
                    for r in result:
                        stub = Stub(unique_id=r['target_id'],
                                    object_type=r['target_table'])
                        sublist[r['source_id']].append(stub)
                    for i in items:
                        setattr(i, tname, sublist[i.unique_id])
                elif isinstance(trait, MetInstance):
                    pass
            for i in items:
                if not i.prev_uid:
                    i.prev_uid = 'origin'
                i._changed = False
            items.sort(key=lambda x: x.last_modified)
            db.commit()
        except Exception as err:
            rollback_and_log(db, err)
        finally:
            close_db_connection(db)
        return items

    def remove(self, object_type, **kwargs):
        """Remove an object from the database"""
        override = kwargs.pop('_override', False)
        if not override:
            msg = 'Are you sure you want to delete the entries? (Y/N)'
            ans = input(msg)
            if not ans[0].lower().startswith('y'):
                print('Aborting')
                return
        object_type = object_type.lower()
        klass = self.subclass_lut.get(object_type, None)
        if not klass:
            raise ValueError('Unknown object type: %s' % object_type)
        object_type = self.tablename_lut[klass]
        kwargs.setdefault('username', getpass.getuser())
        query = 'delete from `%s` where (' % object_type
        clauses = []
        for (key, value) in kwargs.items():
            if not isinstance(value, str):
                clauses.append("%s = %s" % (key, value))
                continue
            if '%%' in value:
                clauses.append('%s = "%s"' % (key, value.replace('%%', '%')))
            elif '%' in value:
                clauses.append('%s like "%s"' % (key, value.replace('*', '%')))
            else:
                clauses.append('%s = "%s"' % (key, value))
        query += ' and '.join(clauses)
        query += ')'
        if not clauses:
            query = query.replace(' where ()', '')
        db = self.get_connection()
        db.begin()
        try:
            if any([isinstance(i, MetList) for i in klass.class_traits().values()]):
                uid_query = query.replace('delete ', 'select unique_id ')
                uids = [i['unique_id'] for i in db.query(uid_query)]
                sub_query = 'delete from `%s` where source_id in ("%s")'
                for (tname, trait) in klass.class_traits().items():
                    table_name = '%s_%s' % (object_type, tname)
                    if not uids or table_name not in db:
                        continue
                    if isinstance(trait, MetList):
                        table_query = sub_query % (table_name, '", "'.join(uids))
                        try:
                            db.query(table_query)
                        except Exception as e:
                            print(e)
            try:
                db.query(query)
            except Exception as e:
                if 'Unknown column' in str(e):
                    keys = [k for k in klass.class_traits().keys()
                            if not k.startswith('_')]
                    raise ValueError('Invalid column name, valid columns: %s' % keys)
                else:
                    raise e
            print('Removed')
            db.commit()
        except Exception as err:
            rollback_and_log(db, err)
        finally:
            close_db_connection(db)

    def remove_objects(self, objects, all_versions=True, **kwargs):
        """Remove a list of objects from the database."""
        if not isinstance(objects, (list, set)):
            objects = [objects]
        if not objects:
            print('No objects selected')
            return
        override = kwargs.pop('_override', False)
        if not override:
            msg = ('Are you sure you want to delete the %s object(s)? (Y/N)'
                   % len(objects))
            ans = input(msg)
            if not ans[0].lower().startswith('y'):
                print('Aborting')
                return
        ids = defaultdict(list)
        username = getpass.getuser()
        attr = 'head_id' if all_versions else 'unique_id'
        db = self.get_connection()
        db.begin()
        try:
            for obj in objects:
                if not override and obj.username != username:
                    continue
                name = self.tablename_lut[obj.__class__]
                ids[name].append(getattr(obj, attr))
                for (tname, trait) in obj.traits().items():
                    if isinstance(trait, MetList):
                        subname = '%s_%s' % (name, tname)
                        ids[subname].append(getattr(obj, attr))
            for (table_name, uids) in ids.items():
                if table_name not in db:
                    continue
                query = 'delete from `%s` where %s in ("'
                query = query % (table_name, attr)
                query += '" , "'.join(uids)
                query += '")'
                db.query(query)
            print(('Removed %s object(s)' % len(objects)))
            db.commit()
        except Exception as err:
            rollback_and_log(db, err)
        finally:
            close_db_connection(db)

def _get_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__()
                                   for g in _get_subclasses(s)]

def format_timestamp(tstamp):
    """Get a formatted representation of a timestamp."""
    try:
        ts = pd.Timestamp.fromtimestamp(int(tstamp))
        return ts.isoformat()
    except Exception:
        return str(tstamp)


class MetList(List):
    allow_none = True
    def validate(self, obj, value):
#         value = super(MetList, self).validate(obj, value)
        value = super().validate(obj, value)
        value = NotifyList(value)

        #value.register_callback(lambda: setattr(obj, '_changed', True))
        callback = functools.partial(setattr, obj, '_changed', True)
        value.register_callback(callback)
        return value


class MetUnicode(CUnicode):
    allow_none = True


class MetFloat(CFloat):
    allow_none = True


class MetInt(CInt):
    allow_none = True


class MetBool(CBool):
    allow_none = True


class Stub(HasTraits):

    unique_id = MetUnicode()
    object_type = MetUnicode()

    def retrieve(self):
        wsi = Workspace.get_instance()
        return wsi.retrieve(self.object_type, username='*', unique_id=self.unique_id)[0]

    def __repr__(self):
        return '%s %s' % (self.object_type.capitalize(),
                          self.unique_id)

    def __str__(self):
        return str(self.unique_id)


class MetInstance(Instance):
    allow_none = True

    def validate(self, obj, value):
        if isinstance(value, (self.klass, Stub)):
            return value
        elif isinstance(value, six.string_types):
            if value:
                return Stub(unique_id=value,
                            object_type=self.klass.__name__)
            else:
                return None
        else:
            self.error(obj, value)


class MetEnum(Enum):
    allow_none = True


def get_from_nersc(user, relative_path):
    """Load a remote data file from NERSC to an H5 file

    Parameters
    ----------
    user : str
        NERSC user account
    relative_path : str
        Path to file from "/project/projectdirs/metatlas/original_data/<user>/"
    """
    import pexpect
    from IPython.display import clear_output

    cmd = 'scp -o StrictHostKeyChecking=no '
    path = "/project/projectdirs/metatlas/original_data/%s/%s"
    path = path % (user, relative_path)
    cmd += '%s@edisongrid.nersc.gov:%s . && echo "Download Complete"'
    cmd = cmd % (user, path)
    print(cmd)
    proc = pexpect.spawn(cmd)
    proc.expect("assword:*")
    passwd = eval(input())
    clear_output()
    proc.send(passwd)
    proc.send('\r')
    proc.expect('Download Complete')
    proc.close()
    return os.path.abspath(os.path.basename(relative_path))


def rollback_and_log(db_connection, err):
    """
    inputs:
        db_connection: a dataset instance in a transaction that needs to be rolled back
        err: exception instance that ended the transaction
    """
    caller_name = inspect.stack()[1][3]
    try:
        db_connection.rollback()
    except AttributeError:
        logger.error("Cannot rollback transaction as db_connection is None.")
    logger.error("Transaction rollback within %s()", caller_name)
    logger.exception(err)
    raise err


def close_db_connection(db_connection):
    """Close database connection without raising exceptions"""
    try:
        db_connection.close()
    except AttributeError:
        logger.warning("AttributeError while closing database connection -- ignoring.")
