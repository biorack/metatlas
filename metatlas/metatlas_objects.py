import getpass
import uuid
import time
import os
import sys
import pprint
import tables
import numpy as np
import six
from collections import defaultdict
from pwd import getpwuid
from tabulate import tabulate

try:
    from traitlets import (
        HasTraits, CUnicode, List, CInt, Instance, Enum,
        CFloat, TraitError, CBool)
except ImportError:
    from IPython.utils.traitlets import (
        HasTraits, CUnicode, List, CInt, Instance, Enum,
        CFloat, TraitError, CBool)
import dataset


POLARITY = ('positive', 'negative', 'alternating')
NERSC_USER = '/project/projectdirs/metatlas/mysql_user.txt'


# Observable List from
# http://stackoverflow.com/a/13259435

def callback_method(func):
    def notify(self, *args, **kwargs):
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

    # Take care to return a new NotifyList if we slice it.
    if sys.version_info[0] < 3:
        __setslice__ = callback_method(list.__setslice__)
        __delslice__ = callback_method(list.__delslice__)

        def __getslice__(self, *args):
            return self.__class__(list.__getslice__(self, *args))

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


class MetList(List):
    allow_none = True

    def validate(self, obj, value):
        value = super(MetList, self).validate(obj, value)
        value = NotifyList(value)

        def callback(*args):
            obj._notify_trait(self.name, value, value)
        value.register_callback(callback)
        return value


class MetUnicode(CUnicode):
    allow_none = True


class MetFloat(CFloat):
    allow_none = True


class MetInt(CInt):
    allow_none = True


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


class Stub(HasTraits):

    unique_id = MetUnicode()
    object_type = MetUnicode()

    def retrieve(self):
        return retrieve(self.object_type, unique_id=self.unique_id)[0]

    def __repr__(self):
        return '%s %s' % (self.object_type.capitalize(),
                          self.unique_id)


class Workspace(object):

    def __init__(self):
        # allow for fallback when not on NERSC
        if os.path.exists(NERSC_USER):
            with open(NERSC_USER) as fid:
                pw = fid.read().strip()
            self.path = 'mysql+pymysql://meta_atlas_admin:%s@scidb1.nersc.gov/meta_atlas' % pw
        else:
            self.path = 'sqlite:///' + getpass.getuser() + '_workspace.db'
        self._db = None
        # handle circular references
        self.seen = dict()

    @property
    def db(self):
        if self._db:
            return self._db
        self._db = dataset.connect(self.path)
        if 'sqlite' in self.path:
            os.chmod(self.path[10:], 0o775)
        return self.db

    def convert_to_double(self, table, entry):
        self.db.query('alter table `%s` modify `%s` double' % (table, entry))

    def insert(self, name, state):
        name = name.lower()
        self.db.create_table(name, primary_id='unique_id',
                             primary_type='String(32)')
        self.db[name].insert(state)

    def find_one(self, table_name, **kwargs):
        return self.db[table_name.lower()].find_one(**kwargs)

    def find(self, table_name, **kwargs):
        return self.db[table_name.lower()].find(**kwargs)

    def save_objects(self, objects, _override=False):
        if not isinstance(objects, (list, set)):
            objects = [objects]
        self._seen = dict()
        self._link_updates = defaultdict(list)
        self._updates = defaultdict(list)
        self._inserts = defaultdict(list)
        for obj in objects:
            self.save(obj, _override)
        for (table_name, updates) in self._link_updates.items():
            if table_name not in self.db:
                continue
            with self.db:
                for (uid, prev_uid) in updates:
                    self.db.query('update `%s` set source_id = "%s" where source_id = "%s"' % (table_name, prev_uid, uid))
        for (table_name, updates) in self._updates.items():
            if '_' not in table_name and table_name not in self.db:
                self.db.create_table(table_name, primary_id='unique_id',
                                     primary_type='String(32)')
            with self.db:
                for (uid, prev_uid) in updates:
                    self.db.query('update `%s` set unique_id = "%s" where unique_id = "%s"' % (table_name, prev_uid, uid))
        for (table_name, inserts) in self._inserts.items():
            if '_' not in table_name and table_name not in self.db:
                self.db.create_table(table_name, primary_id='unique_id',
                                     primary_type='String(32)')
            self.db[table_name].insert_many(inserts)

    def save(self, obj, override=False):
        if obj.unique_id in self._seen:
            return
        if isinstance(obj, Stub):
            return
        name = TABLENAME_LUT[obj.__class__]
        self._seen[obj.unique_id] = True
        changed, prev_uid = obj._update(override)
        state = dict()
        for (tname, trait) in obj.traits().items():
            if tname.startswith('_'):
                continue
            if isinstance(trait, List):
                # handle a list of objects by using a Link table
                # create the link table if necessary
                table_name = '_'.join([name, tname])
                if changed and prev_uid:
                    self._link_updates[table_name].append((obj.unique_id,
                                                           obj.prev_uid))
                value = getattr(obj, tname)
                # do not store this entry in our own table
                if not value:
                    continue
                # create an entry in the table for each item
                # store the item in its own table
                for subvalue in value:
                    self.save(subvalue, override)
                    link = dict(source_id=obj.unique_id,
                                head_id=obj.head_id,
                                target_id=subvalue.unique_id,
                                target_table=subvalue.__class__.__name__.lower() + 's')
                    if changed:
                        self._inserts[table_name].append(link)
            elif isinstance(trait, MetInstance):
                value = getattr(obj, tname)
                # handle a sub-object
                # if it is not assigned, use and empty unique_id
                if value is None:
                    state[tname] = ''
                # otherwise, store the uid and allow the object to store
                # itself
                else:
                    state[tname] = value.unique_id
                    self.save(value, override)
            elif changed:
                value = getattr(obj, tname)
                # store the raw value in this table
                state[tname] = value
        if prev_uid and changed:
            self._updates[name].append((obj.unique_id, obj.prev_uid))
        else:
            state['prev_uid'] = ''
        if changed:
            self._inserts[name].append(state)

# Singleton Workspace object
workspace = Workspace()


def retrieve(object_type, **kwargs):
    """Get objects from the Metatlas object database.

    Parameters
    ----------
    object_type: string
      The type of object to search for (i.e. "Groups").
    **kwargs
      Specific search queries (i.e. name="Sargasso").
      Use '%' for wildcard patterns (i.e. description='Hello%').
      If you want to match a '%' character, use '%%'.

    Returns
    -------
    objects: list
      List of Metatlas Objects meeting the criteria.  Will return the
      latest version of each object.
    """
    object_type = object_type.lower()
    klass = SUBCLASS_LUT.get(object_type, None)
    if object_type not in workspace.db:
        if not klass:
            raise ValueError('Unknown object type: %s' % object_type)
        object_type = TABLENAME_LUT[klass]
    # Example query if group id is given
    # SELECT *
    # FROM tablename
    # WHERE (city = 'New York' AND name like 'IBM%')

    # Example query where unique id and group id are not given
    # (to avoid getting all versions of the same object)
    # http://stackoverflow.com/a/12102288
    # SELECT *
    # from (SELECT * from `groups`
    #       WHERE (name='spam') ORDER BY last_modified)
    # x GROUP BY head_id
    query = 'select * from `%s` where (' % object_type
    clauses = []
    for (key, value) in kwargs.items():
        if not isinstance(value, six.string_types):
            clauses.append("%s = %s" % (key, value))
            continue
        if '%%' in value:
            clauses.append('%s = "%s"' % (key, value.replace('%%', '%')))
        elif '%' in value:
            clauses.append('%s like "%s"' % (key, value.replace('*', '%')))
        else:
            clauses.append('%s = "%s"' % (key, value))
    if 'unique_id' not in kwargs and klass:
        clauses.append('unique_id = head_id')
    query += ' and '.join(clauses) + ')'
    if not clauses:
        query = query.replace(' where ()', '')
    try:
        items = [i for i in workspace.db.query(query)]
    except Exception as e:
        if 'Unknown column' in str(e):
            keys = [k for k in klass.class_traits().keys()
                    if not k.startswith('_')]
            raise ValueError('Invalid column name, valid columns: %s' % keys)
        else:
            raise(e)
    items = [klass(**i) for i in items]
    uids = [i.unique_id for i in items]
    if not items:
        return []
    # get stubs for each of the list items
    for (tname, trait) in items[0].traits().items():
        if isinstance(trait, List):
            table_name = '_'.join([object_type, tname])
            if table_name not in workspace.db:
                for i in items:
                    setattr(i, tname, [])
                continue
            querystr = 'select * from `%s` where source_id in ("' % table_name
            querystr += '" , "'.join(uids)
            result = workspace.db.query(querystr + '")')
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
    return items


def remove(object_type, **kwargs):
    """Remove objects from the Metatlas object database.

    Parameters
    ----------
    object_type: string
      The type of object to remove (i.e. "Groups").
    **kwargs
      Specific search queries (i.e. name="Sargasso").
      Use '%' for wildcard patterns (i.e. description='Hello%').
      If you want to match a '%' character, use '%%'.
    """
    object_type = object_type.lower()
    klass = SUBCLASS_LUT.get(object_type, None)
    if not klass:
        raise ValueError('Unknown object type: %s' % object_type)
    object_type = TABLENAME_LUT[klass]
    kwargs.setdefault('username', getpass.getuser())
    override = kwargs.pop('_override', False)
    # Example query:
    # DELETE *
    # FROM tablename
    # WHERE (city = 'New York' AND name like 'IBM%')
    query = 'delete from `%s` where (' % object_type
    clauses = []
    for (key, value) in kwargs.items():
        if not isinstance(value, six.string_types):
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
    if not override:
        if sys.version.startswith('2'):
            ans = raw_input('Are you sure you want to delete these entries?')
        else:
            ans = input('Are you sure you want to delete these entries?')
        if ans[0].lower() != 'y':
            return
    # check for lists items that need removal
    if any([isinstance(i, MetList) for i in klass.class_traits().values()]):
        uid_query = query.replace('delete ', 'select unique_id ')
        uids = [i['unique_id'] for i in workspace.db.query(uid_query)]
        sub_query = 'delete from `%s` where source_id in ("%s")'
        for (tname, trait) in klass.class_traits().items():
            table_name = '%s_%s' % (object_type, tname)
            if not uids or table_name not in workspace.db:
                continue
            if isinstance(trait, MetList):
                table_query = sub_query % (table_name, '", "'.join(uids))
                print(table_query)
                try:
                    workspace.db.query(table_query)
                except Exception as e:
                    print(e)
    try:
        workspace.db.query(query)
    except Exception as e:
        if 'Unknown column' in str(e):
            keys = [k for k in klass.class_traits().keys()
                    if not k.startswith('_')]
            raise ValueError('Invalid column name, valid columns: %s' % keys)
        else:
            raise(e)


def remove_objects(objects, all_versions=True, **kwargs):
    """Remove objects from the database.

    Parameters
    ----------
    all_versions: boolean, optional
        If True, remove all versions of the object sharing the current
        head_id.
    """
    if not isinstance(objects, (list, set)):
        objects = [objects]
    ids = defaultdict(list)
    username = getpass.getuser()
    override = kwargs.pop('_override', False)
    attr = 'head_id' if all_versions else 'unique_id'
    for obj in objects:
        if not override and obj.username != username:
            continue
        name = TABLENAME_LUT[obj.__class__]
        ids[name].append(getattr(obj, attr))
        # remove list items as well
        for (tname, trait) in obj.traits().items():
            if isinstance(trait, MetList):
                subname = '%s_%s' % (name, tname)
                ids[subname].append(getattr(obj, attr))
    if not override:
        if sys.version.startswith('2'):
            ans = raw_input('Are you sure you want to delete these entries?')
        else:
            ans = input('Are you sure you want to delete these entries?')
        if ans[0].lower() != 'y':
            return
    for (table_name, uids) in ids.items():
        if table_name not in workspace.db:
            continue
        query = 'delete from `%s` where %s in ("'
        query = query % (table_name, attr)
        query += '" , "'.join(uids)
        query += '")'
        workspace.db.query(query)


def store(objects, **kwargs):
    """Store Metatlas objects in the database.

    Parameters
    ----------
    objects: Metatlas object or list of Metatlas Objects
        Object(s) to store in the database.
    """
    workspace.save_objects(objects, **kwargs)


def format_timestamp(tstamp):
    """Get a formatted representation of a timestamp."""
    from pandas import Timestamp
    ts = Timestamp.fromtimestamp(int(tstamp))
    return ts.isoformat()


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
        help_text = trait.get_metadata('help')
        if not help_text:
            help_text = '%s value.' % tname
        help_text = help_text.strip()
        if help_text.endswith('.'):
            help_text = help_text[:-1]
        if trait.get_metadata('readonly'):
            help_text += ' (read only)'
        help_text += '.'
        doc += '    %s\n' % help_text
    cls.__doc__ = doc
    return cls


@set_docstring
class MetatlasObject(HasTraits):

    name = MetUnicode('Untitled', help='Name of the object')
    description = MetUnicode('No description',
                             help='Description of the object')
    unique_id = MetUnicode(help='Unique identifier for the object',
                           readonly=True)
    creation_time = MetInt(help='Unix timestamp at object creation',
                           readonly=True)
    username = MetUnicode(help='Username who created the object',
                          readonly=True)
    last_modified = MetInt(help='Unix timestamp at last object update',
                           readonly=True)
    prev_uid = MetUnicode(help='Unique id of previous version', readonly=True)
    head_id = MetUnicode(help='Unique id of most recent version of this object', readonly=True)
    _loopback_guard = CBool(False, readonly=True)
    _changed = CBool(False, readonly=True)

    def __init__(self, **kwargs):
        """Set the default attributes."""
        kwargs.setdefault('unique_id', uuid.uuid4().hex)
        kwargs.setdefault('head_id', kwargs['unique_id'])
        kwargs.setdefault('username', getpass.getuser())
        kwargs.setdefault('creation_time', int(time.time()))
        kwargs.setdefault('last_modified', int(time.time()))
        super(MetatlasObject, self).__init__(**kwargs)
        self._changed = True
        self.on_trait_change(self._on_update)

    def _update(self, override_user=False):
        """Store the object in the workspace, including child objects.

        Child objects are stored in their own tables, and Lists of
        Child objects are captured in link tables.
        """
        self._loopback_guard = True
        # see if we need to create a new object
        if not override_user and self.username != getpass.getuser():
            self._changed = True
            self.prev_uid = self.unique_id
            self.unique_id = uuid.uuid4().hex
            self.head_id = self.unique_id
            self.username = getpass.getuser()
            self.last_modified = time.time()
            return True, None
        else:
            changed, prev_uid = self._changed, self.prev_uid
            if changed:
                self.last_modified = time.time()
                self.prev_uid = uuid.uuid4().hex
        self._changed = False
        self._loopback_guard = False
        return changed, prev_uid

    def clone(self, recursive=False):
        """Create a new version of this object.

        Parameters
        ----------
        recursive: boolean, optional
            If true, clone all of the descendant objects as well.

        Returns
        -------
        obj: MetatlasObject
            Cloned object.
        """
        obj = self.__class__()
        for (tname, trait) in self.traits().items():
            if tname.startswith('_') or trait.get_metadata('readonly'):
                continue
            val = getattr(self, tname)
            if recursive and isinstance(trait, List):
                val = [v.clone(True) for v in val]
            elif recursive and isinstance(trait, MetInstance) and val:
                val = val.clone(True)
            setattr(obj, tname, val)
        obj.prev_uid = self.unique_id
        obj.head_id = self.unique_id
        obj.unique_id = uuid.uuid4().hex
        return obj

    def edit(self):
        """Create an editor for the object

        Creates labels and editors for all traits.

        If the trait is marked readyonly or is an instance of another Trait,
        it will be disabled.
        If the trait is an Enum, there will be a dropdown.
        Otherwise, there will be a text editor.
        """
        edit_traits(self)

    def show_diff(self, unique_id=None):
        """Show a diff of what has changed between this and previous version.

        Parameters
        ----------
        unique_id: optional, string
            Unique id to compare against (defaults to current entry in db).
        """
        if unique_id is None:
            unique_id = self.unique_id
        obj = retrieve(self.__class__.__name__, unique_id=unique_id)
        if len(obj) != 1:
            print('No change!')
            return
        obj = obj[0]
        msg = []
        for (tname, trait) in self.traits().items():
            if tname.startswith('_') or trait.get_metadata('readonly'):
                continue
            val = getattr(self, tname)
            other = getattr(obj, tname)
            if isinstance(trait, MetInstance):
                if val.unique_id != other.unique_id:
                    msg.append((tname, other.unique_id, val.unique_id))
            elif isinstance(trait, MetList):
                if not len(val) == len(other):
                    msg.append((tname, '%s items' % len(other),
                                '%s items' % len(val)))
                else:
                    for (v, o) in zip(val, other):
                        if not v.unique_id == o.unique_id:
                            msg.append((tname, 'objects changed', ''))
                            break
            elif val != other:
                msg.append((tname, str(other), str(val)))
        print(tabulate(msg))

    def _on_update(self, name):
        """When the model changes, set the update fields.
        """
        if self._loopback_guard or name.startswith('_'):
            return
        self._changed = True

    def __str__(self):
        return self.name + ' (%s)' % self.unique_id

    def __repr__(self):
        names = sorted(self.trait_names())
        names.remove('name')
        names = ['name'] + [n for n in names if not n.startswith('_')]
        state = dict([(n, getattr(self, n)) for n in names])
        state['creation_time'] = format_timestamp(self.creation_time)
        state['last_modified'] = format_timestamp(self.last_modified)
        return pprint.pformat(state)

    def __getattribute__(self, name):
        """Automatically resolve stubs on demand.
        """
        value = super(MetatlasObject, self).__getattribute__(name)
        if isinstance(value, Stub):
            value = value.retrieve()
            setattr(self, name, value)
        elif isinstance(value, list) and value:
            new = []
            changed = False
            for subvalue in value:
                if isinstance(subvalue, Stub):
                    new.append(subvalue.retrieve())
                    changed = True
                else:
                    new.append(subvalue)
            if changed:
                setattr(self, name, new)
                value = new
        return value


def _get_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__()
                                   for g in _get_subclasses(s)]


@set_docstring
class Method(MetatlasObject):
    """
    For each LCMS run, a Method is a consistent description of
    how the sample was prepared and LCMS data was collected.
    """
    protocol_ref = MetUnicode(help='Reference to a published protocol: ' +
        'identical to the protocol used.')
    quenching_method = MetUnicode(help='Description of the method used to ' +
        'stop metabolism.')
    extraction_solvent = MetUnicode(help='Solvent or solvent mixture used to ' +
        'extract metabolites.')
    reconstitution_method = MetUnicode(help='Solvent or solvent mixture ' +
        'the extract is reconstituted in prior to injection for an LCMS Run.')
    mobile_phase_a = MetUnicode(help='Solvent or solvent mixture.')
    mobile_phase_b = MetUnicode(help='Solvent or solvent mixture.')
    temporal_parameters = MetUnicode('List of temporal changes to the' +
        'mixing of mobile_phase_a and mobile_phase_b.')
    column_model = MetUnicode(help='Brand and catalog number of the column')
    column_type = MetUnicode(help='Class of column used.')
    scan_mz_range = MetUnicode(help='Minimum and ' +
        'maximum mz recorded for a run.')
    instrument = MetUnicode(help='Brand and catalog number for the ' +
        'mass spectrometer.')
    ion_source = MetUnicode(help='Method for ionization.')
    mass_analyzer = MetUnicode(help='Method for detection.')
    polarity = MetEnum(POLARITY, 'positive', help='polarity for the run')


@set_docstring
class Sample(MetatlasObject):
    """A Sample is the material that is processed with a Method."""
    pass


def load_lcms_files(mzml_files):
    """Parse mzML files and load them into LcmsRun objects.

    Note: This should be done automatically for runs in
    /project/projectdirs/metatlas/raw_data/<username>

    Parameters
    ----------
    mzml_files: list of str
       List of paths to mzml_files.

    Returns
    -------
    runs: list
       List of LcmsRun objects.
    """
    runs = []
    for fname in mzml_files:
        hdf5_file = fname.replace('.mzML', '.h5')
        if os.path.exists(hdf5_file):
            print('File already exists: %s' % hdf5_file)
            continue
        try:
            from metatlas import LcmsRun, mzml_to_hdf
            hdf_file = mzml_to_hdf(fname)
            user = getpwuid(os.stat(fname).st_uid).pw_name
            filename = os.path.splitext(os.path.basename(fname))[0]
            dirname = os.path.dirname(fname)
            experiment = os.path.basename(dirname)
            description = experiment + ' ' + filename
            ctime = os.stat(fname).st_ctime
            run = LcmsRun(name=filename, description=description,
                          created_by=user,
                          modified_by=user,
                          created=ctime, last_modified=ctime,
                          mzml_file=fname, hdf5_file=hdf_file)
            runs.append(run)
        except Exception as e:
            print(e)
    store(runs)
    return runs


@set_docstring
class LcmsRun(MetatlasObject):
    """An LCMS run is the reference to a file prepared with liquid 
    chromatography and mass spectrometry.

    The msconvert program is used to convert raw data (centroided is prefered)
    to mzML.

    Note: These objects are not intented to be created directly, but by putting
    the files in /project/projectdirs/metatlas/raw_data/<username> or by 
    running `load_lcms_files()`.
    """
    method = MetInstance(Method)
    hdf5_file = MetUnicode(help='Path to the HDF5 file at NERSC')
    mzml_file = MetUnicode(help='Path to the MZML file at NERSC')
    sample = MetInstance(Sample)

    def interact(self, min_mz=None, max_mz=None, polarity=None, ms_level=1):
        """Interact with LCMS data - XIC linked to a Spectrogram plot.

        Parameters
        ----------
        min_mz: float
            Minimum m/z (defaults to file min)
        max_mz: float
            Maximum m/z (defaults to file max)
        polarity: {0, 1}
            Polarity (defaults to neg. if present in file, else pos.)
        ms_level: {0, 1}
            The ms level.
        """
        import matplotlib.pyplot as plt
        from metatlas import get_XIC, get_spectrogram, get_info
        fid = tables.open_file(self.hdf5_file)

        info = get_info(fid)
        if polarity is None:
            if info['ms%s_neg' % ms_level]['nrows']:
                polarity = 0
            else:
                polarity = 1
        if polarity == 0:
            table_name = 'ms%s_neg' % ms_level
        else:
            table_name = 'ms%s_pos' % ms_level
        if min_mz is None:
            min_mz = info[table_name]['min_mz']
        if max_mz is None:
            max_mz = info[table_name]['max_mz']

        rt, irt = get_XIC(fid, min_mz, max_mz, 1, polarity)
        mz, imz = get_spectrogram(fid, rt[0], rt[1], 1, polarity)

        fig, (ax1, ax2) = plt.subplots(ncols=2)
        ax1.plot(rt, irt)
        ax1.set_title('XIC: %0.1f - %0.1f m/z' % (min_mz, max_mz))
        ax1.set_xlabel('Time (min)')
        ax1.set_ylabel('Intensity')
        ax1._vline = ax1.axvline(rt[0])

        ax2.vlines(mz, 0, imz)
        ax2.set_xlabel('Mass (m/z)')
        ax2.set_title('MS%s Spectrogram at %0.1f min' % (ms_level, rt.min()))

        def callback(event):
            if event.inaxes == ax1:
                rt_event = event.xdata
                # get the closest actual RT
                idx = (np.abs(rt - rt_event)).argmin()
                mz, imz = get_spectrogram(fid, rt[idx], rt[idx], 1, polarity)

                ax1._vline.remove()
                ax1._vline = ax1.axvline(rt_event, color='k')

                ax2.clear()
                ax2.vlines(mz, 0, imz)
                ax2.set_xlabel('Mass (m/z)')
                ax2.set_title('Spectrogram at %0.1f min' % rt_event)
                fig.canvas.draw()
        fig.canvas.mpl_connect('button_press_event', callback)
        fig.canvas.draw()


@set_docstring
class ReferenceDatabase(MetatlasObject):
    """External databases (PubChem, KEGG, MetaCyc, KBase, etc)."""
    enabled = CBool(True)


@set_docstring
class FunctionalSet(MetatlasObject):
    """Functional sets of compounds.  For example set called "hexose" would
    include "glucose, galactose, etc".  Functional sets can be sets-of-sets.
    "Sugars" would be a set that contains "Hexoses".
    """
    enabled = CBool(True)
    members = MetList(MetInstance(MetatlasObject))


@set_docstring
class Compound(MetatlasObject):
    """A Compound as a name and description to capture incomplete IDs:
    for example "hexose", "32:0 diacylglycerol".
    For IDs that have high enough confidence for structural assignments an
    InChi string is the ID.
    """
    InChI = MetUnicode(help='IUPAC International Chemical Identifier, optional')
    formula = MetUnicode()
    MonoIsotopic_molecular_weight = MetFloat()
    synonyms = MetUnicode()
    url = MetUnicode(help='Reference database table url')
    reference_xrefs = MetList(MetInstance(ReferenceDatabase),
                           help='Tag a compound with compound ids from ' +
                                'external databases')
    functional_sets = MetList(MetInstance(FunctionalSet))


@set_docstring
class Reference(MetatlasObject):
    """Place holder for future reference sources.
    We expect many in silico methods will soon be robust enough to suggest
    retention times, m/z, and fragmentation.
    MIDAS is a great example of this.
    """
    lcms_run = MetInstance(LcmsRun)
    enabled = CBool(True)
    ref_type = MetUnicode(help='The type of reference')


@set_docstring
class IdentificationGrade(MetatlasObject):
    """
    Each CompoundIdentification will have an identification_grade
    Identification Grades:
    1) High intensity and verifiable by MSMS and RT
    authentic standard
    2) Verifiable by MSMS from database or
    publication
    3) Has fragment ion or neutral loss characteristic
    of a class of compounds
    4) definitive chemical formula and adduct
    5) Significant changing metabolite with MSMS
    suggestion from MIDAS
    6) Significant changing metabolite
    7) Not Significant changing metabolite with
    MSMS suggestion from MIDAS
    8) Not Significant changing metabolite
    """
    pass


ID_GRADES = dict()


class _IdGradeTrait(MetInstance):

    klass = IdentificationGrade

    def validate(self, obj, value):
        global ID_GRADES
        if not value:
            return
        if isinstance(value, self.klass):
            return value
        elif isinstance(value, str):
            if value.upper() in ID_GRADES:
                return ID_GRADES[name.upper()]
            objects = retrieve('identificationgrade', name=value.upper())
            if objects:
                ID_GRADES[value.upper()] = objects[-1]
                return objects[-1]
            else:
                self.error(obj, value)
        else:
            self.error(obj, value)


@set_docstring
class CompoundIdentification(MetatlasObject):
    """A CompoundIdentification links multiple sources of evidence about a 
    compound's identity to an Atlas."""
    compound = MetList(MetInstance(Compound))
    identification_grade = _IdGradeTrait(
        help='Identification grade of the id (can be specified by a letter A-H'
    )
    references = MetList(MetInstance(Reference))

    def select_by_type(self, ref_type):
        """Select references by type.

        Parameters
        ----------
        ref_type: {'mz', 'rt', 'fragmentation'}
          The type of reference.
        """
        if ref_type.lower() in ['mz', 'm/z']:
            return [r for r in self.references if isinstance(r, MzReference)]
        elif ref_type.lower() in ['rt', 'retention_time', 'retention time']:
            return [r for r in self.references if isinstance(r, RtReference)]
        elif ref_type.lower() in ['frag', 'fragmentation']:
            return [r for r in self.references if
                    isinstance(r, FragmentationReference)]
        else:
            raise ValueError('Invalid reference type')


@set_docstring
class Atlas(MetatlasObject):
    """An atlas contains many compound_ids."""
    compound_identifications = MetList(
        MetInstance(CompoundIdentification),
        help='List of Compound Identification objects')


@set_docstring
class Group(MetatlasObject):
    """A Group can be a:
    file,
    group of files, or
    group of groups
    """
    items = MetList(MetInstance(MetatlasObject),
                 help='Can contain other groups or LCMS Runs')


@set_docstring
class MzIntensityPair(MetatlasObject):
    mz = MetFloat()
    intensity = MetFloat()


@set_docstring
class FragmentationReference(Reference):
    polarity = MetEnum(POLARITY, 'positive')
    precursor_mz = MetFloat()
    mz_intensities = MetList(MetInstance(MzIntensityPair),
                          help='list of [mz, intesity] tuples that describe ' +
                               ' a fragmentation spectra')


@set_docstring
class RtReference(Reference):
    RTpeak = MetFloat()
    RTmin = MetFloat()
    RTmax = MetFloat()
    RTUnits = MetEnum(('sec', 'min'), 'sec')


@set_docstring
class MzReference(Reference):
    """Source of the assertion that a compound has a given m/z and
    other properties directly tied to m/z.
    """ 
    mz = MetFloat()
    mz_tolerance = MetFloat()
    mz_tolerance_units = MetEnum(('ppm', 'Da'), 'ppm')
    detected_polarity = MetEnum(POLARITY, 'positive')
    adduct = MetUnicode(help='Optional adduct')
    modification = MetUnicode(help='Optional modification')
    observed_formula = MetUnicode(help='Optional observed formula')


def find_invalid_runs(**kwargs):
    """Find invalid runs.
    """
    override = kwargs.pop('_override', False)
    if not override:
        kwargs.setdefault('username', getpass.getuser())
    all_runs = retrieve('lcmsruns', **kwargs)
    invalid = []
    for run in all_runs:
        if (not os.path.exists(run.hdf5_file) or
            not os.path.exists(run.mzml_file) or
            not os.stat(run.hdf5_file).st_size):
                invalid.append(run)
    return invalid


SUBCLASS_LUT = dict()
TABLENAME_LUT = dict()
for klass in _get_subclasses(MetatlasObject):
    name = klass.__name__.lower()
    SUBCLASS_LUT[name] = klass
    if name.endswith('s'):
        SUBCLASS_LUT[name + 'es'] = klass
        TABLENAME_LUT[klass] = name + 'es'
    else:
        SUBCLASS_LUT[name + 's'] = klass
        TABLENAME_LUT[klass] = name + 's'


def edit_traits(obj):
    """Create an IPython widget editor for a Traits object"""
    try:
        from ipywidgets import Text, Dropdown, HBox, VBox
    except ImportError:
        from IPython.html.widgets import Text, Dropdown, HBox, VBox
    from IPython.display import display
    names = sorted(obj.trait_names())
    names.remove('name')
    names = ['name'] + [n for n in names if not n.startswith('_')]
    items = [Text('', disabled=True)]
    for name in names:
        if name.startswith('_'):
            continue
        try:
            value = getattr(obj, name)
        except TraitError:
            value = None
        trait = obj.traits()[name]
        if name in ['created', 'last_modified']:
            value = format_timestamp(value)
        if (trait.get_metadata('readonly') or
                isinstance(trait, Instance) or value is None):
            items.append(Text(str(value), disabled=True))

        elif isinstance(trait, Enum):
            # create a closure around "name" for the on_trait_change
            # callback
            def create_dropdown(name):
                dd = Dropdown(value=value, options=trait.values)

                def callback(dummy, value):
                    setattr(obj, name, value)
                dd.on_trait_change(callback, 'value')
                items.append(dd)

            create_dropdown(name)
        else:
            # create a closure around "name for the on_trait_change
            # callback
            def create_textbox(name):
                textbox = Text(str(value))

                def callback(dummy, value):
                    try:
                        setattr(obj, name, value)
                    except Exception:
                        textbox.color = 'red'
                textbox.on_trait_change(callback, 'value')
                items.append(textbox)

            create_textbox(name)

    labels = [Text(name, disabled=True) for name in names]
    labels = [Text(obj.__class__.__name__, disabled=True)] + labels
    display(HBox(children=[VBox(children=labels), VBox(children=items)]))


if __name__ == '__main__':
    m1 = Group(name='spam')
    store(m1)
    m1.description = 'baz'
    store(m1)
    print(retrieve('group', name='spam'))
