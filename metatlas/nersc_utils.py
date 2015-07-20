"""Selected functions from https://github.com/kbase/transform/blob/master/lib/biokbase/Transform/script_utils.py
For use on Nersc Ipython instance.
"""
import sys
import os
import dataset
import uuid
import pickle

SHOCK_URL = None
NERSC_WORKSPACE = '/project/projectdirs/metatlas/workspace'


class Workspace(object):

    def __init__(self):
        path = os.environ['USER'] + '.db'
        # allow for fallback when not on NERC
        if os.path.exists(NERSC_WORKSPACE):
            path = os.path.join(NERSC_WORKSPACE, path)
        self.db = dataset.connect('sqlite:///%s' % path)
        os.chmod(path, 0o775)

    def get_object_uid(self, name):
        objects = self.db[name].all()
        if objects:
            return list(objects)[-1]['objid']
        else:
            raise ValueError('No object found called "%s"' % name)

    def save_object(self, obj):
        objid = uuid.uuid4().hex
        self.db['objids'].insert(dict(name=obj['name'], objid=objid))

        obj['objid'] = objid
        obj['data'] = pickle.dumps(obj['data'])
        self.db[obj['name']].insert(obj)

    def get_object_from_ref(self, objid):
        obj = self.db['objids'].find_one(objid=objid)
        if obj:
            obj = self.db[obj['name']].find_one(objid=objid)
            return pickle.loads(obj['data'])

    def get_object(self, name):
        objects = list(self.db[name].all())
        if objects:
            obj = list(objects)[-1]
            return pickle.loads(obj['data'])
        else:
            raise ValueError('No object found called "%s"' % name)

# Singleton workspace object
WORKSPACE = Workspace()


def save_ws_object(obj):
    """Save an object to the workspace

    Parameters
    ----------
    obj : dict
        Object with the fields: type, data, and name.
        The type must be the full typespec
        (e.g. 'MetaboliteAtlas.Compound-0.3')

    Returns
    -------
    id : str
        Object workspace id
    """
    return WORKSPACE.save_object(obj)


def get_object_uid(name):
    return WORKSPACE.get_object_uid(name)


def get_object(name):
    return WORKSPACE.get_object(name)


def get_object_from_ref(ref):
    return WORKSPACE.get_object_from_ref(ref)


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
    if sys.version.startswith('3'):
        passwd = input()
    else:
        passwd = raw_input()
    clear_output()
    proc.send(passwd)
    proc.send('\r')
    proc.expect('Download Complete')
    proc.close()
    return os.path.abspath(os.path.basename(relative_path))
