""" mzml to h5 file conversion """

import argparse
import fcntl
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

from datetime import datetime
from subprocess import check_output

from metatlas.datastructures.metatlas_objects import LcmsRun, store, retrieve
from metatlas.io.mzml_loader import mzml_to_hdf
from metatlas.io.mzml_loader import VERSION_TIMESTAMP
from metatlas.io.system_utils import send_mail


readonly_files = {}  # username (or uid) | a set of files associated with them
other_errors = {}  # info with user | list of error messages
patt = re.compile(r".+\/raw_data\/(?P<username>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")

MEMBERS_CMD = "getent group metatlas | cut -d: -f4"
USERS = sorted(subprocess.check_output(MEMBERS_CMD, shell=True, text=True).strip().split(','))


def move_file(src, dest):
    """ move file and create directories if needed """
    assert os.path.isfile(src)
    dest_dir = dest if os.path.isdir(dest) else os.path.dirname(dest)
    os.makedirs(dest_dir, exist_ok=True)
    shutil.move(src, dest)


def _file_name_to_username(file_name):
    initials_field = os.path.basename(file_name).split("_")[1].lower()
    initials = initials_field.split("-")[-1]
    for user in USERS:
        if user.startswith(initials):
            return user
    return None


def _write(message, out_fh):
    out_fh.write(message + "\n")
    out_fh.flush()


def _write_stdout(message):
    _write(message, sys.stdout)


def _write_stderr(message):
    _write(message, sys.stderr)


def get_acqtime_from_mzml(mzml_file):
    start_time = None
    with open(mzml_file, "r", encoding="utf-8") as mzml:
        for line in mzml:
            if 'start_time' in line:
                start_time = line.split('start_time="')[1].split('"')[0].replace('T', ' ').rstrip('Z')
                break
    if start_time is not None and '-infinity' not in start_time:
        date_object = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
        utc_timestamp = int(time.mktime(date_object.timetuple()))
    else:
        utc_timestamp = int(0)
    return utc_timestamp


def convert(ind, fname):
    """Helper function, converts a single file"""
    _write_stdout(f"({ind+1}): {fname}")

    # Get relevant information about the file.
    username = _file_name_to_username(fname)
    info = patt.match(os.path.abspath(fname))
    if info:
        info = info.groupdict()
    else:
        _write_stdout(f"Invalid path name: {fname}")
        return
    dirname = os.path.dirname(fname)

    # Convert to HDF and store the entry in the database.
    try:
        hdf5_file = fname.replace('mzML', 'h5')
        _write_stderr(f"hdf5file is: {hdf5_file}")
        acquisition_time = get_acqtime_from_mzml(fname)
        mzml_to_hdf(fname, hdf5_file, True)
        os.chmod(hdf5_file, 0o640)
        description = info['experiment'] + ' ' + info['path']
        ctime = os.stat(fname).st_ctime
        # Add this to the database unless it is already there
        try:
            runs = retrieve('lcmsrun', username='*', mzml_file=fname)
        except Exception:
            runs = []
        if not runs:
            run = LcmsRun(name=info['path'],
                          description=description,
                          username=username,
                          experiment=info['experiment'],
                          creation_time=ctime,
                          last_modified=ctime,
                          mzml_file=fname,
                          hdf5_file=hdf5_file,
                          acquisition_time=acquisition_time)
            store(run)
    except Exception as e:
        if 'exists but it can not be written' in str(e):
            if username not in readonly_files:
                readonly_files[username] = set()
            readonly_files[username].add(dirname)
        else:
            msg = traceback.format_exception(*sys.exc_info())
            msg.insert(0, f"Cannot convert {fname}")
            dat = username
            if dat not in other_errors:
                other_errors[username] = []
            other_errors[username].append('\n'.join(msg))
            fail_path = fname.replace('raw_data', 'conversion_failures')
            _write_stderr(f"Moving file\n{fname}\nto\n{fail_path}\n")
            move_file(fname, fail_path)
        _write_stderr(str(e))
        try:
            os.remove(hdf5_file)
        except:
            pass


def update_metatlas(directory):
    """
    Converts all files to HDF in metatlas. Emails the user if there was
    any kind of error with converting a file.
    """
    mzml_files = check_output(f'find {directory} -name "*.mzML"', shell=True)
    mzml_files = mzml_files.decode('utf-8').splitlines()

    # Find valid h5 files newer than the format version timestamp.
    delta = int((time.time() - VERSION_TIMESTAMP) / 60)
    check = f'find {directory} -name "*.h5" -mmin -{delta} -size +2k'
    valid_files = set(check_output(check, shell=True).decode('utf-8').splitlines())
    new_files = [file for file in mzml_files if file.replace('.mzML', '.h5') not in valid_files]

    if new_files:
        _write_stdout(f"Found {len(new_files)} files")
        for ind, ffff in enumerate(new_files):
            convert(ind, ffff)
        if readonly_files:
            for (username, dirnames) in readonly_files.items():
                body = ("Please log in to NERSC and run 'chmod g+rwXs' on the "
                        "following directories:\n%s" % ('\n'.join(dirnames)))
                send_mail('Metatlas Files are Inaccessible', username, body)
        if other_errors:
            for (username, errors) in other_errors.items():
                body = ('Errored files found while loading in Metatlas files:\n\n%s' %
                        '\n********************************\n'.join(errors))
                send_mail('Errors loading Metatlas files', username, body)
    _write_stdout('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")
    args = parser.parse_args()
    _write_stdout(str(args))
    update_metatlas(args.directory[0])
