import os
import pwd
import re
import shutil
import sys
import time
import tables
from collections import defaultdict
from subprocess import Popen, PIPE
from datetime import datetime, time as dtime

from metatlas.mzml_loader import FORMAT_VERSION
from metatlas import LcmsRun, mzml_to_hdf, store, retrieve

ADMIN = 'silvest'


def send_mail(subject, username, body, force=False):
    """Send the mail only once per day."""
    now = datetime.now()
    if force or dtime(00, 00) <= now.time() <= dtime(00, 10):
        for name in [username, ADMIN]:
            msg = 'mail -s "%s" %s@nersc.gov <<< "%s"' % (subject, name, body)
            sys.stdout.write(msg + '\n')
            sys.stdout.flush()
            p = Popen(["bash"], stdin=PIPE)
            p.communicate(msg)


def check_file_validity(hdf5_file):
    # Check for invalid file
    try:
        fid = tables.open_file(hdf5_file)
    except Exception:
        return False
    # Check for current version.
    try:
        format_version = fid.get_node_attr('/', 'format_version')
    except AttributeError:
        return False
    finally:
        fid.close()

    if format_version.decode('utf-8') != FORMAT_VERSION:
        return False
    else:
        return True


def update_metatlas(directory):
    if directory.endswith('/raw_data'):
        return

    new_files = []
    readonly_files = defaultdict(set)
    other_errors = defaultdict(list)
    for root, directories, filenames in os.walk(directory):
        for fname in filenames:
            if fname.endswith('.mzML'):
                fname = os.path.join(root, fname)
                hdf5_file = fname.replace('.mzML', '.h5')
                if not os.path.exists(hdf5_file):
                    new_files.append(fname)
                else:
                    if not check_file_validity(hdf5_file):
                        new_files.append(fname)

    patt = re.compile(r".+\/raw_data\/(?P<username>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")

    sys.stdout.write('Found %s files\n' % len(new_files))
    sys.stdout.flush()

    for (ind, fname) in enumerate(new_files):
        info = patt.match(os.path.abspath(fname))
        if info:
            info = info.groupdict()
        else:
            sys.stdout.write("Invalid path name: %s\n" % fname)
            sys.stdout.flush()
            continue
        dirname = os.path.dirname(fname)
        try:
            username = pwd.getpwuid(os.stat(fname).st_uid).pw_name
        except OSError:
            try:
                username = pwd.getpwuid(os.stat(dirname).st_uid).pw_name
            except Exception:
                username = info['username']

        sys.stdout.write('(%s of %s): %s\n' % (ind + 1, len(new_files), fname))
        sys.stdout.flush()

        try:
            os.chmod(fname, 0o660)
        except Exception as e:
            sys.stderr.write(str(e) + '\n')
            sys.stderr.flush()

        # copy the original file to a pasteur backup
        if os.environ['USER'] == 'pasteur':
            pasteur_path = fname.replace('raw_data', 'pasteur_backup')
            dname = os.path.dirname(pasteur_path)
            if not os.path.exists(dname):
                os.makedirs(dname)
            try:
                shutil.copy(fname, pasteur_path)
            except IOError as e:
                readonly_files[username].add(dirname)
                continue

        # Check if another process already converted this file.
        hdf5_file = fname.replace('.mzML', '.h5')
        if os.path.exists(hdf5_file):
            if check_file_validity(hdf5_file):
                msg = '%s was already converted by another process \n' % fname
                sys.stderr.write(msg)
                sys.stderr.flush()
                continue
            if time.time() - os.stat(hdf5_file).st_mtime < 100:
                msg = '%s is being converted by another process \n' % fname
                sys.stderr.write(msg)
                sys.stderr.flush()
                continue

        # convert to HDF and store the entry in the database
        try:
            mzml_to_hdf(fname, hdf5_file, True)
            os.chmod(hdf5_file, 0o660)
            description = info['experiment'] + ' ' + info['path']
            ctime = os.stat(fname).st_ctime
            # Add this to the database unless it is already there
            runs = retrieve('lcmsrun', username='*', mzml_file=fname)
            if not len(runs):
                run = LcmsRun(name=info['path'], description=description,
                              username=info['username'],
                              experiment=info['experiment'],
                              creation_time=ctime, last_modified=ctime,
                              mzml_file=fname, hdf5_file=hdf5_file)
                store(run)
        except Exception as e:
            if 'exists but it can not be written' in str(e):
                readonly_files[username].add(dirname)
            else:
                other_errors[info['username']].append(str(e))
            sys.stderr.write(str(e) + '\n')
            sys.stderr.flush()

    from metatlas.metatlas_objects import find_invalid_runs
    invalid_runs = find_invalid_runs(_override=True)

    if readonly_files:
        for (username, dirnames) in readonly_files.items():
            body = ("Please log in to NERSC and run 'chmod 777' on the "
                   "following directories:\n%s" % ('\n'.join(dirnames)))
            send_mail('Metatlas Files are Inaccessible', username, body)
    if invalid_runs:
        grouped = defaultdict(list)
        for run in invalid_runs:
            grouped[run.username].append(run.mzml_file)
        for (username, filenames) in grouped.items():
            body = 'You have runs that are not longer accessible\n'
            body += 'To remove them from the database, run the following on ipython.nersc.gov:\n\n'
            body += 'from metatlas.metatlas_objects import find_invalid_runs, remove_objects\n'
            body += 'remove_objects(find_invalid_runs())\n\n'
            body += 'The invalid runs are:\n%s' % ('\n'.join(filenames))
            send_mail('Metatlas Runs are Inaccessible', username, body)
    if other_errors:
        for (username, errors) in other_errors.items():
            body = 'Errored files found while loading in Metatlas files:\n%s' % '\n'.join(errors)
            send_mail('Errors loading Metatlas files', username, body)
    sys.stdout.write('Done!\n')
    sys.stdout.flush()

    send_mail('Run complete', ADMIN, directory, True)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")

    args = parser.parse_args()
    sys.stdout.write(str(args) + '\n')
    sys.stdout.flush()
    update_metatlas(args.directory[0])
