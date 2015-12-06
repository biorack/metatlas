import os
import pwd
import re
import shutil
from collections import defaultdict
from subprocess import Popen, PIPE
from datetime import datetime, time


def send_mail(subject, username, body):
    """Send the mail only once per day."""
    #now = datetime.now()
    #if time(00, 00) <= now.time() <= time(00, 10):
    #    # send it to silvest for now
    body += '\nwas %s' % username
    username = 'silvest'
    msg = 'mail -s "%s" %s@nersc.gov <<< "%s"' % (subject, username, body)
    p = Popen(["bash"], stdin=PIPE)
    p.communicate(msg)


def update_metatlas(directory):

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

    patt = re.compile(r".+\/raw_data\/(?P<username>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")

    print('Found %s new files' % len(new_files))

    for fname in new_files:
        info = patt.match(os.path.abspath(fname))
        if info:
            info = info.groupdict()
        else:
            print("Invalid path name", fname)
            continue
        dirname = os.path.dirname(fname)
        try:
            username = pwd.getpwuid(os.stat(fname).st_uid).pw_name
        except OSError:
            try:
                username = pwd.getpwuid(os.stat(dirname).st_uid).pw_name
            except Exception:
                username = info['username']

        print(fname)
        try:
            os.chmod(fname, 0o660)
        except Exception as e:
            print(e)

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

        # convert to HDF and store the entry in the database
        try:
            from metatlas import LcmsRun, mzml_to_hdf, store
            hdf5_file = mzml_to_hdf(fname)
            os.chmod(hdf5_file, 0o660)
            description = info['experiment'] + ' ' + info['path']
            ctime = os.stat(fname).st_ctime
            run = LcmsRun(name=info['path'], description=description,
                          username=info['username'],
                          creation_time=ctime, last_modified=ctime,
                          mzml_file=fname, hdf5_file=hdf5_file)
            store(run)
        except Exception as e:
            if 'exists but it can not be written' in str(e):
                readonly_files[username].add(dirname)
            else:
                other_errors[info['username']].append(str(e))
            print(e)

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
            body += 'To remove them all, run the following on ipython.nersc.gov:\n'
            body += 'from metatlas.metatlas_objects import find_invalid_runs, remove\n'
            body += 'remove(find_invalid_runs())\n\n'
            body += 'The invalid runs are:\n%s' % ('\n'.join(filenames))
            send_mail('Metatlas Runs are Inaccessible', username, body)
    if other_errors:
        for (username, errors) in other_errors.items():
            body = 'Errored files found while loading in Metatlas files:\n%s' % '\n'.join(errors)
            send_mail('Errors loading Metatlas files', username, body)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")

    args = parser.parse_args()
    print(args)
    update_metatlas(args.directory[0])
