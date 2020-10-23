from __future__ import absolute_import
import os
import fcntl
import pwd
import re
import shutil
import sys
import time
import random
import smtplib
import traceback
import time

from collections import defaultdict
from subprocess import check_output
from datetime import datetime, time as dtime

from metatlas.mzml_loader import VERSION_TIMESTAMP
from metatlas import LcmsRun, mzml_to_hdf, store, retrieve

ADMIN = 'bpb'


def send_mail(subject, username, body, force=False):
    """Send the mail only once per day."""
    now = datetime.now()
    if force or dtime(00, 00) <= now.time() <= dtime(00, 10):
        sender = 'pasteur@nersc.gov'
        receivers = ['%s@nersc.gov' % username, '%s@nersc.gov' % ADMIN]
        message = """\
From: %s
To: %s
Subject: %s

%s
        """ % (sender, ", ".join(receivers), subject, body)
        try:
            smtpObj = smtplib.SMTP('localhost')
            smtpObj.sendmail(sender, receivers, message)
            sys.stdout.write("Successfully sent email to %s\n" % username)
            sys.stdout.flush()
        except smtplib.SMTPException:
            sys.stderr.write("Error: unable to send email to %s\n" % username)
            sys.stdout.flush()

def get_acqtime_from_mzml(mzml_file):
    startTimeStamp=None
    with open(mzml_file) as mzml:
        for line in mzml:
            if 'startTimeStamp' in line:
                startTimeStamp = line.split('startTimeStamp="')[1].split('"')[0].replace('T',' ').rstrip('Z')
                break
#     print startTimeStamp
    if not '-infinity' in startTimeStamp:
        date_object = datetime.strptime(startTimeStamp, '%Y-%m-%d %H:%M:%S')
        utc_timestamp = int(time.mktime(date_object.timetuple()))
    else:
        utc_timestamp = int(0)
    return utc_timestamp

def update_metatlas(directory):
    readonly_files = defaultdict(set)
    other_errors = defaultdict(list)
    directory = os.path.abspath(directory)

    # Sleep a random amount of time to avoid running at the same time as
    # other processes.
    time.sleep(random.random() * 2)
    mzml_files = check_output('find %s -name "*.mzML"' % directory, shell=True)
    mzml_files = mzml_files.decode('utf-8').splitlines()

    # Find valid h5 files newer than the format version timestamp.
    delta = int((time.time() - VERSION_TIMESTAMP) / 60)
    check = 'find %s -name "*.h5" -mmin -%s -size +2k' % (directory, delta)
    valid_files = check_output(check, shell=True).decode('utf-8').splitlines()
    valid_files = set(valid_files)

    new_files = []
    for mzml_file in mzml_files:
        if mzml_file.replace('.mzML', '.h5') not in valid_files:
            new_files.append(mzml_file)



    patt = re.compile(r".+\/raw_data\/(?P<username>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")

    sys.stdout.write('Found %s files\n' % len(new_files))
    sys.stdout.flush()


    for (ind, fname) in enumerate(new_files):
        sys.stdout.write('(%s of %s): %s\n' % (ind + 1, len(new_files), fname))
        sys.stdout.flush()

        # Get relevant information about the file.
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

        # Change to read only.
        try:
            os.chmod(fname, 0o660)
        except Exception as e:
            sys.stderr.write(str(e) + '\n')
            sys.stderr.flush()

        # Copy the original file to a pasteur backup.
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

        # Get a lock on the mzml file to prevent interference.
        try:
            fid = open(fname, 'r')
            fcntl.flock(fid, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except IOError:
            fid.close()
            msg = '%s already converting in another process\n' % fname
            sys.stderr.write(msg)
            sys.stderr.flush()
            continue

        # Convert to HDF and store the entry in the database.
        try:
            hdf5_file = fname.replace('mzML', 'h5')
            
            #Get Acquisition Time Here
            acquisition_time = get_acqtime_from_mzml(fname)
            mzml_to_hdf(fname, hdf5_file, True)
            os.chmod(hdf5_file, 0o660)
            description = info['experiment'] + ' ' + info['path']
            ctime = os.stat(fname).st_ctime
            # Add this to the database unless it is already there
            try:
                runs = retrieve('lcmsrun', username='*', mzml_file=fname)
            except Exception:
                runs = list()
            if not len(runs):
                run = LcmsRun(name=info['path'], description=description,
                              username=info['username'],
                              experiment=info['experiment'],
                              creation_time=ctime, last_modified=ctime,
                              mzml_file=fname, hdf5_file=hdf5_file, acquisition_time = acquisition_time)
                store(run)
        except Exception as e:
            if 'exists but it can not be written' in str(e):
                readonly_files[username].add(dirname)
            else:
                msg = traceback.format_exception(*sys.exc_info())
                msg.insert(0, 'Cannot convert %s' % fname)
                other_errors[info['username']].append('\n'.join(msg))
            sys.stderr.write(str(e) + '\n')
            sys.stderr.flush()
            try:
                os.remove(hdf5_file)
            except:
                pass
        finally:
            fid.close()

    # Handle errors.
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
            send_mail('Metatlas Runs are Invalid', username, body)
    if other_errors:
        for (username, errors) in other_errors.items():
            body = 'Errored files found while loading in Metatlas files:\n\n%s' % '\n********************************\n'.join(errors)
            send_mail('Errors loading Metatlas files', username, body)
    sys.stdout.write('Done!\n')
    sys.stdout.flush()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")
    args = parser.parse_args()
    sys.stdout.write(str(args) + '\n')
    sys.stdout.flush()
    update_metatlas(args.directory[0])
