""" mzml to h5 file conversion """

import argparse
import functools
import logging
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

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s", level=logging.INFO)

readonly_files = {}  # username (or uid) | a set of files associated with them
other_errors = {}  # info with user | list of error messages
patt = re.compile(r".+\/raw_data\/(?P<username>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")

MEMBERS_CMD = "getent group metatlas | cut -d: -f4"
ALL_USERS = set(sorted(subprocess.check_output(MEMBERS_CMD, shell=True, text=True).strip().split(",")))
REMOVE_USERS = {
    "msdata",
    "jaws",
    "jgi_dna",
    "vrsingan",
    "wjholtz",
    "mjblow",
    "greensi",
    "annau",
    "jfroula",
    "pasteur",
}
USERS = tuple(ALL_USERS - REMOVE_USERS)
DEFAULT_USERNAME = "smkosina"
EXPLICIT_USERNAMES = {"ag": "agolini", "ao": "arosborn"}


@functools.lru_cache
def _initials_to_username(initials):
    if initials == "":
        return None
    if initials in EXPLICIT_USERNAMES:
        return EXPLICIT_USERNAMES[initials]
    for user in USERS:
        if user.startswith(initials):
            return user
    for user in USERS:
        pat = re.compile(f"^{initials[0]}[a-z]{initials[1]}")
        if pat.match(user):
            return user
    return None


def move_file(src, dest):
    """move file and create directories if needed"""
    assert os.path.isfile(src)
    dest_dir = dest if os.path.isdir(dest) else os.path.dirname(dest)
    os.makedirs(dest_dir, exist_ok=True)
    shutil.move(src, dest)


def _file_name_to_username(file_name, default):
    """extract initials from filename and convert to nersc username"""
    initials_field = os.path.basename(file_name).split("_")[1].lower()
    for initials in initials_field.split("-"):  # from left to right
        username = _initials_to_username(initials.replace("_", ""))
        if username is not None:
            return username
    return default


def get_acqtime_from_mzml(mzml_file):
    start_time = None
    with open(mzml_file, "r", encoding="utf-8") as mzml:
        for line in mzml:
            if "start_time" in line:
                start_time = line.split('start_time="')[1].split('"')[0].replace("T", " ").rstrip("Z")
                break
    if start_time is not None and "-infinity" not in start_time:
        date_object = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
        utc_timestamp = int(time.mktime(date_object.timetuple()))
    else:
        utc_timestamp = int(0)
    return utc_timestamp


def mzml_to_h5_and_add_to_db(mzml_file_name: str) -> None:
    """converts a single file and inserts a record in lcmsruns table"""
    logger.info("Converting mzML file %s", mzml_file_name)

    pat = re.compile(r".+\/raw_data\/(?P<sub_dir>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")
    mzml_file_name = os.path.abspath(mzml_file_name)
    try:
        info = pat.match(mzml_file_name).groupdict()
    except AttributeError:
        logger.error("Invalid path name: %s", mzml_file_name)
        return
    try:
        hdf5_file = mzml_file_name.replace("mzML", "h5")
        logger.info("Generating h5 file: %s", hdf5_file)
        mzml_to_hdf(mzml_file_name, hdf5_file, True)
        try:
            runs = retrieve("lcmsrun", username="*", mzml_file=mzml_file_name)
        except Exception:
            runs = []
        if not runs:
            username = _file_name_to_username(mzml_file_name, DEFAULT_USERNAME)
            ctime = os.stat(mzml_file_name).st_ctime
            logger.info("LCMS run not in DB, inserting new entry.")
            run = LcmsRun(
                name=info["path"],
                description=f"{info['experiment']} {info['path']}",
                username=username,
                experiment=info["experiment"],
                creation_time=ctime,
                last_modified=ctime,
                mzml_file=mzml_file_name,
                hdf5_file=hdf5_file,
                acquisition_time=get_acqtime_from_mzml(mzml_file_name),
            )
            store(run)
    except Exception as e:
        logger.error("During file conversion: %s", str(e))
        if "exists but it can not be written" in str(e):
            dirname = os.path.dirname(mzml_file_name)
            logger.error("Cannot write to file within directory %s", dirname)
        else:
            fail_path = mzml_file_name.replace("raw_data", "conversion_failures")
            logger.error("Moving mzml file to %s", fail_path)
            move_file(mzml_file_name, fail_path)
        try:
            os.remove(hdf5_file)
        except:
            pass


def convert(ind, fname):
    """Helper function, converts a single file"""
    logger.info("Converting file number %d: %s", ind + 1, fname)

    # Get relevant information about the file.
    username = _file_name_to_username(fname, DEFAULT_USERNAME)
    info = patt.match(os.path.abspath(fname))
    if info:
        info = info.groupdict()
    else:
        logger.error("Invalid path name: %s", fname)
        return
    dirname = os.path.dirname(fname)

    # Convert to HDF and store the entry in the database.
    try:
        hdf5_file = fname.replace("mzML", "h5")
        logger.info("Generating h5 file: %s", hdf5_file)
        mzml_to_hdf(fname, hdf5_file, True)
        os.chmod(hdf5_file, 0o660)  # this can be changed to 0o440 once everyone is on the current code
        # Add this to the database unless it is already there
        try:
            runs = retrieve("lcmsrun", username="*", mzml_file=fname)
        except Exception:
            runs = []
        if not runs:
            ctime = os.stat(fname).st_ctime
            logger.info("LCMS run not in DB, inserting new entry.")
            run = LcmsRun(
                name=info["path"],
                description=f"{info['experiment']} {info['path']}",
                username=username,
                experiment=info["experiment"],
                creation_time=ctime,
                last_modified=ctime,
                mzml_file=fname,
                hdf5_file=hdf5_file,
                acquisition_time=get_acqtime_from_mzml(fname),
            )
            store(run)
    except Exception as e:
        logger.error("During file conversion: %s", str(e))
        if "exists but it can not be written" in str(e):
            logger.error("Cannot write to file within directory %s", dirname)
            if username not in readonly_files:
                readonly_files[username] = set()
            readonly_files[username].add(dirname)
        else:
            msg = traceback.format_exception(*sys.exc_info())
            msg.insert(0, f"Cannot convert {fname}")
            dat = username
            if dat not in other_errors:
                other_errors[username] = []
            other_errors[username].append("\n".join(msg))
            fail_path = fname.replace("raw_data", "conversion_failures")
            logger.error("Moving mzml file to %s", fail_path)
            move_file(fname, fail_path)
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
    mzml_files = mzml_files.decode("utf-8").splitlines()

    # Find valid h5 files newer than the format version timestamp.
    delta = int((time.time() - VERSION_TIMESTAMP) / 60)
    check = f'find {directory} -name "*.h5" -mmin -{delta} -size +2k'
    valid_files = set(check_output(check, shell=True).decode("utf-8").splitlines())
    new_files = [file for file in mzml_files if file.replace(".mzML", ".h5") not in valid_files]

    if new_files:
        logger.info("Found %d files", len(new_files))
        for ind, ffff in enumerate(new_files):
            convert(ind, ffff)
        if readonly_files:
            for (username, dirnames) in readonly_files.items():
                logger.info("Sending email to %s about inaccessible files.", username)
                body = (
                    "Please log in to NERSC and run 'chmod g+rwXs' on the "
                    "following directories:\n%s" % ("\n".join(dirnames))
                )
                send_mail("Metatlas Files are Inaccessible", username, body)
        if other_errors:
            for (username, errors) in other_errors.items():
                logger.info("Sending email to %s about conversion error.", username)
                body = (
                    "Errored files found while loading in Metatlas files:\n\n%s"
                    % "\n********************************\n".join(errors)
                )
                send_mail("Errors loading Metatlas files", username, body)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")
    args = parser.parse_args()
    logger.info("Monitoring directory: %s", args.directory[0])
    update_metatlas(args.directory[0])
    logger.info("Done! - file_converter.py run has completed.")
