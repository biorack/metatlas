from watchdog.utils.dirsnapshot import DirectorySnapshot, DirectorySnapshotDiff
import pickle
import os
import re


TEMP_FILE_NAME = os.path.join(os.path.dirname(os.path.dirname(__file__)), "snap.shot")


class SnapshotWrapper:

    @classmethod
    def load(cls, file_name):
        print(file_name)
        f = open(file_name, "rb")
        return pickle.load(f)

    def __init__(self, directory, file_name):
        self.directory = directory
        self.file_name = file_name
        self.snapshot = DirectorySnapshot(directory)

    def save(self, to_file=None):
        if to_file is None:
            snap_file = open(self.file_name, "wb")
        else:
            snap_file = open(to_file, "wb")
        pickle.dump(self, snap_file)
        snap_file.close()
        return self

    def get_added_files(self):
        new_snap = DirectorySnapshot(self.directory)
        old_snap = self.snapshot
        return DirectorySnapshotDiff(old_snap, new_snap).files_created

    def update(self):
        self.snapshot = DirectorySnapshot(self.directory)
        return self


def update_metatlas(directory):
    if not os.path.exists(TEMP_FILE_NAME):
        snw = SnapshotWrapper(directory, TEMP_FILE_NAME)
    else:
        snw = SnapshotWrapper.load(TEMP_FILE_NAME)
    new_files = snw.get_added_files()

    patt = re.compile(r".+\/raw_data\/(?P<username>[^/]+)\/(?P<experiment>[^/]+)\/(?P<path>.+)")

    for fname in new_files:
        info = patt.match(os.path.abspath(fname))
        if info:
            info = info.groupdict()
        else:
            print("Invalid path name", fname)
            continue
        print(fname)

        # convert to HDF and store the entry in the database
        try:
            from metatlas import LcmsRun, mzml_to_hdf
            hdf_file = mzml_to_hdf(fname)
            description = info['experiment'] + ' ' + info['path']
            ctime = os.stat(fname).st_ctime
            run = LcmsRun(name=info['path'], description=description,
                          created_by=info['username'],
                          modified_by=info['username'],
                          created=ctime, last_modified=ctime,
                          mzml_file=fname, hdf_file=hdf_file)
            run.store()
        except TypeError as e:
            print(e)

    snw.update().save()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")

    args = parser.parse_args()
    print(args)
    update_metatlas(args.directory[0])
