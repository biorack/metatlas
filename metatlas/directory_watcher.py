import os
import re
import sys
import shutil


def update_metatlas(directory):

    new_files = []
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
            shutil.copy(fname, pasteur_path)

        # todo: search the database at the end for unreachable files

        # convert to HDF and store the entry in the database
        try:
            from metatlas import LcmsRun, mzml_to_hdf
            hdf5_file = mzml_to_hdf(fname)
            os.chmod(hdf5_file, 0o660)
            description = info['experiment'] + ' ' + info['path']
            ctime = os.stat(fname).st_ctime
            run = LcmsRun(name=info['path'], description=description,
                          created_by=info['username'],
                          modified_by=info['username'],
                          created=ctime, last_modified=ctime,
                          mzml_file=fname, hdf5_file=hdf5_file)
            run.store()
        except Exception as e:
            print(e)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Watchdog to monitor directory for new files")
    parser.add_argument("directory", type=str, nargs=1, help="Directory to watch")

    args = parser.parse_args()
    print(args)
    update_metatlas(args.directory[0])
