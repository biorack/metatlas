from __future__ import print_function
import pymzml
import os
import pwd
import datetime
import sys
import requests
import tables

DEBUG = False


class Spectrum(tables.IsDescription):
    mz = tables.Float32Col(pos=0)
    scan_time = tables.Float32Col(pos=1)
    i = tables.Float32Col(pos=2)
    polarity = tables.UInt8Col(pos=3)
    ms_level = tables.UInt8Col(pos=4)
    precursor_MZ = tables.Float32Col(pos=5)
    precursor_intensity = tables.Float32Col(pos=6)
    collision_energy = tables.Float32Col(pos=7)


def mzml_to_hdf(in_file_name, out_file_name=None):
    """Converts in_file (mzml) to binary and stores it in out_file
    """
    if not out_file_name:
        out_file_name = in_file_name.replace('.mzML', '.h5')

    FILTERS = tables.Filters(complib='lzo', complevel=1)
    out_file = tables.open_file(out_file_name, "w", filters=FILTERS)

    table = out_file.create_table('/', 'spectra', description=Spectrum)
    if DEBUG:
        print("STATUS: Converting %s to %s (mzML to HDF)" %
              (in_file_name, out_file_name), end='')

    # Extra accessions for pymzml to read
    extraAccessions = [
        ('MS:1000016', ['value', 'unitName']),  # scan start time
        ('MS:1000129', ['name']),  # negative scan
        ('MS:1000130', ['name']),  # positive scan
        ('MS:1000744', ['name', 'value']),  # selected ion m/z
        ('MS:1000042', ['name', 'value']),  # peak intensity
        ('MS:1000045', ['name', 'value']),  # collision energy
    ]

    mzml_reader = pymzml.run.Reader(in_file_name,
                                    extraAccessions=extraAccessions)

    for spectrum in mzml_reader:
        try:
            polarity = 1 if 'positive scan' in spectrum.keys() else 0
            ms_level = spectrum['ms level']
            # check if scan start time exists. some thermo spectra
            #  are missing this value
            if 'scan start time' in spectrum.keys():
                scan_time = spectrum['scan start time'][0]
            else:
                scan_time = spectrum['MS:1000016'][0]
        except KeyError:
            continue

        # These parameters are constant for each "scan"
        precursor_MZ = 0.0
        precursor_intensity = 0.0
        collision_energy = 0.0
        if ms_level == 2:
            collision_energy = spectrum['collision energy'][1]
            if 'peak intensity' in spectrum.keys():
                precursor_intensity = spectrum['peak intensity'][1]
            precursor_MZ = spectrum['selected ion m/z'][0]

        # scan_time is also constant for each "scan"

        for mz, i in spectrum.peaks:
            # this should be changed to append numpy arrays instead of values
            # take the append out of the [mz, intensity] loop
            table.append([((float(mz), float(scan_time), float(i),
                            int(polarity), int(ms_level),
                            float(precursor_MZ), float(precursor_intensity),
                            float(collision_energy)))])
        table.flush()
        if DEBUG:
            sys.stdout.write('.')
            sys.stdout.flush()

    out_file.set_node_attr('/', "upload_date", datetime.datetime.utcnow())
    out_file.set_node_attr('/', "uploaded_by",
                           pwd.getpwuid(os.getuid())[0])

    if DEBUG:
        print('\nSaving file')
    out_file.close()
    if DEBUG:
        print("STATUS: Finished mzML to HDF conversion")

    return out_file


def get_test_data():
    dname = os.path.dirname(__file__)
    path = os.path.join(dname, 'test.mzML')

    url = ("https://drive.google.com/uc?"
           "export=download&id=0B2pT935MmTv2TlZxcTBkdGczWHM")

    if not os.path.exists(path):
        # NOTE the stream=True parameter
        print('Downloading: %s\n' % url, file=sys.stderr)
        r = requests.get(url, stream=True)
        with open(path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
        print('Download complete\n', file=sys.stderr)

    return path


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Load mzml files to HDF")
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name", required=False)
    parser.add_argument("-i", "--input", type=str,
                        help="Input mzML file", required=True)
    parser.add_argument("-d", "--debug", help="Sets debug mode",
                        action="store_true")

    args = parser.parse_args()

    # Toggles debug mode base on --debug flag
    DEBUG = args.debug

    mzml_to_hdf(args.input, args.output)
