from __future__ import print_function

import argparse
import os
import pwd
import datetime
import sys

import requests
import tables
import pymzml

DEBUG = False


class Spectrum(tables.IsDescription):
    mz = tables.Float32Col(pos=0)
    rt = tables.Float32Col(pos=1)
    i = tables.Float32Col(pos=2)
    polarity = tables.UInt8Col(pos=3)
    ms_level = tables.UInt8Col(pos=4)
    precursor_MZ = tables.Float32Col(pos=5)
    precursor_intensity = tables.Float32Col(pos=6)
    collision_energy = tables.Float32Col(pos=7)


def mzml_to_hdf(in_file_name, out_file_name=None, debug=False):
    """Converts in_file (mzml) to binary and stores it in out_file
    """
    if not out_file_name:
        out_file_name = in_file_name.replace('.mzML', '.h5')

    FILTERS = tables.Filters(complib='lzo', complevel=1)
    out_file = tables.open_file(out_file_name, "w", filters=FILTERS)

    table = out_file.create_table('/', 'spectra', description=Spectrum)

    debug = debug or DEBUG
    if debug:
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

    for (ind, spectrum) in enumerate(mzml_reader):
        try:
            polarity = 'positive scan' in spectrum or 'MS:1000130' in spectrum
            ms_level = spectrum['ms level']
            rt = spectrum.get('scan start time', spectrum['MS:1000016'])[0]
        except (KeyError, TypeError):
            continue

        precursor_MZ = 0.0
        precursor_intensity = 0.0
        collision_energy = 0.0

        if ms_level == 2:
            try:
                collision_energy = spectrum.get('collision energy',
                                                spectrum['MS:1000045'])[1]
                precursor_intensity = spectrum.get('peak intensity',
                                                   spectrum['MS:1000042'])[1]
                precursor_MZ = spectrum.get('selected ion m/z',
                                            spectrum['MS:1000744'])[0]
            except (KeyError, TypeError):
                continue

        mylist = []
        for mz, i in spectrum.peaks:
            mylist.append((mz, rt, i, polarity, ms_level,
                           precursor_MZ, precursor_intensity,
                           collision_energy))
        table.append(mylist)
        table.flush()

        if debug and not (ind % 100):
            sys.stdout.write('.')
            sys.stdout.flush()

    out_file.set_node_attr('/', "upload_date", datetime.datetime.utcnow())
    out_file.set_node_attr('/', "uploaded_by",
                           pwd.getpwuid(os.getuid())[0])

    if debug:
        print('\nSaving file')
    out_file.close()
    if debug:
        print("STATUS: Finished mzML to HDF conversion")

    return out_file


def get_test_data():
    dname = os.path.dirname(__file__)
    path = os.path.join(dname, 'test.mzML')

    # SargossoDepth/021715_QC_6_neg.mzML
    url = ("https://drive.google.com/uc?"
           "export=download&id=0B2pT935MmTv2WDMtSGhvNmxkU2M")

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


def main():
    parser = argparse.ArgumentParser(description="Load mzml files to HDF")
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name", required=False)
    parser.add_argument("-i", "--input", type=str,
                        help="Input mzML file", required=True)
    parser.add_argument("-d", "--debug", help="Sets debug mode",
                        action="store_true")

    args = parser.parse_args()

    mzml_to_hdf(args.input, args.output, args.debug)


if __name__ == '__main__':  # pragma : no cover
    main()
