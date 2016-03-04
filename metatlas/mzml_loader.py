from __future__ import print_function

import argparse
import os
import pwd
import datetime
import sys
import traceback

import tables
import pymzml

from metatlas import __version__

DEBUG = False
FORMAT_VERSION = 5
# time.time() when the latest version of the format was created
VERSION_TIMESTAMP = 1451062032


class MS1Data(tables.IsDescription):
    mz = tables.Float32Col(pos=0)
    i = tables.Float32Col(pos=1)
    rt = tables.Float32Col(pos=2)


class MS2Data(MS1Data):
    # the rest are only relevant for ms2 spectra
    precursor_MZ = tables.Float32Col(pos=5)
    precursor_intensity = tables.Float32Col(pos=6)
    collision_energy = tables.Float32Col(pos=7)


class ScanInfo(tables.IsDescription):
    rt = tables.Float32Col(pos=0)
    polarity = tables.Int16Col(pos=1)
    ms_level = tables.Int16Col(pos=2)
    max_mz = tables.Float32Col(pos=3)
    min_mz = tables.Float32Col(pos=4)
    # the rest are only relevant for ms2 spectra
    precursor_MZ = tables.Float32Col(pos=5)
    precursor_intensity = tables.Float32Col(pos=6)
    collision_energy = tables.Float32Col(pos=7)


def read_spectrum(spectrum, index):
    """Read a single spectrum

    Parameters
    ----------
    spectrum : mzML spectrum
        The mzML spectrum to parse.

    Returns
    -------
    out : list of tuples
        List of values associated with each peak in the spectrum.
    """
    polarity = 'MS:1000130' in spectrum
    ms_level = spectrum['ms level']
    rt, units = spectrum['MS:1000016']
    if units != 'minute':
        rt /= 60

    collision_energy = spectrum.get('MS:1000045', 0)
    precursor_intensity = spectrum.get('MS:1000042', 0)
    precursor_MZ = spectrum.get('MS:1000744', 0)
    min_mz = spectrum.get('lowest m/z value', 0)
    max_mz = spectrum.get('highest m/z value', 0)

    info = (rt, polarity, ms_level, min_mz, max_mz, precursor_MZ,
            precursor_intensity, collision_energy)

    if ms_level == 1:
        data = [(mz, i, rt) for (mz, i) in spectrum.peaks]
    else:
        data = [(mz, i, rt, precursor_MZ, precursor_intensity,
                 collision_energy) for (mz, i) in spectrum.peaks]

    return data, info


def mzml_to_hdf(in_file_name, out_file_name=None, debug=False):
    """Converts in_file (mzml) to binary and stores it in out_file
    """
    debug = debug or DEBUG
    if debug:
        sys.stdout.write("STATUS: Converting %s to %s (mzML to HDF)" %
              (in_file_name, out_file_name))
        sys.stdout.flush()

    # Extra accessions for pymzml to read
    extraAccessions = [
        ('MS:1000016', ['value', 'unitName']),  # scan start time
        ('MS:1000129', ['name']),  # negative scan
        ('MS:1000130', ['name']),  # positive scan
        ('MS:1000744', ['value']),  # selected ion m/z
        ('MS:1000042', ['value']),  # peak intensity
        ('MS:1000045', ['value']),  # collision energy
        ('MS:1000528', ['value']),  # lowest observed m/z
        ('MS:1000527', ['value']),  # highest observed m/z
        ('MS:1000529', ['value']),  # serial number
    ]

    try:
        mzml_reader = pymzml.run.Reader(in_file_name,
                                        extraAccessions=extraAccessions)
    except Exception as e:
        sys.stderr.write('\nMzml error: %s\n' % e)
        sys.stderr.flush()
        raise TypeError('Not a valid mzML file: "%s"' % in_file_name)

    if not out_file_name:
        out_file_name = in_file_name.replace('.mzML', '.h5')

    FILTERS = tables.Filters(complib='blosc', complevel=1)
    out_file = tables.open_file(out_file_name, "w", filters=FILTERS)
    try:
        _convert(out_file, mzml_reader, debug)
    except Exception as e:
        sys.stderr.write('\nConversion error:\n')
        traceback.print_exception(*sys.exc_info())
        sys.stderr.flush()
        sys.stdout.flush()
        raise
    finally:
        out_file.close()
    return out_file_name


def _convert(out_file, mzml_reader, debug):
    ms1_neg = out_file.create_table('/', 'ms1_neg', description=MS1Data)
    ms1_pos = out_file.create_table('/', 'ms1_pos', description=MS1Data)
    ms2_neg = out_file.create_table('/', 'ms2_neg', description=MS2Data)
    ms2_pos = out_file.create_table('/', 'ms2_pos', description=MS2Data)
    info_table = out_file.create_table('/', 'info', description=ScanInfo)

    got_first = False

    for (ind, spectrum) in enumerate(mzml_reader):

        if got_first and spectrum['id'] == 1:
            # check for a repeat
            break
        try:
            data, info = read_spectrum(spectrum, ind)
        except (KeyError, TypeError):
            continue
        except Exception as e:
            sys.stdout.write('Read spectrum error: %s\n' % e)
            sys.stdout.flush()
            continue

        if not data:
            continue

        got_first = True

        if info[2] == 1:  # ms level
            if not info[1]:  # polarity
                table = ms1_neg
            else:
                table = ms1_pos
        elif not info[1]:
            table = ms2_neg
        else:
            table = ms2_pos

        table.append(data)
        table.flush

        info_table.append([info])

        if debug and not (ind % 100):
            sys.stdout.write('.')
            sys.stdout.flush()

    info_table.flush()

    for name in ['ms1_neg', 'ms2_neg', 'ms1_pos', 'ms2_pos']:
        table = out_file.get_node('/' + name)
        table.cols.mz.create_csindex()
        table.copy(sortby='mz', newname=name + '_mz')
        table.cols.mz.remove_index()

    out_file.set_node_attr('/', "upload_date", datetime.datetime.utcnow())
    out_file.set_node_attr('/', "uploaded_by",
                           pwd.getpwuid(os.getuid())[0])

    serial = mzml_reader.param.get('MS:1000529', 'Unknown')
    out_file.set_node_attr('/', 'instrument_serial_number', serial)
    out_file.set_node_attr('/', 'format_version', FORMAT_VERSION)
    out_file.set_node_attr('/', 'metatlas_version', __version__)

    if debug:
        sys.stdout.write('\nSaving file\n')
        sys.stdout.flush()
    out_file.close()
    if debug:
        sys.stdout.write("STATUS: Finished mzML to HDF conversion\n")
        sys.stdout.flush()


def get_test_data():
    dname = os.path.dirname(__file__)

    urls = dict(basic="https://www.dropbox.com/s/j54q5amle7nyl5h/021715_QC_6_neg.mzML?dl=1",
               mix32_64="https://www.dropbox.com/s/3w83p0gjpnghqzs/QExactive_FastPolaritySwitching_Mixture_of_32_and_64_bit.mzML.xml?dl=1",
               ms_ms="https://www.dropbox.com/s/wzq7ykc436cj4ic/QExactive_Targeted_MSMS_Pos.mzML.xml?dl=1",
               wrong_fmt="https://www.dropbox.com/s/59ypkfhjgvzplm4/QExactive_Wrong_FileFormat.mzXML.xml?dl=1")

    paths = dict()
    for (name, url) in urls.items():
        path = os.path.join(dname, 'test_%s.mzML' % name)
        if name == 'wrong_fmt':
            path = os.path.join(dname, 'test.mzXML')
        if not os.path.exists(path):
            # NOTE the stream=True parameter
            print('Downloading: %s\n' % url, file=sys.stderr)
            os.system('wget %s -O %s' % (url, path))
            print('Download complete: %s bytes\n' % os.stat(path).st_size,
                  file=sys.stderr)
        else:
            print("File already exists: %s" % path)
        paths[name] = path

    return paths


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
