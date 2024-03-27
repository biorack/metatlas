from __future__ import print_function
from __future__ import absolute_import

import sys
import traceback
import pandas as pd
import pymzml


import argparse
import os
import pwd
import datetime
import traceback

import tables

import warnings
warnings.filterwarnings("ignore", module="plotly")
# pymzml has plotly in it and if you don't have a setup in your
# home directory, you get a warning.
from pymzml.run import Reader

from metatlas import __version__

DEBUG = False
FORMAT_VERSION = 5
# time.time() when the latest version of the format was created
VERSION_TIMESTAMP = 1451062032


class MS1Data(tables.IsDescription):
    mz = tables.Float32Col(pos=0)
    i = tables.Float32Col(pos=1)
    rt = tables.Float32Col(pos=2)
    polarity = tables.Int16Col(pos=3)


class MS2Data(MS1Data):
    # the rest are only relevant for ms2 spectra
    precursor_MZ = tables.Float32Col(pos=4)
    precursor_intensity = tables.Float32Col(pos=5)
    collision_energy = tables.Float32Col(pos=6)


# class ScanInfo(tables.IsDescription):
#     rt = tables.Float32Col(pos=0)
#     polarity = tables.Int16Col(pos=1)
#     ms_level = tables.Int16Col(pos=2)
#     max_mz = tables.Float32Col(pos=3)
#     min_mz = tables.Float32Col(pos=4)
#     # the rest are only relevant for ms2 spectra
#     precursor_MZ = tables.Float32Col(pos=5)
#     precursor_intensity = tables.Float32Col(pos=6)
#     collision_energy = tables.Float32Col(pos=7)

def read_spectrum(spectrum):
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
    if spectrum['negative scan'] is True:
        polarity = 0
    elif spectrum['positive scan'] is True:
        polarity = 1
    else:
        raise('EXITING: polarity is unknown')


    ms_level = spectrum.ms_level
    rt = spectrum.scan_time_in_minutes()


#     min_mz = spectrum.get('lowest m/z value', 0)
#     max_mz = spectrum.get('highest m/z value', 0)

    if ms_level == 1:
        info = [rt, polarity, ms_level]
        data = [(mz, i, rt, polarity) for (mz, i) in spectrum.peaks('centroided')]
    else:
        info = [rt, polarity, ms_level]
        prec = spectrum.selected_precursors
        if len(prec) != 1:
            raise('EXITING: Not only one precursor ion selected')
        else:
            prec = prec[0]

        collision_energy = spectrum['collision energy']
        precursor_intensity = prec['i']
        precursor_mz = prec['mz']
        data = [(mz, i, rt, polarity, precursor_mz, precursor_intensity,
                 collision_energy) for (mz, i) in spectrum.peaks('centroided')]
    return data, info


def mzml_to_df(in_file_name):
    """
    Converts in_file (mzml) to binary and stores it in a dictionary of dataframes
    """

#     # Extra accessions for pymzml to read
#     extraAccessions = [
#         ('MS:1000016', ['value', 'unitName']),  # scan start time
#         ('MS:1000129', ['name']),  # negative scan
#         ('MS:1000130', ['name']),  # positive scan
#         ('MS:1000744', ['value']),  # selected ion m/z
#         ('MS:1000042', ['value']),  # peak intensity
#         ('MS:1000045', ['value']),  # collision energy
#         ('MS:1000528', ['value']),  # lowest observed m/z
#         ('MS:1000527', ['value']),  # highest observed m/z
#         ('MS:1000529', ['value']),  # serial number
#     ]
    ms1_columns = ['mz','i','rt','polarity']
    ms2_columns = ['mz','i','rt','polarity','precursor_mz','precursor_intensity','collision_energy']

    try:
        mzml_reader = pymzml.run.Reader(in_file_name)
#                                         extraAccessions=extraAccessions)
    except Exception as e:
        sys.stderr.write('\nMzml error: %s\n' % e)
        sys.stderr.flush()
        raise TypeError('Not a valid mzML file: "%s"' % in_file_name)

    try:
        dict_of_lists = _convert(mzml_reader,out_file=None)
        for name in ['ms1_neg', 'ms2_neg', 'ms1_pos', 'ms2_pos']:
            if 'ms1' in name:
                if dict_of_lists[name]:
                    #flatten list of lists
                    dict_of_lists[name] = [item for sublist in dict_of_lists[name] for item in sublist]
                    dict_of_lists[name] = pd.DataFrame(dict_of_lists[name],columns=ms1_columns)
                    dict_of_lists[name].loc[~dict_of_lists[name]['polarity']==True,'polarity'] = 'negative'
                    dict_of_lists[name].loc[dict_of_lists[name]['polarity']==True,'polarity'] = 'positive'
            else:
                if dict_of_lists[name]:
                    #flatten list of lists
                    dict_of_lists[name] = [item for sublist in dict_of_lists[name] for item in sublist]
                    dict_of_lists[name] = pd.DataFrame(dict_of_lists[name],columns=ms2_columns)
                    dict_of_lists[name].loc[~dict_of_lists[name]['polarity']==True,'polarity'] = 'negative'
                    dict_of_lists[name].loc[dict_of_lists[name]['polarity']==True,'polarity'] = 'positive'

    except Exception as e:
        sys.stderr.write('\nConversion error:\n')
        traceback.print_exception(*sys.exc_info())
        sys.stderr.flush()
        sys.stdout.flush()
        raise
    return dict_of_lists #returns a dictionary of dataframes


def mzml_to_hdf(in_file_name, out_file=None, debug=False, remove_easyic_signal=True):
    """Converts in_file (mzml) to binary and stores it in out_file"""
    debug = debug or DEBUG
    if debug:
        sys.stdout.write("STATUS: Converting %s to %s (mzML to HDF)" %
              (in_file_name, out_file))
        sys.stdout.flush()
#     # Extra accessions for pymzml to read
#     extraAccessions = [
#         ('MS:1000016', ['value', 'unitName']),  # scan start time
#         ('MS:1000129', ['name']),  # negative scan
#         ('MS:1000130', ['name']),  # positive scan
#         ('MS:1000744', ['value']),  # selected ion m/z
#         ('MS:1000042', ['value']),  # peak intensity
#         ('MS:1000045', ['value']),  # collision energy
#         ('MS:1000528', ['value']),  # lowest observed m/z
#         ('MS:1000527', ['value']),  # highest observed m/z
#         ('MS:1000529', ['value']),  # serial number
#     ]

    try:
        mzml_reader = Reader(in_file_name)
#                                         extraAccessions=extraAccessions)
    except Exception as e:
        sys.stderr.write('\nMzml error: %s\n' % e)
        sys.stderr.flush()
        raise TypeError('Not a valid mzML file: "%s"' % in_file_name)

    if out_file is None:
        out_file = in_file_name.replace('.mzML', '.h5')

    try:
        _convert(mzml_reader,out_file=out_file, debug=debug)
    except Exception as e:
        sys.stderr.write('\nConversion error:\n')
        traceback.print_exception(*sys.exc_info())
        sys.stderr.flush()
        sys.stdout.flush()
        raise

    return out_file

def remove_easyic_signal(spectrum_data, easyic_mzs={1:202.07770, 0:202.07880}, mz_tolerance=0.001):
    """Remove peaks from spectral data that match the internal calibrant m/z.
    
    Depending on Easy-IC method setting on modern Thermo Orbitrap Mass Spectrometers, it may be necessary
    to remove the fluoranthene lock mass signal.
    """
    
    targeted_mz = easyic_mzs[spectrum_data[0][3]]

    cleaned_spectrum = []
    for peak_data in spectrum_data:
        
        if abs(peak_data[0] - targeted_mz) > mz_tolerance:
            cleaned_spectrum.append(peak_data)
    
    return cleaned_spectrum

def _convert(mzml_reader, out_file=None, filter_easyic_signal=True, debug=None):
    """Convert and save mzML spectra to h5 file.
    
    Parameters
    ----------
    remove_easyic_signal : Bool
        Whether or not to remove lock mass signal during conversion.
    """
    if out_file is not None:
        FILTERS = tables.Filters(complib='blosc', complevel=1)
        out_file = tables.open_file(out_file, "w", filters=FILTERS)
        ms1_neg = out_file.create_table('/', 'ms1_neg', description=MS1Data)
        ms1_pos = out_file.create_table('/', 'ms1_pos', description=MS1Data)
        ms2_neg = out_file.create_table('/', 'ms2_neg', description=MS2Data)
        ms2_pos = out_file.create_table('/', 'ms2_pos', description=MS2Data)
    else:
        out_dict = {}
        out_dict['ms1_neg'] = []
        out_dict['ms1_pos'] = []
        out_dict['ms2_neg'] = []
        out_dict['ms2_pos'] = []
#         info_table = out_file.create_table('/', 'info', description=ScanInfo)

    got_first = False

    for spectrum in mzml_reader:

#         if got_first and spectrum['id'] == 1:
#             # check for a repeat
#             break
        try:
            data,info = read_spectrum(spectrum)

        except (KeyError, TypeError):
            continue
        except Exception as e:
            sys.stdout.write('Read spectrum error: %s\n' % e)
            sys.stdout.flush()
            continue

        if not data:
            continue
            
        if filter_easyic_signal:
            data = remove_easyic_signal(data)

#         got_first = True
        if out_file is not None:
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
        else:
            if info[2] == 1:  # ms level
                if not info[1]:  # polarity
                    k = 'ms1_neg'
                else:
                    k = 'ms1_pos'
            elif not info[1]:
                k = 'ms2_neg'
            else:
                k = 'ms2_pos'
            out_dict[k].append(data)


    if out_file is not None:
        for name in ['ms1_neg', 'ms2_neg', 'ms1_pos', 'ms2_pos']:
            table = out_file.get_node('/' + name)
            table.cols.mz.create_csindex()
            table.copy(sortby='mz', newname=name + '_mz')
            table.cols.mz.remove_index()#             info_table.append([info])
        out_file.set_node_attr('/', 'format_version', FORMAT_VERSION)
        out_file.set_node_attr('/', 'metatlas_version', __version__)

        if debug:
            sys.stdout.write('\nSaving file\n')
            sys.stdout.flush()

        out_file.close()
        if debug:
            sys.stdout.write("STATUS: Finished mzML to HDF conversion\n")
            sys.stdout.flush()

#         if debug and not (ind % 100):
#             sys.stdout.write('.')
#             sys.stdout.flush()

#     info_table.flush()


#     out_file.set_node_attr('/', "upload_date", datetime.datetime.utcnow())
#     out_file.set_node_attr('/', "uploaded_by",
#                            pwd.getpwuid(os.getuid())[0])

#     serial = mzml_reader.param.get('MS:1000529', 'Unknown')
#     out_file.set_node_attr('/', 'instrument_serial_number', serial)
    if out_file is None:
        return out_dict

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
            os.system('wget -nv %s -O %s' % (url, path))
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

    mzml_to_hdf(args.input, out_file=args.output, debug=args.debug)


if __name__ == '__main__':  # pragma : no cover
    main()



# def _convert(mzml_reader):
#     # dict_of_dfs = {'ms1_neg',pd.DataFrame(),'ms1_neg',pd.DataFrame(),
#     #                 'ms2_neg',pd.DataFrame(),'ms2_neg',pd.DataFrame()}

#     got_first = False

#     dict_of_lists = {'ms1_pos':[], 'ms1_neg':[], 'ms2_pos':[], 'ms2_neg':[]}

#     for (ind, spectrum) in enumerate(mzml_reader):

#         if got_first and spectrum['id'] == 1:
#             # check for a repeat
#             break
#         try:
#             data, ms_level, polarity = read_spectrum(spectrum, ind)
#         except (KeyError, TypeError):
#             continue
#         except Exception as e:
#             sys.stdout.write('Read spectrum error: %s\n' % e)
#             sys.stdout.flush()
#             continue

#         if not data:
#             continue

#         got_first = True

#         if ms_level == 1:  # ms level
#             if not polarity:  # polarity
#                 table = 'ms1_neg'
#             else:
#                 table = 'ms1_pos'
#         elif not polarity:
#             table = 'ms2_neg'
#         else:
#             table = 'ms2_pos'

#         dict_of_lists[table].append(data)

# #     serial = mzml_reader.param.get('MS:1000529', 'Unknown')
# #     dict_of_lists['run_info'] = {'instrument_serial_number':serial}

#     return dict_of_lists
