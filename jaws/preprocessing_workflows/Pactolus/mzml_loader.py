from __future__ import print_function

import argparse
import os
import pwd
import datetime
import sys
import traceback

import pandas as pd
import pymzml

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
    precursor_mz = spectrum.get('MS:1000744', 0)
    min_mz = spectrum.get('lowest m/z value', 0)
    max_mz = spectrum.get('highest m/z value', 0)

    if ms_level == 1:
        data = [(mz, i, rt, polarity) for (mz, i) in spectrum.peaks]
    else:
        data = [(mz, i, rt, polarity, precursor_mz, precursor_intensity,
                 collision_energy) for (mz, i) in spectrum.peaks]

    return data, ms_level, polarity


def mzml_to_df(in_file_name):
    """
    Converts in_file (mzml) to binary and stores it in a dictionary of dataframes
    """

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
    ms1_columns = ['mz','i','rt','polarity']
    ms2_columns = ['mz','i','rt','polarity','precursor_mz','precursor_intensity','collision_energy']

    try:
        mzml_reader = pymzml.run.Reader(in_file_name,
                                        extraAccessions=extraAccessions)
    except Exception as e:
        sys.stderr.write('\nMzml error: %s\n' % e)
        sys.stderr.flush()
        raise TypeError('Not a valid mzML file: "%s"' % in_file_name)

    try:
        dict_of_lists = _convert(mzml_reader)
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


def _convert(mzml_reader):
    # dict_of_dfs = {'ms1_neg',pd.DataFrame(),'ms1_neg',pd.DataFrame(),
    #                 'ms2_neg',pd.DataFrame(),'ms2_neg',pd.DataFrame()}

    got_first = False

    dict_of_lists = {'ms1_pos':[], 'ms1_neg':[], 'ms2_pos':[], 'ms2_neg':[]}

    for (ind, spectrum) in enumerate(mzml_reader):

        if got_first and spectrum['id'] == 1:
            # check for a repeat
            break
        try:
            data, ms_level, polarity = read_spectrum(spectrum, ind)
        except (KeyError, TypeError):
            continue
        except Exception as e:
            sys.stdout.write('Read spectrum error: %s\n' % e)
            sys.stdout.flush()
            continue

        if not data:
            continue

        got_first = True

        if ms_level == 1:  # ms level
            if not polarity:  # polarity
                table = 'ms1_neg'
            else:
                table = 'ms1_pos'
        elif not polarity:
            table = 'ms2_neg'
        else:
            table = 'ms2_pos'

        dict_of_lists[table].append(data)

    serial = mzml_reader.param.get('MS:1000529', 'Unknown')
    dict_of_lists['run_info'] = {'instrument_serial_number':serial}

    return dict_of_lists
