import copy
import logging
import math
import os.path
import re
import sys

from collections import defaultdict
from pathlib import Path
from textwrap import wrap

import dill
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six
import tables
from typing import List, Any
from matchms import Spectrum

from metatlas.datastructures import metatlas_objects as metob
from metatlas.io import h5_query as h5q
from metatlas.io import write_utils

logger = logging.getLogger(__name__)

MetatlasDataset = List[List[Any]]  # avoiding a circular import

def sort_atlas_table(input_atlas: str, column1: str, column2: str, istd_atlas: bool) -> str:
    """
    Reads in the atlas table, sorts it based on two columns, and writes the sorted table to a new file.

    Parameters:
    - input_atlas: Path to the input atlas file.
    - column1: The first column to sort by (ascending).
    - column2: The second column to sort by (descending).
    - istd_atlas: Whether this is an internal standard atlas with isotopic labeling.

    Returns:
    - Path to the sorted output file as a string.
    """
    # Use file extension to determine separator and output file name
    if input_atlas.lower().endswith('.csv'):
        sep = ','
        output_file = input_atlas[:-4] + '_sorted.csv'
    elif input_atlas.lower().endswith('.tsv'):
        sep = '\t'
        output_file = input_atlas[:-4] + '_sorted.tsv'
    elif input_atlas.lower().endswith('.txt'):
        sep = '\t'
        output_file = input_atlas[:-4] + '_sorted.txt'
    elif input_atlas.lower().endswith('.tab'):
        sep = '\t'
        output_file = input_atlas[:-4] + '_sorted.tab'
    else:
        raise ValueError("Input atlas file must be a .csv, .txt, .tsv, or .tab file.")

    df = pd.read_csv(input_atlas, sep=sep)

    # Check columns exist
    for col in [column1, column2]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in atlas file.")

    # Sort the DataFrame
    if istd_atlas:
        if 'label' in df.columns:
            if not df['label'].str.contains('unlabeled').any():
                logger.error("The designation 'unlabeled' does not appear in the 'label' column. Only set 'istd_atlas' to True if this is an internal standard atlas with isotopic labeling.")
                raise ValueError("'unlabeled' not found in 'label' column.")
        else:
            logger.warning("The 'label' column is missing. Not sorting with heavy isotope compound first even though 'istd_atlas' is set to True.")
    sorted_df = df.sort_values(by=[column1, column2], ascending=[True, False])

    logger.info('Writing sorted atlas to: ' + output_file)
    sorted_df.to_csv(output_file, sep=sep, index=False)

    return output_file

def create_msms_dataframe(df):
    """
    create a dataframe organized into spectra from a raw dataframe of points
    """
    #removed polarity and hdf5_file
    if 'precursor_MZ' in df.columns:
        grouped = df.groupby(['precursor_MZ','rt','precursor_intensity','collision_energy']).aggregate(lambda x: tuple(x))
    elif 'precursor_mz' in df.columns:
        grouped = df.groupby(['precursor_mz','rt','precursor_intensity','collision_energy']).aggregate(lambda x: tuple(x))
    grouped.mz = grouped.mz.apply(list)
    grouped.i = grouped.i.apply(list)
    grouped = grouped.reset_index()
    grouped.loc[:, 'spectrum'] = list(map(lambda x, y: (x, y), grouped['mz'], grouped['i']))
    grouped.loc[:, 'spectrum'] = grouped['spectrum'].apply(lambda x: list(zip(x[0], x[1])))
    grouped.drop(['mz','i'], axis=1, inplace=True)
    return grouped


def arrange_ms2_data(metatlas_dataset: MetatlasDataset, do_centroid: bool) -> pd.DataFrame:
    """
    Reformat MS2 data in metatlas dataset for efficient scoring.
    """

    file_names = get_file_names(metatlas_dataset)
    compound_names = get_compound_names(metatlas_dataset)[0]

    msms_data = []
    for file_idx, filename in enumerate(file_names):
        for compound_idx in range(len(compound_names)):

            file_compound_data = metatlas_dataset[file_idx][compound_idx]
            if 'data' not in file_compound_data['data']['msms']:
                continue

            if 'rt' not in file_compound_data['data']['msms']['data']:
                continue

            cid = file_compound_data['identification']
            cid_name = cid.compound[0].name
            adduct = cid.mz_references[0].adduct
            precursor_mz = cid.mz_references[0].mz
            mz_tolerance = cid.mz_references[0].mz_tolerance
            rt_min = cid.rt_references[0].rt_min
            rt_max = cid.rt_references[0].rt_max

            scan_rts = np.unique(file_compound_data['data']['msms']['data']['rt'])
            if scan_rts.shape[0] < 1:
                continue

            inchi_key = file_compound_data['identification'].compound[0].inchi_key

            for scan_rt in scan_rts:

                scan_mask = file_compound_data['data']['msms']['data']['rt'] == scan_rt
                mzs = file_compound_data['data']['msms']['data']['mz'][scan_mask]
                intensities = file_compound_data['data']['msms']['data']['i'][scan_mask]
                measured_precursor_mz = file_compound_data['data']['msms']['data']['precursor_MZ'][scan_mask][0]
                measured_precursor_intensity = file_compound_data['data']['msms']['data']['precursor_intensity'][scan_mask][0]

                spectrum = np.array([mzs, intensities])
                matchms_spectrum = Spectrum(spectrum[0], spectrum[1], metadata={'precursor_mz': measured_precursor_mz})

                msms_data.append({'file_name': filename, 'msms_scan': scan_rt,
                                  'measured_precursor_mz': measured_precursor_mz, 'measured_precursor_intensity': measured_precursor_intensity,
                                  'precursor_mz': precursor_mz, 'name': cid_name, 'adduct': adduct, 'inchi_key': inchi_key,
                                  'matchms_spectrum': matchms_spectrum, 'query_spectrum': spectrum, 'cid_pmz_tolerance': mz_tolerance,
                                  'cid_rt_min': rt_min, 'cid_rt_max': rt_max, 'compound_idx': compound_idx})

    return pd.DataFrame(msms_data)


def compare_EIC_to_BPC_for_file(metatlas_dataset,file_index,yscale = 'linear'):
    """
    Plot the base peak chromatogram overlaid with extracted
    ion chromatograms for features
    Input
    metatlas_dataset: a list of lists containing all atlas, file, and feature information
    file_index: the integer index of which file to plot
    """
    full_file_names = get_file_names(metatlas_dataset,full_path=True)
    base_file_names = get_file_names(metatlas_dataset,full_path=False)
    bpc = get_bpc(full_file_names[file_index])
    plt.ioff()
    fig = plt.figure()
    plt.plot(bpc.rt,bpc.i,'k-')
    for d in metatlas_dataset[file_index]:
        plt.plot(d['data']['eic']['rt'],d['data']['eic']['intensity'],'r.',alpha=0.2)
    ax = plt.gca()
    ax.set_yscale(yscale)
    ax.set_title('\n'.join(wrap(base_file_names[file_index],50)))
    ax.set_xlabel('Retention Time (min)')
    ax.set_ylabel('Intensity')
    plt.close(fig)
    return fig


def get_data_for_atlas_df_and_file(input_tuple):
    my_file, group, atlas_df, atlas = input_tuple[:4]
    extra_time = input_tuple[4] if len(input_tuple) >= 5 else 0.5
    extra_mz = input_tuple[5] if len(input_tuple) == 6 else 0.0
    if atlas.compound_identifications == []:
        return tuple([])
    df_container = remove_ms1_data_not_in_atlas(atlas_df, df_container_from_metatlas_file(my_file))
    dict_ms1_summary, dict_eic, dict_ms2 = get_data_for_atlas_and_lcmsrun(atlas_df, df_container,
                                                                          extra_time, extra_mz)
    row = []
    for i in range(atlas_df.shape[0]):
        result = {'atlas_name': atlas.name, 'atlas_unique_id': atlas.unique_id, 'lcmsrun': my_file,
                  'group': group, 'identification': copy.deepcopy(atlas.compound_identifications[i])}
        result['data'] = {'msms': {}, 'eic': dict_eic[i] if dict_eic else None}
        result['data']['ms1_summary'] = dict_ms1_summary[i] if dict_ms1_summary else None
        if dict_ms2:
            if len(dict_ms2[i]) > 0:
                result['data']['msms']['data'] = {key: np.asarray(val) for key, val in dict_ms2[i].items()}
        else:
            result['data']['msms']['data'] = []
        row.append(result)
    return tuple(row)


def get_bpc(filename,dataset='ms1_pos',integration='bpc'):
    """
    Gets the basepeak chromatogram for a file.
    filename: File can be either a metatlas lcmsrun object or a full path to an hdf5file
    dataset: ms1_pos, ms1_neg, ms2_pos, ms2_neg

    Returns:
    A pandas dataframe with the value at the maximum intensity at each retention time
    """
    df_container = df_container_from_metatlas_file(filename)
    if integration=='bpc':
        bpc = df_container[dataset].sort_values('i', ascending=False).groupby('rt', as_index=False).first().sort_values('rt',ascending=True)
        return bpc
    else:
        tic = df_container[dataset][['rt','i']].groupby('rt',as_index=False).sum().sort_values('rt',ascending=True)
        return tic

def df_container_from_metatlas_file(my_file):
    """

    """
    data_df = pd.DataFrame()

    # if my_file is a string then it's a file name - get its data
    if isinstance(my_file, six.string_types):
        filename = my_file
    else:
    # assume its a metatlas lcmsrun object
        filename = my_file.hdf5_file

    pd_h5_file = pd.HDFStore(filename, 'r')
    try:
        keys = list(pd_h5_file.keys(include='native'))
    except:
        keys = list(pd_h5_file.keys())
    pd_h5_file.close()
    df_container = {}
    for k in keys:
        if ('ms' in k) and not ('_mz' in k):
            new_df = pd.read_hdf(filename,k)
            # print new_df.keys()
            # if 'rt' in new_df.keys():
                # print 'rounding'
            # new_df.rt = new_df.rt.round()
            df_container[k[1:]] = new_df
    return df_container


def fast_nearest_interp(xi, x, y):
    """Assumes that x is monotonically increasing!!."""
    if len(y) == 1:
        return y.tolist()*len(xi)
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]


def remove_ms1_data_not_in_atlas(atlas_df, data):
    for polarity, name in [('positive', 'ms1_pos'), ('negative', 'ms1_neg')]:
        has_current_polarity = atlas_df.detected_polarity == polarity
        if any(has_current_polarity):
            atlas_mz = atlas_df[has_current_polarity].mz.copy().sort_values().values
            max_mz_tolerance = atlas_df[has_current_polarity].mz_tolerance.max()
            if data[name].shape[0] > 1:
                original_mz = data[name].mz.values
                nearest_mz = fast_nearest_interp(original_mz, atlas_mz, atlas_mz)
                data[name]['ppm_difference'] = abs(original_mz - nearest_mz) / original_mz * 1e6
                query_str = 'ppm_difference < %f' % max_mz_tolerance
                data[name] = data[name].query(query_str)
    return data


def extract(data, ids, default=None):
    """
    inputs:
        data: hierarchical data structure consisting of lists, dicts, and objects with attributes.
        ids: a list of idices, key names, and attribute names
        default: optional object
    output:
        the value stored at the location indicated by the ids list or default if the location
        does not exist.

    Strings in ids are first tried as key name and if no such key name exists, then they are
    tried as attribute names. To designate that a member of ids should be used as an attribute
    and not a key name, make it a tuple with the attribute name string as the first member, such
    as: ('attribute_name',). If you want to make it more explict to the reader, you can add a
    second member to the tuple, which will not be used, such as ('attribute_name', 'as attribute')
    """
    if data is None:
        return default
    if len(ids) == 0:
        return data
    try:
        if isinstance(ids[0], tuple):
            sub_data = getattr(data, ids[0][0])
        else:
            try:
                sub_data = data[ids[0]]
            except TypeError:
                sub_data = getattr(data, ids[0])
    except (AttributeError, IndexError, KeyError):
        return default
    else:
        return extract(sub_data, ids[1:], default)
        
def set_nested(data, ids, value):
    """
    inputs:
        data: hierarchical data structure consisting of lists, dicts, and objects with attributes.
        ids: a list of idices, key names, and attribute names
        value: object
    output:
        modifies data in place so that the value is stored at the location indicated by the ids list

    Strings in ids are first tried as key name and if no such key name exists, then they are
    tried as attribute names. To designate that a member of ids should be used as an attribute
    and not a key name, make it a tuple with the attribute name string as the first member, such
    as: ('attribute_name',). If you want to make it more explict to the reader, you can add a
    second member to the tuple, which will not be used, such as ('attribute_name', 'as attribute')
    """
    if len(ids) == 0:
        raise ValueError('ids cannot be empty')
    if len(ids) == 1:
        if isinstance(ids[0], tuple):
            setattr(data, ids[0][0], value)
        elif isinstance(ids[0], str) and hasattr(data, ids[0]):
            setattr(data, ids[0], value)
        else:
            data[ids[0]] = value  # works for list or dict
    else:
        if isinstance(ids[0], tuple):
            set_nested(getattr(data, ids[0][0]), ids[1:], value)
        elif isinstance(ids[0], str) and hasattr(data, ids[0]):
            set_nested(getattr(data, ids[0]), ids[1:], value)
        else:
            set_nested(data[ids[0]], ids[1:], value)


def make_atlas_df(atlas):
    """
    inputs:
        atlas: metatlas.datastructures.metatlas_objects.Atlas
    output:
        pandas DataFrame with one row per CompoundIdentification in atlas and each row also includes
        the first RtReference, MzReference, and Compound from the CompoundIdentification
    """
    mzs = [extract(aci, ['mz_references', 0], metob.MzReference()) for aci in atlas.compound_identifications]
    rts = [extract(aci, ['rt_references', 0], metob.RtReference()) for aci in atlas.compound_identifications]
    compounds = [extract(aci, ['compound', 0], metob.Compound()) for aci in atlas.compound_identifications]

    ci_df = metob.to_dataframe(atlas.compound_identifications)
    ci_df.rename(columns={'name': 'label'}, inplace=True)
    compound_df = metob.to_dataframe(compounds)
    compound_df.rename(columns={'name': 'compound_name', 'description': 'compound_description'}, inplace=True)
    atlas_df = pd.concat([metob.to_dataframe(rts), metob.to_dataframe(mzs), compound_df, ci_df], axis=1)

    atlas_keys = ['inchi_key', 'compound_name', 'rt_max', 'rt_min', 'rt_peak', 'rt_units',
                  'detected_polarity', 'mz', 'mz_tolerance', 'mz_tolerance_units',
                  'mono_isotopic_molecular_weight', 'pubchem_compound_id', 'synonyms', 'inchi', 'adduct',
                  'label', 'ms1_notes', 'ms2_notes', 'identification_notes']
    return atlas_df[atlas_keys]


def transfer_identification_data_to_atlas(data, atlas, ids_list=None):
    """
    inputs:
        data: metatlas_dataset object containing compound identification attribute data to transfer
        atlas: metatlas.datastructures.metatlas_objects.Atlas
    outputs:
        returns atlas with attribute data from data[0] added in
        overwrites if attribute exists in both atlas and data[0] but does not delete attribute value
        in atlas if they do not exist in data[0]
    """
    if ids_list is None:
        ids_list = [['ms1_notes'], ['ms2_notes'], ['identification_notes'], ['rt_reference', 0, 'rt_min'],
                    ['rt_reference', 0, 'rt_max'], ['rt_reference', 0, 'rt_peak']]

    out = atlas.clone(recursive=True)
    for aci, mdci in zip(out.compound_identifications, [x['identification'] for x in data[0]]):
        for ids in ids_list:
            from_data = extract(mdci, ids)
            if from_data is None:
                continue
            set_nested(aci, ids, from_data)
    return out


def get_data_for_mzrt(row, data_df_pos, data_df_neg, extra_time=0.5, use_mz='mz', extra_mz=0.0):
    mz_min = '%s >= %5.4f' % (use_mz, row.mz - row.mz*row.mz_tolerance / 1e6 - extra_mz)
    rt_min = 'rt >= %5.4f' % (row.rt_min - extra_time)
    rt_max = 'rt <= %5.4f' % (row.rt_max + extra_time)
    mz_max = '%s <= %5.4f' % (use_mz, row.mz + row.mz*row.mz_tolerance / 1e6 + extra_mz)
    ms1_query_str = f"({mz_min} & {rt_min} & {rt_max} & {mz_max})"
    data_df = data_df_pos if row.detected_polarity == 'positive' else data_df_neg
    if len(data_df) == 0:
        return pd.Series(dtype=np.float64)
    all_df = data_df.query(ms1_query_str)
    return pd.Series({'padded_feature_data': all_df,
                      'in_feature': (all_df.rt >= row.rt_min) & (all_df.rt <= row.rt_max)})


def get_ms1_summary(row):
    # A DataFrame of all points typically padded by "extra time"
    all_df = row.padded_feature_data
    # slice out ms1 data that is NOT padded by extra_time
    ms1_df = all_df[(row.in_feature)]
    num_ms1_datapoints = ms1_df.shape[0]
    has_data = num_ms1_datapoints > 0
    if has_data:
        ms1_peak_df = ms1_df.loc[ms1_df['i'].idxmax()]
        peak_area = sum(ms1_df.i)
    return pd.Series({
        'num_ms1_datapoints': num_ms1_datapoints,
        'mz_peak': ms1_peak_df.mz if has_data else np.nan,
        'rt_peak': ms1_peak_df.rt if has_data else np.nan,
        'mz_centroid': sum(ms1_df.mz * ms1_df.i) / peak_area if has_data else np.nan,
        'rt_centroid': sum(ms1_df.rt * ms1_df.i) / peak_area if has_data else np.nan,
        'peak_height': ms1_peak_df.i if has_data else np.nan,
        'peak_area': peak_area if has_data else np.nan
    })


def get_ms2_data(row):
    #A DataFrame of all points typically padded by "extra time"
    all_df = row.padded_feature_data

    #slice out ms2 data that is NOT padded by extra_time
    #ms2_df = all_df[(row.in_feature == True)]#[['collision_energy','i','mz','polarity','precursor_MZ','precursor_intensity','rt']]

    #Allow for extra_time for ms2 data
    ms2_df = all_df

    num_ms2_datapoints = ms2_df.shape[0]

    return_df = pd.Series({'ms2_datapoints':ms2_df.T,
                            'num_ms2_datapoints':num_ms2_datapoints})
    return return_df


def prefilter_ms1_dataframe_with_boundaries(data_df, rt_max, rt_min, mz_min, mz_max, extra_time=0.5, extra_mz=0.01):
    if (data_df.shape[0] == 0) | (math.isnan(rt_max)):
        return []
    return data_df.query(f"rt <= {rt_max+extra_time:5.4f} & rt >= {rt_min-extra_time:5.4f} "
                         f"& mz >= {mz_min-extra_mz:5.4f} & mz <= {mz_max+extra_mz:5.4f}")


def get_ms1_eic(row):
    #A DataFrame of all points typically padded by "extra time"
    all_df = row.padded_feature_data
    ms1_df = all_df[['i','mz','rt']]
    ms1_df = ms1_df.sort_values('rt',ascending=True)
    if ms1_df.shape[1] == 0:
        ms1_df = pd.DataFrame({'i','mz','rt'})
    return pd.Series({'eic':ms1_df.T})


def retrieve_most_intense_msms_scan(data):
    urt,idx = np.unique(data['rt'],return_index=True)
    sx = np.argsort(data['precursor_intensity'][idx])[::-1]
    prt = data['rt'][idx[sx]]
    pmz = data['precursor_MZ'][idx[sx]]
    pintensity = data['precursor_intensity'][idx[sx]]
    #setup data format for searching
    msms_data = {}
    msms_data['spectra'] = []
    msms_data['precursor_mz'] = []
    msms_data['precursor_intensity'] = []
    idx = np.argwhere((data['precursor_MZ'] == pmz[0]) & (data['rt'] == prt[0] )).flatten()
    arr = np.array([data['mz'][idx], data['i'][idx]]).T
    msms_data['spectra'] = arr
    msms_data['precursor_mz'] = pmz
    msms_data['precursor_intensity'] = pintensity
    return msms_data


def get_data_for_atlas_and_lcmsrun(atlas_df, df_container, extra_time, extra_mz):
    '''
    Accepts
    an atlas dataframe made by make_atlas_df
    a metatlas lcms file dataframe made by df_container_from_metatlas_file

    Returns python dictionaries of ms1, eic, and ms2 results for each compound in the atlas dataframe.
    '''
    # filtered the ms2 and ms1 pos and neg frames in the container by rt and mz extreme points.
    filtered = {}
    for level in ['ms1', 'ms2']:
        for polarity in ['positive', 'negative']:
            mode = f'{level}_{polarity[:3]}'
            pol = atlas_df[atlas_df.detected_polarity == polarity]
            params = [pol.rt_max.max(), pol.rt_min.min(), 0, pol.mz.max()+1, extra_time, extra_mz]
            filtered[mode] = prefilter_ms1_dataframe_with_boundaries(df_container[mode], *params)

    def get_feature_data(atlas_df, pos_df, neg_df, use_mz='mz'):
        return atlas_df.apply(
            lambda x: get_data_for_mzrt(x, pos_df, neg_df, extra_time, use_mz, extra_mz), axis=1
        )
    ms1_features = get_feature_data(atlas_df, filtered['ms1_pos'], filtered['ms1_neg'])
    if ms1_features.shape[1] == 0:
        return None, None, None
    ms2_features = get_feature_data(atlas_df, filtered['ms2_pos'], filtered['ms2_neg'], use_mz='precursor_MZ')
    return get_ms1_summary_data(ms1_features), get_eic_data(ms1_features), get_ms2_dict(ms2_features)


def get_ms2_dict(ms2_feature_data_df):
    """ extract a dict of ms2 data from the ms2 dataframe """
    ms2_data = ms2_feature_data_df.apply(get_ms2_data, axis=1)
    return [row.ms2_datapoints.T.to_dict(orient='list') if 'ms2_datapoints' in list(row.keys()) else []
            for _, row in ms2_data.iterrows()]


def get_ms1_summary_data(ms1_feature_data_df):
    """ extract a list of ms1 data from the ms1 dataframe """
    ms1_summary = ms1_feature_data_df.apply(get_ms1_summary, axis=1)
    return [dict(row) for _, row in ms1_summary.iterrows()]


def get_eic_data(ms1_feature_data_df):
    """ extract a list of eic data from the ms1 dataframe """
    ms1_eic = ms1_feature_data_df.apply(get_ms1_eic, axis=1)
    dict_eic = [row.eic.T.to_dict(orient='list') for _, row in ms1_eic.iterrows()]
    for _, value in enumerate(dict_eic):
        value['intensity'] = value.pop('i')  # rename the "i" to "intensity"
    return dict_eic


def get_unique_scan_data(data):
    """
    Input:
    data - numpy nd array containing MSMS data

    Output:
    rt - retention time of scan
    pmz - precursor m/z of scan
    Both are sorted by descending precursor ion intensity

    for data returned from h5query.get_data(),
    return the retention time and precursor m/z
    sorted by descending precursor ion intensity
    """
    urt,idx = np.unique(data['rt'],return_index=True)
    sx = np.argsort(data['precursor_intensity'][idx])[::-1]
    prt = data['rt'][idx[sx]]
    pmz = data['precursor_MZ'][idx[sx]]
    pintensity = data['precursor_intensity'][idx[sx]]
    return prt,pmz,pintensity


def get_non_redundant_precursor_list(prt,pmz,rt_cutoff,mz_cutoff):
    """
    Input:
    rt - retention time of scan
    pmz - precursor m/z of scan
    Both are sorted by descending precursor ion intensity
    rt_cutoff -
    mz_cutoff -

    Output:
    list_of_prt - list of
    list_of_pmz - list of
    """

    list_of_pmz = [] #contains list of precursor m/z [pmz1,pmz2,...,pmz_n]
    list_of_prt = [] #contains list of precursor rt [prt1,prt2,...,prt_n]

    for i in range(len(prt)):
        if len(list_of_pmz) == 0:
            # none are in the list yet; so there is nothing to check
            list_of_pmz.append(pmz[i])
            list_of_prt.append(prt[i])
        else:
            # check if new rt qualifies for inclusion
            if min(abs(list_of_prt - prt[i])) > rt_cutoff or min(abs(list_of_pmz - pmz[i])) > mz_cutoff:
                list_of_pmz.append(pmz[i])
                list_of_prt.append(prt[i])
    return list_of_prt,list_of_pmz


def organize_msms_scan_data(data,list_of_prt,list_of_pmz,list_of_pintensity):
    #setup data format for searching
    msms_data = {}
    msms_data['spectra'] = []
    msms_data['precursor_mz'] = []
    msms_data['precursor_rt'] = []
    msms_data['precursor_intensity'] = []
    for i,(prt,pmz,pintensity) in enumerate(zip(list_of_prt,list_of_pmz,list_of_pintensity)):
        idx = np.argwhere((data['precursor_MZ'] == pmz) & (data['rt'] == prt )).flatten()
        arr = np.array([data['mz'][idx], data['i'][idx]]).T
        msms_data['spectra'].append(arr)
        msms_data['precursor_mz'].append(pmz)
        msms_data['precursor_rt'].append(prt)
        msms_data['precursor_intensity'].append(pintensity)
    return msms_data


def get_data_for_a_compound(mz_ref,rt_ref,what_to_get,h5file,extra_time):
    """
    A helper function to query the various metatlas data selection
    commands for a compound defined in an experimental atlas.

    Parameters
    ----------
    MZref : a MetAtlas Object for a m/z reference Class
        this contains the m/z, m/z tolerance, and tolerance units to slice the m/z dimension
    RTref : a MetAtlas Object for a retention time reference Class
        this contains the rt min, max, peak, and units to slice the retention time dimension
    what_to_get : a list of strings
        this contains one or more of [ 'ms1_summary', 'eic', '2dhist', 'msms' ]
    h5_file : str
        Path to input_file
    polarity : int
        [0 or 1] for negative or positive ionzation

    Returns
    -------
    """
    #TODO : polarity should be handled in the experiment and not a loose parameter

    #get a pointer to the hdf5 file
    fid = tables.open_file(h5file) #TODO: should be a "with open:"

    if mz_ref.detected_polarity  == 'positive':
        polarity = 1
    else:
        polarity = 0


    mz_theor = mz_ref.mz
    if mz_ref.mz_tolerance_units  == 'ppm': #convert to ppm
        ppm_uncertainty = mz_ref.mz_tolerance
    else:
        ppm_uncertainty = mz_ref.mz_tolerance / mz_ref.mz * 1e6

#     if 'min' in rt_ref.rt_units: #convert to seconds
#     rt_min = rt_ref.rt_min / 60
#     rt_max = rt_ref.rt_max / 60
#     else:
    rt_min = rt_ref.rt_min
    rt_max = rt_ref.rt_max

    mz_min = mz_theor - mz_theor * ppm_uncertainty / 1e6
    mz_max = mz_theor + mz_theor * ppm_uncertainty / 1e6

    return_data = {}

    if 'ms1_summary' in what_to_get:
        #Get Summary Data

        #First get MS1 Raw Data
        ms_level=1
        return_data['ms1_summary'] = {}
        try:
            ms1_data = h5q.get_data(fid,
                                     ms_level=1,
                                     polarity=polarity,
                                     min_mz=mz_min,
                                     max_mz=mz_max,
                                     min_rt=rt_min,
                                     max_rt=rt_max)

            return_data['ms1_summary']['polarity'] = polarity
            return_data['ms1_summary']['mz_centroid'] = np.sum(np.multiply(ms1_data['i'],ms1_data['mz'])) / np.sum(ms1_data['i'])
            return_data['ms1_summary']['rt_centroid'] = np.sum(np.multiply(ms1_data['i'],ms1_data['rt'])) / np.sum(ms1_data['i'])
            idx = np.argmax(ms1_data['i'])
            return_data['ms1_summary']['mz_peak'] = ms1_data['mz'][idx]
            return_data['ms1_summary']['rt_peak'] = ms1_data['rt'][idx]
            return_data['ms1_summary']['peak_height'] = ms1_data['i'][idx]
            return_data['ms1_summary']['peak_area'] = np.sum(ms1_data['i'])

        except:
            return_data['ms1_summary']['polarity'] = []
            return_data['ms1_summary']['mz_centroid'] = []
            return_data['ms1_summary']['rt_centroid'] = []
            return_data['ms1_summary']['mz_peak'] = []
            return_data['ms1_summary']['rt_peak'] = []
            return_data['ms1_summary']['peak_height'] = []
            return_data['ms1_summary']['peak_area'] = []


    if 'eic' in what_to_get:
        #Get Extracted Ion Chromatogram
        # TODO : If a person calls for summary, then they will already have the MS1 raw data
        return_data['eic'] = {}
        try:
            rt,intensity = h5q.get_chromatogram(fid, mz_min, mz_max, ms_level=ms_level, polarity=polarity, min_rt = rt_min - extra_time, max_rt = rt_max + extra_time)
            return_data['eic']['rt'] = rt
            return_data['eic']['intensity'] = intensity
            return_data['eic']['polarity'] = polarity

        except:
            return_data['eic']['rt'] = []
            return_data['eic']['intensity'] = []
            return_data['eic']['polarity'] = []

    if '2dhist' in what_to_get:
        #Get 2D histogram of intensity values in m/z and retention time
        mzEdges = np.logspace(np.log10(100),np.log10(1000),10000)
#         mzEdges = np.linspace(mz_theor - 3, mz_theor + 30,100) #TODO : number of mz bins should be an optional parameter
        rtEdges = np.linspace(rt_min,rt_max,100) #TODO : number of rt bins should be an optional parameter. When not provided, it shoulddefauly to unique bins
        ms_level = 1 #TODO : ms_level should be a parameter
        return_data['2dhist'] = {}
        return_data['2dhist'] = h5q.get_heatmap(fid,mzEdges,rtEdges,ms_level,polarity)
        return_data['2dhist']['polarity'] = polarity

    if 'msms' in what_to_get:
        #Get Fragmentation Data
        ms_level=2
        return_data['msms'] = {}
        try:
            fragmentation_data = h5q.get_data(fid,
                                     ms_level=ms_level,
                                     polarity=polarity,
                                     min_mz=0,
                                     max_mz=mz_theor+2,#TODO : this needs to be a parameter
                                     min_rt=rt_min,
                                     max_rt=rt_max,
                                     min_precursor_MZ=mz_min - 0.015,
                                     max_precursor_MZ=mz_max + 0.015) #Add the 0.01 because Thermo doesn't store accurate precursor m/z
        #                     min_precursor_intensity=0, #TODO : this needs to be a parameter
        #                     max_precursor_intensity=0,#TODO : this needs to be a parameter
        #                     min_collision_energy=0,#TODO : this needs to be a parameter
        #                     max_collision_energy=0)#TODO : this needs to be a parameter
    #         prt,pmz = get_unique_scan_data(fragmentation_data)
    #         rt_cutoff = 0.23
    #         mz_cutoff = 0.05
    #         list_of_prt,list_of_pmz = get_non_redundant_precursor_list(prt,pmz,rt_cutoff,mz_cutoff)
    #         return_data['msms']['data'] = organize_msms_scan_data(fragmentation_data,list_of_prt,list_or_pmz)
            return_data['msms']['most_intense_precursor'] = retrieve_most_intense_msms_scan(fragmentation_data)
            return_data['msms']['data'] = fragmentation_data
            return_data['msms']['polarity'] = polarity
        except:
            return_data['msms']['most_intense_precursor'] = []
            return_data['msms']['data'] = []
            return_data['msms']['polarity'] = []

    fid.close() #close the file
    return return_data


def get_dill_data(fname):
    """
    Parameters
    ----------
    fname: dill file name

    Returns a list containing the data present in the dill file
    -------
    """
    if os.path.exists(fname):
        with open(fname, 'r') as handle:
            try:
                return dill.load(handle)
            except IOError as err:
                print(("I/O error({0}): {1}".format(err.errno, err.strerror)))
            except:  # handle other exceptions such as attribute errors
                print(("Unexpected error:", sys.exc_info()[0]))
    return list()


def get_group_names(data):
    """
    Parameters
    ----------
    data: either a file name (str) or a list generated from loading the dill file

    Returns list containing the group names present in the dill file
    -------
    """

    # if data is a string then it's a file name - get its data
    if isinstance(data, six.string_types):
        data = get_dill_data(data)

    group_names = list()
    for i,d in enumerate(data):
        group_names.append(d[0]['group'].name)

    return group_names

def get_group_shortnames(data):
    """
    Parameters
    ----------
    data: either a file name (str) or a list generated from loading the dill file

    Returns list containing the group short names present in the dill file
    -------
    """

    # if data is a string then it's a file name - get its data
    if isinstance(data, six.string_types):
        data = get_dill_data(data)

    group_shortnames = list()
    for i,d in enumerate(data):
        group_shortnames.append(d[0]['group'].short_name)

    return group_shortnames


def get_file_names(data,full_path=False):
    """
    Parameters
    ----------
    data: either a file name (str) or a list generated from loading the dill file
    full_path: True/False returns full path to hdf5 file or just base filename
    Returns list containing the hdf file names present in the dill file
    -------
    """

    # if data is a string then it's a file name - get its data
    if isinstance(data, six.string_types):
        data = get_dill_data(data)

    file_names = list()
    for i,d in enumerate(data):
        if full_path:
            file_names.append(d[0]['lcmsrun'].hdf5_file)
        else:
            file_names.append(os.path.basename(d[0]['lcmsrun'].hdf5_file))

    return file_names


def get_compound_names(data,use_labels=False):
    """
    Parameters
    ----------
    data: either a file name (str) or a list generated from loading the dill file

    Returns a tuple of lists containing the compound names and compound objects present in the dill file
    -------
    """
    # if data is a string then it's a file name - get its data
    if isinstance(data, six.string_types):
        data = get_dill_data(data)
    compound_names = list()
    compound_objects = list()
    if len(data) == 0:
        return (compound_names, compound_objects)
    for i,d in enumerate(data[0]):
        compound_objects.append(d['identification'])
        if use_labels:
            _str = d['identification'].name
        else:
            if len(d['identification'].compound) > 0:
                _str = d['identification'].compound[0].name
            else:
                _str = d['identification'].name
        _str = _str.split('///')[0]
        newstr = '%s_%s_%s_%s_%.4f_%.2f'%(str(i).zfill(4),_str,d['identification'].mz_references[0].detected_polarity,
                d['identification'].mz_references[0].adduct,d['identification'].mz_references[0].mz,
                d['identification'].rt_references[0].rt_peak)
        newstr = re.sub(r'\.', 'p', newstr)  # 2 or more in regexp
        newstr = re.sub(r'[\[\]]', '', newstr)
        newstr = re.sub('[^A-Za-z0-9+-]+', '_', newstr)
        newstr = re.sub('i_[A-Za-z]+_i_', '', newstr)
        if newstr[0] in ['_', '-']:
            newstr = newstr[1:]
        if newstr[-1] == '_':
            newstr = newstr[:-1]
        newstr = re.sub('[^A-Za-z0-9]{2,}', '', newstr) #2 or more in regexp
        compound_names.append(newstr)
    # If duplicate compound names exist, then append them with a number
    D = defaultdict(list)
    for i,item in enumerate(compound_names):
        D[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1}
    for k in D.keys():
        for i,f in enumerate(D[k]):
            compound_names[f] = '%s%d'%(compound_names[f],i)
    return (compound_names, compound_objects)


def make_data_sources_tables(groups, myatlas, output_loc: Path, polarity=None, overwrite=True):
    """
    polarity must be one of None, 'POS', 'NEG' or will throw ValueError
    """
    if polarity not in [None, 'POS', 'NEG']:
        raise ValueError("Polarity parameter must be one of None, 'POS', or 'NEG'.")
    prefix = f"{polarity}_" if polarity else ""
    output_dir = output_loc / f"{prefix}data_sources"
    atlas_path = output_dir / f"{prefix}atlas_metadata.tab"
    write_utils.export_dataframe(metob.to_dataframe([myatlas]), atlas_path, "atlas metadata",
                                 overwrite, sep='\t', float_format="%.8e")
    groups_path = output_dir / f"{prefix}groups_metadata.tab"
    write_utils.export_dataframe(metob.to_dataframe(groups), groups_path, "groups metadata",
                                 overwrite, sep='\t')

    atlas_df = make_atlas_df(myatlas)
    atlas_df['label'] = [cid.name for cid in myatlas.compound_identifications]
    atlas_df_path = output_dir / f"{myatlas.name}_originalatlas.tab"
    write_utils.export_dataframe(atlas_df, atlas_df_path, "atlas dataframe", overwrite, sep='\t', float_format="%.6e")

    group_path_df = pd.DataFrame(columns=['group_name', 'group_path', 'file_name'])
    loc_counter = 0
    for group in groups:
        for run in group.items:
            group_path_df.loc[loc_counter, 'group_name'] = group.name
            group_path_df.loc[loc_counter, 'group_path'] = os.path.dirname(run.mzml_file)
            group_path_df.loc[loc_counter, 'file_name'] = run.mzml_file
            loc_counter += 1

    group_path_path = output_dir / f"{prefix}groups.tab"
    write_utils.export_dataframe(group_path_df, group_path_path, "group-file mapping",
                                 overwrite, sep='\t', index=False)
