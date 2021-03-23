from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os.path
import sys
import copy
import tables
from metatlas.datastructures import metatlas_objects as metob
import pandas as pd
from textwrap import wrap
import matplotlib.pyplot as plt
import six
from six.moves import map
from six.moves import range
from six.moves import zip


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
    grouped['spectrum'] = list(map(lambda x,y:(x,y),grouped['mz'],grouped['i']))
    grouped['spectrum'] = grouped['spectrum'].apply(lambda x: list(zip(x[0],x[1])))
    grouped.drop(['mz','i'], axis=1, inplace=True)
    return grouped

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
    my_file = input_tuple[0]
    my_group = input_tuple[1]
    atlas_df = input_tuple[2]
    myAtlas = input_tuple[3]
    extra_time = 0.5
    extra_mz = 0.0
    if len(input_tuple) == 6:
        extra_time = input_tuple[4]
        extra_mz = input_tuple[5]
    elif len(input_tuple) == 5:
        extra_time = input_tuple[4]
    
    df_container = df_container_from_metatlas_file(my_file)

    df_container = remove_ms1_data_not_in_atlas(atlas_df,df_container)
    dict_ms1_summary,dict_eic,dict_ms2 = get_data_for_atlas_and_lcmsrun(atlas_df,df_container,extra_time, extra_mz)
    row = []
    for i in range(atlas_df.shape[0]):
        result = {}
        result['atlas_name'] = myAtlas.name
        result['atlas_unique_id'] = myAtlas.unique_id
        result['lcmsrun'] = my_file
        result['group'] = my_group
        temp_compound = copy.deepcopy(myAtlas.compound_identifications[i])
        result['identification'] = temp_compound
        result['data'] = {}
        if dict_eic:
            result['data']['eic'] = dict_eic[i]
        else:
            result['data']['eic'] = None
        if dict_ms1_summary:
            result['data']['ms1_summary'] = dict_ms1_summary[i]
        else:
            result['data']['ms1_summary'] = None

        result['data']['msms'] = {}
        if dict_ms2:
            if len(dict_ms2[i])>0:#dict_ms2[i]['mz']:
                for k in dict_ms2[0].keys():
                    dict_ms2[i][k] = np.asarray(dict_ms2[i][k])
        #                 if temp_compound.mz_references[0].observed_polarity == 'positive':
        #                     dict_ms2[i]['polarity'] = dict_ms2[i]['mz'] * 0.0 + 1.0
        #                 else:
        #                     dict_ms2[i]['polarity'] = dict_ms2[i]['mz'] * 0.0
                result['data']['msms']['data'] = dict_ms2[i]
        else:
            result['data']['msms']['data'] = []
        row.append(result)
    return row

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
        
    pd_h5_file  = pd.HDFStore(filename)
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
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

def remove_ms1_data_not_in_atlas(atlas_df,data):
    things_to_do = [('positive','ms1_pos'),('negative','ms1_neg')]
    for thing in things_to_do:
        if sum(atlas_df.detected_polarity == thing[0])>0:
            atlas_mz = atlas_df[atlas_df.detected_polarity == thing[0]].mz.copy()
            atlas_mz = atlas_mz.sort_values()
            atlas_mz = atlas_mz.values
            max_mz_tolerance = atlas_df[atlas_df.detected_polarity == thing[0]].mz_tolerance.max()
            if data[thing[1]].shape[0]>1:
                original_mz = data[thing[1]].mz.values
                nearest_mz = fast_nearest_interp(original_mz,atlas_mz,atlas_mz)
                data[thing[1]]['ppm_difference'] = abs(original_mz - nearest_mz) / original_mz * 1e6
                query_str = 'ppm_difference < %f'%(max_mz_tolerance)
                data[thing[1]] = data[thing[1]].query(query_str)
    return data

def make_atlas_df(atlas):
    mz = []
    rt = []
    atlas_compound = []
    label = []
    for compound in atlas.compound_identifications:
        label.append(compound.name)
        if compound.mz_references:
            mz.append(compound.mz_references[0])
        else:
            mz.append(metob.MzReference())
        if compound.rt_references:
            rt.append(compound.rt_references[0])
        else:
            rt.append(metob.RtReference())
        if compound.compound:
            atlas_compound.append(compound.compound[0])
        else:
            atlas_compound.append(metob.Compound())

    compound_df = metob.to_dataframe(atlas_compound)
    compound_df.rename(columns = {'name':'compound_name','description':'compound_description'},inplace=True)
    #.rename(columns = {'name':'compound_name'}, inplace = True)
    atlas_df = pd.concat([metob.to_dataframe(rt),metob.to_dataframe(mz), compound_df],axis=1)
    # atlas_df['label'] = label

    atlas_keys = [u'inchi_key','compound_name',u'rt_max', u'rt_min', u'rt_peak',u'rt_units', u'detected_polarity', u'mz', u'mz_tolerance',u'mz_tolerance_units','mono_isotopic_molecular_weight','pubchem_compound_id','synonyms','inchi','adduct']

    # atlas_keys = [u'label','compound_name','compound_description',u'synonyms', u'num_free_radicals', u'number_components', u'permanent_charge', u'rt_max', u'rt_min', u'rt_peak',
    #        u'rt_units', u'detected_polarity', u'mz', u'mz_tolerance',u'mz_tolerance_units',
    #         u'inchi', u'inchi_key', u'neutralized_2d_inchi', u'neutralized_2d_inchi_key', u'neutralized_inchi',
    #        u'neutralized_inchi_key',u'chebi_id', u'hmdb_id', u'img_abc_id', u'kegg_id',u'lipidmaps_id', u'metacyc_id',
    #        u'mono_isotopic_molecular_weight', u'pubchem_compound_id', u'kegg_url', u'chebi_url', u'hmdb_url', u'lipidmaps_url', u'pubchem_url',u'wikipedia_url',  u'source']

    atlas_df = atlas_df[atlas_keys]
    print((atlas_df.shape,len(label)))
    atlas_df['label'] = label
    return atlas_df

def get_data_for_mzrt(row,data_df_pos,data_df_neg,extra_time = 0.5,use_mz = 'mz',extra_mz = 0.0):
    min_mz = '(%s >= %5.4f & '%(use_mz,row.mz - row.mz*row.mz_tolerance / 1e6 - extra_mz)
    rt_min = 'rt >= %5.4f & '%(row.rt_min - extra_time)
    rt_max = 'rt <= %5.4f & '%(row.rt_max + extra_time)
    max_mz = '%s <= %5.4f)'%(use_mz,row.mz + row.mz*row.mz_tolerance / 1e6 + extra_mz)
    ms1_query_str = '%s%s%s%s'%(min_mz,rt_min,rt_max,max_mz)
    if row.detected_polarity == 'positive':
        if len(data_df_pos)>0:
            all_df = data_df_pos.query(ms1_query_str)
        else:
            return pd.Series()
    else:
        if len(data_df_neg)>0:
            all_df = data_df_neg.query(ms1_query_str)
        else:
            return pd.Series()
    return_df = pd.Series({'padded_feature_data':all_df,'in_feature':(all_df.rt >= row.rt_min) & (all_df.rt <= row.rt_max)})
    return return_df

def get_ms1_summary(row):
    #A DataFrame of all points typically padded by "extra time"
    all_df = row.padded_feature_data
    
    #slice out ms1 data that is NOT padded by extra_time
    ms1_df = all_df[(row.in_feature == True)]#[['i','mz','polarity','rt']]

    num_ms1_datapoints = ms1_df.shape[0]
    if num_ms1_datapoints > 0:
        idx = ms1_df.i.idxmax()
        ms1_peak_df = ms1_df.loc[ms1_df['i'].idxmax()]
        mz_peak = ms1_peak_df.mz
        rt_peak = ms1_peak_df.rt
        mz_centroid = sum(ms1_df.mz * ms1_df.i) / sum(ms1_df.i)
        rt_centroid = sum(ms1_df.rt * ms1_df.i) / sum(ms1_df.i)
        peak_height = ms1_peak_df.i
        peak_area = sum(ms1_df.i)
    else:
        mz_peak = np.nan
        rt_peak = np.nan
        mz_centroid = np.nan
        rt_centroid = np.nan
        peak_height = np.nan
        peak_area = np.nan
        
    return_df = pd.Series({ 'num_ms1_datapoints':num_ms1_datapoints,
                            'mz_peak':mz_peak,
                            'rt_peak':rt_peak,
                            'mz_centroid':mz_centroid,
                            'rt_centroid':rt_centroid,
                            'peak_height':peak_height,
                            'peak_area':peak_area})
    
    return return_df

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


def prefilter_ms1_dataframe_with_boundaries(data_df, rt_max, rt_min, mz_min, mz_max, extra_time = 0.5, extra_mz = 0.01):
    import math
    if (data_df.shape[0]==0) | (math.isnan(rt_max)):
        return []
    prefilter_query_str = 'rt <= %5.4f & rt >= %5.4f & mz >= %5.4f & mz <= %5.4f'%(rt_max+extra_time, rt_min-extra_time, mz_min-extra_mz, mz_max+extra_mz)
    new_df = data_df.query(prefilter_query_str)
    return new_df

def get_ms1_eic(row):
    #A DataFrame of all points typically padded by "extra time"
    all_df = row.padded_feature_data
    ms1_df = all_df[['i','mz','rt']]
    ms1_df = ms1_df.sort_values('rt',ascending=True)
    if ms1_df.shape[1] == 0:
        ms1_df = pd.DataFrame({'i','mz','rt'})
    return pd.Series({'eic':ms1_df.T})


def retrieve_most_intense_msms_scan(data):
    import numpy as np
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
    
    #filtered the ms2 and ms1 pos and neg frames in the container by rt and mz extreme points.
    filtered_ms1_pos = prefilter_ms1_dataframe_with_boundaries(df_container['ms1_pos'],
                                                               atlas_df[atlas_df.detected_polarity == 'positive'].rt_max.max(),
                                                               atlas_df[atlas_df.detected_polarity == 'positive'].rt_min.min(),
                                                               0,
                                                               #atlas_df[atlas_df.detected_polarity == 'positive'].mz.min()-1,
                                                               atlas_df[atlas_df.detected_polarity == 'positive'].mz.max()+1,
                                                               extra_time = extra_time,
                                                               extra_mz = extra_mz)
    filtered_ms1_neg = prefilter_ms1_dataframe_with_boundaries(df_container['ms1_neg'],
                                                           atlas_df[atlas_df.detected_polarity == 'negative'].rt_max.max(),
                                                           atlas_df[atlas_df.detected_polarity == 'negative'].rt_min.min(),
                                                           0,
                                                           #atlas_df[atlas_df.detected_polarity == 'negative'].mz.min()-1,
                                                           atlas_df[atlas_df.detected_polarity == 'negative'].mz.max()+1,
                                                           extra_time = extra_time,
                                                           extra_mz = extra_mz)
    filtered_ms2_pos = prefilter_ms1_dataframe_with_boundaries(df_container['ms2_pos'],
                                                           atlas_df[atlas_df.detected_polarity == 'positive'].rt_max.max(),
                                                           atlas_df[atlas_df.detected_polarity == 'positive'].rt_min.min(),
                                                           0,
                                                           #atlas_df[atlas_df.detected_polarity == 'positive'].mz.min()-1,
                                                           atlas_df[atlas_df.detected_polarity == 'positive'].mz.max()+1,
                                                           extra_time = extra_time,
                                                           extra_mz = extra_mz)
    
    filtered_ms2_neg = prefilter_ms1_dataframe_with_boundaries(df_container['ms2_neg'],
                                                           atlas_df[atlas_df.detected_polarity == 'negative'].rt_max.max(),
                                                           atlas_df[atlas_df.detected_polarity == 'negative'].rt_min.min(),
                                                           0,
                                                           #atlas_df[atlas_df.detected_polarity == 'negative'].mz.min()-1,
                                                           atlas_df[atlas_df.detected_polarity == 'negative'].mz.max()+1,
                                                           extra_time = extra_time,
                                                           extra_mz = extra_mz)
    

    ms1_feature_data = atlas_df.apply(lambda x: get_data_for_mzrt(x,filtered_ms1_pos,filtered_ms1_neg, extra_time=extra_time, extra_mz = extra_mz),axis=1)
    ms1_summary = ms1_feature_data.apply(get_ms1_summary,axis=1)
    #if ms1_summary.size == 0:
    #    return [],[],[]
    if ms1_feature_data.shape[1] == 0:
        return None,None,None
    else:
        ms1_eic = ms1_feature_data.apply(get_ms1_eic,axis=1)
    #print ms1_eic
        ms2_feature_data = atlas_df.apply(lambda x: get_data_for_mzrt(x,filtered_ms2_pos,filtered_ms2_neg,use_mz = 'precursor_MZ', extra_mz = extra_mz, extra_time=extra_time),axis=1)
        ms2_data = ms2_feature_data.apply(get_ms2_data,axis=1)
        dict_ms1_summary = [dict(row) for i,row in ms1_summary.iterrows()]
    
    dict_eic = []
    for i,row in ms1_eic.iterrows():
        dict_eic.append(row.eic.T.to_dict(orient='list'))
            
    #rename the "i" to "intensity".
    for i,d in enumerate(dict_eic):
        dict_eic[i]['intensity'] = dict_eic[i].pop('i')
    
    dict_ms2 = []
    for i,row in ms2_data.iterrows():
        if 'ms2_datapoints' in list(row.keys()):
            dict_ms2.append(row.ms2_datapoints.T.to_dict(orient='list'))
        else:
            dict_ms2.append([])
    
    return dict_ms1_summary,dict_eic,dict_ms2




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

def retrieve_most_intense_msms_scan(data):
    import numpy as np
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
    import numpy as np
    from metatlas.io import h5_query as h5q
    import tables
    
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
    import dill

    data = list()
    
    if os.path.exists(fname):
        with open(fname,'r') as f:
            try:
                data = dill.load(f)
            except IOError as e:
                print(("I/O error({0}): {1}".format(e.errno, e.strerror)))
            except:  # handle other exceptions such as attribute errors
                print(("Unexpected error:", sys.exc_info()[0]))


    return data


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
    import os.path

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
    from collections import defaultdict
    import re

    # if data is a string then it's a file name - get its data
    if isinstance(data, six.string_types):
        data = get_dill_data(data)

    compound_names = list()
    compound_objects = list()

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
        newstr = re.sub('\.', 'p', newstr) #2 or more in regexp

        newstr = re.sub('[\[\]]','',newstr)
        newstr = re.sub('[^A-Za-z0-9+-]+', '_', newstr)
        newstr = re.sub('i_[A-Za-z]+_i_', '', newstr)
        if newstr[0] == '_':
            newstr = newstr[1:]
        if newstr[0] == '-':
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

def make_data_sources_tables(groups, myatlas, output_loc, polarity=None):
    """
    polarity must be one of None, 'POS', 'NEG' or will throw ValueError
    """
    if polarity and not polarity in ['POS', 'NEG']:
        raise ValueError
    prefix = polarity + '_' if polarity else ''
    if not os.path.exists(output_loc):
        os.mkdir(output_loc)
    output_dir = os.path.join(output_loc,'data_sources')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    metob.to_dataframe([myatlas]).to_csv(os.path.join(output_dir, prefix+'atlas_metadata.tab'), sep='\t')
    metob.to_dataframe(groups).to_csv(os.path.join(output_dir, prefix+'groups_metadata.tab'), sep='\t')
    
    atlas_df = make_atlas_df(myatlas)
    atlas_df['label'] = [cid.name for cid in myatlas.compound_identifications]
    atlas_df.to_csv(os.path.join(output_dir,prefix+myatlas.name+'_originalatlas.tab'), sep='\t')

    group_path_df = pd.DataFrame(columns=['group_name','group_path','file_name'])
    loc_counter = 0
    for g in groups:
        for f in g.items:
            group_path_df.loc[loc_counter, 'group_name'] = g.name
            group_path_df.loc[loc_counter, 'group_path'] = os.path.dirname(f.mzml_file)
            group_path_df.loc[loc_counter, 'file_name'] = f.mzml_file
            loc_counter += 1

    group_path_df.to_csv(os.path.join(output_dir,prefix+'groups.tab'), sep='\t', index=False)
