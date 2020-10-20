import os
import numpy as np
import pandas as pd
from scipy import interpolate
import time

#pandas columns that are "objects", but you are 100% sure contain strings
# will through this warning.  There is no way to set them as strings.
# pandas will permanently keep a reference to the object even if you 
# set it as a string.  Bottom line: make it a string and ignore the error.
# it prints to hdf5 just fine.
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

"""

"""

# SUMMARY OF TIMING AND MEMORY TESTING AT NERSC
# 10 threads
# on denovo, about 30 seconds each, a little over 9 minutes
# on denovo using scratchb, 20 seconds each, 7.43 minutes
# on cori, about 20 seconds each, 6.52 minutes staging all files to $SCRATCG=6.6 minutes

# 20 threads
# on denovo, about 30 seconds each, a little over 6 minutes
# on denovo scratch, 5.1 minutes
# on cori, didn't run on jupyter-dev, on cori withscratch = 5.4 minutes, 5.5 using realtime scratch, burst-buffer 5.23 minutes

# 30 threads
# on denovo, about 35 seconds each, the job never finished
# on cori, burst buffer: 4.1962 minutes repeated: 4.157 minutes
# on cori, scratch: 4.31 minutes

def setup_file_slicing_parameters(atlas,filenames,extra_time=0.1,ppm_tolerance=20,polarity='positive',project_dir=False,base_dir = '/project/projectdirs/metatlas/projects/',overwrite=True):
    """
    Make parameters that have to be setup to run the fast feature finding process.
    
    This function is called first when doing feature selection. It standardizes 
    all necessary inputs and files so downstream functions get a consistent place
    to work.
    
    Args:
        atlas (pandas dataframe): with [label,mz,rt_min,rt_max,rt_peak]. optional 
        parameters are fine too.

        filenames (list): full paths to hdf5 files

        extra_time (float): default=0.1 Time to get in addition to rt-min/max window (for making nice EICs)
        custom is to store metatlas hdf5 files in minutes, but always double check

        ppm_tolerance (float): default=20 Calibration is sometimes a problem. 20
        is safe.

        polarity (str): default='positive' or 'negative'

        project_dir (bool/str): default=False if user doesn't want to save their results
        or "path to output files"

        only relevant if project_dir is not False
        base_dir (str): '/project/projectdirs/metatlas/projects/'
        other paths that have been tried, but only a few percent faster than project:
            scratch_dir = os.environ['SCRATCH']
            scratch_dir = os.environ['DW_JOB_STRIPED']
    
    Returns:
        input_data list(dict): a list of python dictionaries with the following attributes
        for each lcmsrun to process:
            outfile (str): path to output hdf5 file of feature signals
            lcmsrun (str): lcmsrun to process
            atlas (pandas dataframe): atlas with necessary attributes for feature slicing,
            polarity (str): passthrough of input polarity string
            

    """

    """
    setup atlas table and define extra_time and ppm_tolerance in the atlas
    for compound atlases, label isn't stritly necessary add it if not provided
    """
    if not 'label' in atlas.columns:
        atlas['label'] = range(atlas.shape[0])
    atlas['extra_time'] = extra_time
    atlas['ppm_tolerance'] = ppm_tolerance

    """
    Group together m/z values that are within ppm_tolerance.
    This gives an index to acknowledge that there are multiple features with nearly 
    equal m/z.  Assigning it here speeds up the file slicing and feature selection 
    down the road.
    """
    atlas['group_index'] = group_consecutive(atlas['mz'].values[:],
                                             stepsize=ppm_tolerance,
                                             do_ppm=True)
    
    """
    define output directory
    """
    if project_dir is not False: #user doesn't want to save their results
        output_dir = os.path.join(base_dir,project_dir)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

    """
    get lcmsruns to process and build fullpath to output files
    """

    """
    setup input dictionary that will be the get_data input for each file
    """
    input_data = []
    for i,f in enumerate(filenames):
        #strip off the path and extension from the filename
        file_frag = ''.join(os.path.basename(f).split('.')[:-1])
        if len(file_frag)>0:
            output_filename = '%s_features.h5'%file_frag
            if project_dir is not False: #user doesn't want to save their results
                outfile = os.path.join(output_dir,output_filename)
            else:
                outfile = None
            input_data.append({'outfile':outfile,'lcmsrun':f,'atlas':atlas,'polarity':polarity}) # 'ppm_tolerance':ppm_tolerance,'extra_time':extra_time,,'start_time':time.time() 'file_index':i,

    """
    wipe out all the files and put the atlas in each one
    """
    if overwrite==True:
        for i in input_data:
            if i['outfile'] is not None: #user doesn't want to save their results
                with pd.HDFStore(i['outfile'],mode='w',complib='zlib',complevel=9) as f:
                    f.put('atlas',atlas,data_columns=True)

    return input_data


def group_consecutive(data,stepsize=10.0,do_ppm=True):
    """
    split a numpy array where consecutive elements are greater than stepsize
    can be ppm or value
    if idx is not None then returns indices of elements otherwise returns values
    
    The main use case is an unsorted list of m/zs as "data" and optionally 
    their numerical index as "idx". Typically a user would want to retrieve 
    the group indices in the original order that they provided their list 
    of m/zs.
    
    usage:
    
       
    """
    if type(data) is np.ndarray:
        # cool way to sort and unsort array:
        idx_sorted = data.argsort()
        sort_w_unsort = np.column_stack((np.arange(idx_sorted.size),idx_sorted))
        # sort_w_unsort[:,0] are the original indices of data
        # sort_w_unsort[:,1] are the sorted indices of data
        data_sorted = data[sort_w_unsort[:,1]]
        # np.argsort(sort_w_unsort[:,1]) returns the indices to map the sorted data back to the original
        # data_unsorted = data_sorted[np.argsort(sort_w_unsort[:,1])]
        
        if do_ppm:
            d = np.diff(data_sorted) / data_sorted[:-1] * 1e6
        else:
            d = np.diff(data_sorted)
        
        # make groups of the array
        data_groups = np.split(data_sorted, np.where(d > 2.0*stepsize)[0]+1)
        # replace each group of values with group index
        for i,data_slice in enumerate(data_groups):
            data_groups[i] = data_groups[i]*0 + i
        group_indices = np.concatenate(data_groups)
        # reorder the group indices
        group_indices = group_indices[np.argsort(sort_w_unsort[:,1])]
        return group_indices.astype(int)#
    else:
        print('not a numpy array. convert it and sort it first')
        
def map_mzgroups_to_data(mz_atlas,mz_group_indices,mz_data):
    """
    mz_atlas: m/z values from atlas
    mz_group_indices: integer index from "group_consecutive"
    mz_data: m/z values from raw data
    
    """
    from scipy import interpolate

    f = interpolate.interp1d(mz_atlas,np.arange(mz_atlas.size),kind='nearest',bounds_error=False,fill_value='extrapolate') #get indices of all mz values in the atlas
    idx = f(mz_data)   # iterpolate to find the nearest mz in the data for each mz in an atlas
    idx = idx.astype('int')
#     d = 1e6#np.abs(mz_data - mz_atlas[idx]) / mz_data * 1.0e6
#     output_mat = np.column_stack((d,))
    return mz_group_indices[idx]#output_mat

def df_container_from_metatlas_file(filename,desired_key=None):
    """
    
    """
    data_df = pd.DataFrame()

    pd_h5_file  = pd.HDFStore(filename)
        
    keys = pd_h5_file.keys()
    pd_h5_file.close()
    df_container = {}
    if desired_key is not None:
        return pd.read_hdf(filename,desired_key)
    else:
        for k in keys:
            if ('ms' in k) and not ('_mz' in k):
                new_df = pd.read_hdf(filename,k)
                df_container[k[1:]] = new_df
    return df_container






def group_duplicates(df,group_col,make_string=False,precision={'i':0,'mz':4,'rt':2}):
    """
    takes in a list of grouping columns and turns the rest into arrays
    """
    
    all_cols = np.asarray(df.columns)
    #get the index of the grouping term as array
    idx_group = np.argwhere(all_cols == group_col).flatten()
    #get the indices of all other terms as array
    idx_list = np.argwhere(all_cols != group_col).flatten()
    cols = all_cols[idx_list]
    
    # create a sorted numpy array (sorted by column=group_col)
    a = df.sort_values(group_col).values.T

    #get the indices of the first instance of each unique identifier
    ukeys, index = np.unique(a[idx_group,:],return_index=True)

    #split the other rows of the array into separate arrays using the 
    #unique index
    arrays = np.split(a[idx_list,:],index[1:],axis=1)

    #make a list of dicts with column headings as keys
    #if there are not multiple items then return value
    #If there are multiple items then return list
    
#     ucpds = [dict([(c,aa) if len(aa)>1 else (c,aa[0]) for c,aa in zip(cols,a)]) for a in arrays ]
    ucpds = [dict([(c,aa) for c,aa in zip(cols,a)]) for a in arrays ]

    #make a dataframe from the list of dicts
    df2 = pd.DataFrame(ucpds,index=ukeys)
    
    #make strings of array columns if you want to save it in anything useful
    if make_string==True:
        for c in cols:
#             df2[c] = df2[c].apply(lambda x: np.array2string(x, precision=5, separator=','))
            if c in precision.keys():
                pre_str = '{:.%df}'%precision[c]
            else:
                pre_str = '{:.4f}'
            df2[c] = df2[c].apply(lambda x: [pre_str.format(n) for n in x.tolist()])
#             df2[c] = df2[c].apply(lambda x: str(x.tolist()))
            
    df2.index = df2.index.set_names(group_col)
    df2.reset_index(inplace=True)

    #return dataframe
    return df2




def get_atlas_data_from_file(filename,atlas,desired_key='ms1_pos'):#,bundle=True,make_string=False):
    msdata = df_container_from_metatlas_file(filename,desired_key=desired_key)
    if 'ms2' in desired_key:
        # throw away all the intensity duplication here to make merging faster
        # this has the expense of having to remerge it later.
        msdata = msdata[['rt','precursor_MZ']].drop_duplicates('rt')
        msdata = msdata.rename(columns={'precursor_MZ':'mz'})
    
    g = map_mzgroups_to_data(atlas['mz'].values[:],
                               atlas['group_index'].values[:],
                               msdata['mz'].values[:])
    msdata['group_index'] = g#[:,1]
#     msdata['group_index_ppm'] = g[:,0]
    df = pd.merge(atlas,msdata,left_on='group_index',right_on='group_index',how='outer',suffixes=('_atlas','_data'))
    
    #grab all datapoints including "extra"
    mz_condition = abs(df['mz_data']-df['mz_atlas'])/df['mz_atlas']*1e6<df['ppm_tolerance']
    rt_min_condition = df['rt']>=(df['rt_min']-df['extra_time'])
    rt_max_condition = df['rt']<=(df['rt_max']+df['extra_time'])
    df = df[(mz_condition) & (rt_min_condition) & (rt_max_condition)]
    
    #label datapoints that are within the bounds of the feature vs "extra"
    df['in_feature'] = True
    if df['extra_time'].max()>0.0:
        cond_rt = (df['rt']<df['rt_min']) | (df['rt']>df['rt_max'])
        df.loc[cond_rt,'in_feature'] = False       
    
    #above, the df has mz_data and mz_atlas.  we don't need to differentiate anymore so:
    df = df.rename(columns={'mz_data':'mz'})
    if 'ms2' in desired_key:
        # keep in mind we don't have intensity or scan attributes
        df = df[['label','rt','in_feature']]
        # you've got to add it back in; so reload original file
        msdata = df_container_from_metatlas_file(filename,desired_key=desired_key)
        # This will merge back into the MSMS data 
        # the missing intensity and scan attributes
        mcols = ['rt','i','mz','precursor_MZ','precursor_intensity','collision_energy']
        df = pd.merge(df,msdata[mcols],left_on='rt',right_on='rt',how='left')
        return df.reset_index(drop=True)
    else:
        df = df[['label','rt','mz','i','in_feature']]
        return df.reset_index(drop=True)


def calculate_ms1_summary(row):
    """
    Calculate summary properties for features from data
    """
    d = {}
    #Before doing this make sure "in_feature"==True has already occured
    d['num_datapoints'] = row['i'].count()
    d['peak_area'] = row['i'].sum()
    idx = row['i'].idxmax()
    d['peak_height'] = row.loc[idx,'i']
    d['mz_centroid'] = sum(row['i']*row['mz'])/d['peak_area']
    d['rt_peak'] = row.loc[idx,'rt']
    return pd.Series(d)

# def calculate_ms1_summary(df):
#     a = df[['label','rt','mz','i','in_feature']].values
#     labels, row_pos = np.unique(a[:, 0], return_inverse=True) #these are feature labels
#     rt, col_pos = np.unique(a[:, 1], return_inverse=True) #these are rt values

#     pivot_table = np.zeros((len(labels), len(rt),3), dtype=float)
#     pivot_table[row_pos, col_pos] = a[:, [2,3,4]]
#     eic = pd.DataFrame(index=labels,data=pivot_table[:,:,1],columns=rt)
#     emzc = pd.DataFrame(index=labels,data=pivot_table[:,:,0],columns=rt)
#     efeaturec = pd.DataFrame(index=labels,data=pivot_table[:,:,2],columns=rt)
#     in_feature = efeaturec.values.astype(int)
#     intensity = np.multiply(eic.values,in_feature)
#     mz = np.multiply(emzc.values,in_feature)
#     rt = np.asarray(eic.columns)
#     labels = eic.index.tolist()
#     idx_max = np.argmax(intensity,axis=1)
#     df = pd.DataFrame(index=labels)
#     df['num_datapoints']=in_feature.sum(axis=1)
#     df['peak_area']=intensity.sum(axis=1)
#     df['peak_height']=np.diag(intensity[:,idx_max]) #I shouldn't have to do this and must be doing numpy slicing wrong!
#     df['mz_centroid']=np.divide(np.sum(np.multiply(mz,intensity),axis=1),intensity.sum(axis=1))
#     df['rt_peak']=rt[idx_max]
#     return df

def get_data(input_data,return_data=False,save_file=True):
    """
    Required Inputs a Dict that has these attributes:
    {'file_index':i, #a numerical index that helps with bookkeeping
    'outfile':outfile, #the hdf5 container to store the results
    'lcmsrun':new_file, #the hdf5 file corresponding to an lcms run
    'atlas':atlas, #the atlas dataframe containing minimally: [mz, rt_min,rt_max,rt_peak)]
    'ppm_tolerance':ppm_tolerance, #ppm tolerance in m/z
    'extra_time':extra_time} #time to add to the collected data beyond rt_min and rt_max
    
    The goal is to write to a file,
    ms1_data:
    
    ms1_summary:
    
    ms2_data:
    
    Returns a dictionary
    """
    out_data = {} #setup a container to store any data to return to the user otherwise save it to file
    polarity_short_string = input_data['polarity'][:3]
    d = get_atlas_data_from_file(input_data['lcmsrun'],input_data['atlas'],desired_key='ms1_%s'%polarity_short_string)#,bundle=True,make_string=True)
    if return_data is True:
        out_data['atlas'] = input_data['atlas']
        out_data['ms1_data'] = d
    if save_file is True:
        with pd.HDFStore(input_data['outfile'],mode='a',complib='zlib',complevel=9) as f:
            f.put('ms1_data',d,data_columns=True)
    
    d = d[d['in_feature']==True].groupby('label').apply(calculate_ms1_summary).reset_index()
    if d.shape[0]==0: #there isn't any data!
        for c in ['num_datapoints','peak_area','peak_height','mz_centroid','rt_peak']:
            d[c] = 0
        
    if return_data is True:
        out_data['ms1_summary'] = d
    if save_file is True:
        with pd.HDFStore(input_data['outfile'],mode='a',complib='zlib',complevel=9) as f:
            f.put('ms1_summary',d,data_columns=True)
    
#     input_data['atlas']['extra_time'] = 0.0 # set extratime here to be zero for msms getting
    d = get_atlas_data_from_file(input_data['lcmsrun'],input_data['atlas'],desired_key='ms2_%s'%polarity_short_string)#,bundle=True,make_string=True)
    if return_data is True:
        out_data['ms2_data'] = d
    if save_file is True:
        with pd.HDFStore(input_data['outfile'],mode='a',complib='zlib',complevel=9) as f:
            f.put('ms2_data',d,data_columns=True)
    
    if return_data is True:
        return out_data
