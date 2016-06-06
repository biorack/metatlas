import numpy as np
import os.path
import sys
import tables
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
    from metatlas import h5_query as h5q
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
                print "I/O error({0}): {1}".format(e.errno, e.strerror)
            except:  # handle other exceptions such as attribute errors
                print "Unexpected error:", sys.exc_info()[0]


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
    if isinstance(data, basestring):
        data = get_dill_data(data)

    group_names = list()
    for i,d in enumerate(data):
        group_names.append(d[0]['group'].name)

    return group_names


def get_file_names(data):
    """
    Parameters
    ----------
    data: either a file name (str) or a list generated from loading the dill file

    Returns list containing the hdf file names present in the dill file
    -------
    """
    import os.path

    # if data is a string then it's a file name - get its data
    if isinstance(data, basestring):
        data = get_dill_data(data)

    file_names = list()
    for i,d in enumerate(data):
        file_names.append(os.path.basename(d[0]['lcmsrun'].hdf5_file))

    return file_names


def get_compound_names(data):
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
    if isinstance(data, basestring):
        data = get_dill_data(data)

    compound_names = list()
    compound_objects = list()

    for i,d in enumerate(data[0]):
        compound_objects.append(d['identification'])
        if len(d['identification'].compound) > 0:
            _str = d['identification'].compound[0].name
        else:
            _str = d['identification'].name
        newstr = '%s_%s_%s_%5.2f'%(_str,d['identification'].mz_references[0].detected_polarity,
                d['identification'].mz_references[0].adduct,d['identification'].rt_references[0].rt_peak)
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
