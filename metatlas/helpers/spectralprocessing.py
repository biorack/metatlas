import re
import numpy as np

def remove_ms_vector_noise(ms_v, threshold = 500):
    return ms_v[:,ms_v[1,:] > threshold]

def sort_by_mz(ms_v):
    """
    Returns ms vector sorted by m/z with associated intensity
    
    :param ms_v: numpy 2d array, ms_v[0] is m/z values and ms_v[1] is intensity values
    
    :return: numpy 2d array
    """
    
    return ms_v[:,ms_v[0].argsort()]


def filter_frag_refs(metatlas_dataset, frag_refs, compound_idx, file_idx, condition):
    """
    Translates high-level condition into bitwise condition used to filter frag_refs 
    before returning dataframe of frag_refs meeting condition
    
    Keywords and their effects:
    inchi_key, 
        metatlas_dataset inchi_key matches frag_refs inchi_key
    neutralized_inchi_key, 
        metatlas_dataset neutralized_inchi_key matches frag_refs neutralized_inchi_key
    neutralized_2d_inchi_key, 
        metatlas_dataset neutralized_2d_inchi_key matches frag_refs neutralized_2d_inchi_key
    polarity,
        metatlas_dataset polarity matches frag_refs polarity
    rt, 
        median of metatlas_dataset measured retention time falls between min and max metatlas_dataset reference retention time
    precursor_mz, 
        absolute difference of metatlas_dataset precursor_mz and frag_refs precursor_mz
    collision_energy, 
        metatlas_dataset collision_energy matches frag_refs collision_energy
    
    Example:
    Filter by matching inchi_key or matching precursor_mz within .005 tolerance 
    and measured retention time within reference window and polarity
    '(inchi_key or (precursor_mz <= .005)) and rt and polarity'
    
    NOTE:
    All comparisons done at high-level likely require parentheses for bitwise condition to function. If in doubt, add more parentheses.
    
    :param metatlas_dataset:
    :param frag_refs: pandas dataframe, columns include 'inchi_key', 'neutralized_inchi_key', 'neutralized_2d_inchi_key', 'polarity', 'collision_energy', 'precursor_mz'
    :param file_idx: index of file in metatlas_dataset
    :param compound_idx: index of compound in metatlas_dataset
    :param condition: string, high-level condition to be evaluated
    
    :return: subset of frag_refs meeting condition
    """
    
    if 'data' in metatlas_dataset[file_idx][compound_idx]['data']['msms'].keys() and metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['mz'].size > 0:
        
        condition_list = list(filter(None, re.split(r"(\w+|[()])", condition)))
        
        for i in range(len(condition_list)):
            #Identifiers
            if metatlas_dataset[file_idx][compound_idx]['identification'].compound:
                if condition_list[i] == "inchi_key":
                    condition_list[i] = "(metatlas_dataset[file_idx][compound_idx]['identification'].compound[0].inchi_key == frag_refs['inchi_key'])"
                
                if condition_list[i] == "neutralized_inchi_key":
                    condition_list[i] = "(metatlas_dataset[file_idx][compound_idx]['identification'].compound[0].neutralized_inchi_key == frag_refs['neutralized_inchi_key'])"
                    
                if condition_list[i] == "neutralized_2d_inchi_key":
                    condition_list[i] = "(metatlas_dataset[file_idx][compound_idx]['identification'].compound[0].neutralized_2d_inchi_key == frag_refs['neutralized_2d_inchi_key'])"
                    
            else:
                if condition_list[i] == "inchi_key":
                    condition_list[i] = "(frag_refs.index != frag_refs.index)"
                
                if condition_list[i] == "neutralized_inchi_key":
                    condition_list[i] = "(frag_refs.index != frag_refs.index)"
                    
                if condition_list[i] == "neutralized_2d_inchi_key":
                    condition_list[i] = "(frag_refs.index != frag_refs.index)"
                    
            #Physical attributes
            if condition_list[i] == "polarity":
                condition_list[i] = "(metatlas_dataset[file_idx][compound_idx]['identification'].mz_references[0].detected_polarity == frag_refs['polarity'])"
                
            if condition_list[i] == "rt":
                    condition_list[i] = "(metatlas_dataset[file_idx][compound_idx]['identification'].rt_references[0].rt_min \
                                         <= np.median(metatlas_dataset[file_idx][compound_idx]['data']['eic']['rt']) <= \
                                         metatlas_dataset[file_idx][compound_idx]['identification'].rt_references[0].rt_max)"
            
            if condition_list[i] == "precursor_mz":
                    condition_list[i] = "(np.abs(metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['precursor_MZ'][0] - frag_refs['precursor_mz']))"
            
            if condition_list[i] == "collision_energy":
                    condition_list[i] = "(metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['collision_energy'][0] == frag_refs['collision_energy'])"
                    
            #Logic
            if condition_list[i] == "and":
                    condition_list[i] = "&"
            if condition_list[i] == "or":
                    condition_list[i] = "|"
            if condition_list[i] == "not":
                    condition_list[i] = "~"
                    
        return frag_refs[eval(''.join(condition_list))]
    
    else:
        return frag_refs[(frag_refs.index != frag_refs.index)]
    

def find_all_mz_matches(mz_v1, mz_v2, mz_tolerance):
    """
    Taken from pactolus and tweaked.
    """
    
    start_idxs = np.searchsorted(mz_v2, mz_v1-mz_tolerance)

    end_idxs = len(mz_v2) - np.searchsorted(-mz_v2[::-1], -(mz_v1+mz_tolerance))

    matches = [range(start, end) for start, end in zip(start_idxs, end_idxs)]

    unique_matches = np.unique(np.concatenate(matches))

    return matches, unique_matches
    
    
def match_ms_vectors(ms_v1, ms_v2, mz_tolerance, resolve_by):
    
    assert(resolve_by == 'distance' or resolve_by == 'intensity')
    
    matches = find_all_mz_matches(ms_v1[0], ms_v2[0], mz_tolerance)[0]
    
    if resolve_by == 'distance':
        match_matrix = np.fromfunction(np.vectorize(lambda i,j: np.absolute(ms_v1[0,int(i)] - ms_v2[0,int(j)]) if int(j) in matches[int(i)] else np.inf), (ms_v1[0].size, ms_v2[0].size)).astype(float)
    if resolve_by == 'intensity':
        match_matrix = np.fromfunction(np.vectorize(lambda i,j: ms_v1[1,int(i)] * ms_v2[1,int(j)] if int(j) in matches[int(i)] else -np.inf), (ms_v1[0].size, ms_v2[0].size)).astype(float)
        
    ms_v1_to_mv_v2 = np.array([np.nan for i in range(ms_v1[0].size)], dtype=float)
    ms_v2_to_ms_v1 = np.array([np.nan for i in range(ms_v2[0].size)], dtype=float)
    
    ms_v1_resolved = set([])
    ms_v2_resolved = set([])
    
    flat_match_matrix = match_matrix.ravel()
 
    if resolve_by == 'distance':
        masked_flat_match_matrix = np.ma.array(flat_match_matrix, mask=np.isinf(flat_match_matrix), fill_value=np.inf)
        sorted_match_indices = np.dstack(np.unravel_index(np.ma.argsort(masked_flat_match_matrix), (ms_v1[0].size, ms_v2[0].size)))[0]
    if resolve_by == 'intensity':
        masked_flat_match_matrix = np.ma.array(flat_match_matrix, mask=np.isinf(flat_match_matrix), fill_value=-np.inf)
        sorted_match_indices = np.dstack(np.unravel_index(np.ma.argsort(masked_flat_match_matrix, endwith=False)[::-1], (ms_v1[0].size, ms_v2[0].size)))[0]        
    
    for ms_v1_idx, ms_v2_idx in sorted_match_indices:
        if np.isinf(match_matrix[ms_v1_idx, ms_v2_idx]):
            break
        
        if ms_v1_idx not in ms_v1_resolved and ms_v2_idx not in ms_v2_resolved:
            ms_v1_resolved.add(ms_v1_idx)
            ms_v2_resolved.add(ms_v2_idx)
            
            ms_v1_to_mv_v2[ms_v1_idx] = ms_v2_idx
            ms_v2_to_ms_v1[ms_v2_idx] = ms_v1_idx
        
    return ms_v1_to_mv_v2, ms_v2_to_ms_v1


def partition_ms_vectors(ms_v1, ms_v2, mz_tolerance, resolve_by):
    """
    :param ms_v1: numpy 2d array, ms_v1[0] is m/z and ms_v1[1] is intensities sorted by m/z
    :param ms_v2: numpy 2d array, ms_v2[0] is m/z and ms_v2[1] is intensities sorted by m/z
    :param mz_tolerance: float, limits matches to be within +/- mz_tolerance
    :param resolve_by: ['intensity'], method to resolve conflicting matches into one to one matches
    
    :return ms_v1_matches: numpy 2d array, 
    :return ms_v2_matches: numpy 2d array,
    :return ms_v1_nonmatches: numpy 2d array, 
    :return ms_v2_nonmatches: numpy 2d array, 
    """
    
    #Find one to one matches within mz_tolerance
    ms_v1_to_ms_v2, ms_v2_to_ms_v1 = match_ms_vectors(ms_v1, ms_v2, mz_tolerance, resolve_by)

    #Assign matching ms vector components to same dimensions
    ms_v1_matches = ms_v1[:,ms_v2_to_ms_v1[~np.isnan(ms_v2_to_ms_v1)].astype(int)]
    ms_v2_matches = ms_v2[:,ms_v1_to_ms_v2[~np.isnan(ms_v1_to_ms_v2)].astype(int)]
    
    #Assign nonmatching ms vector components to differing dimensions
    ms_v1_nonmatches = ms_v1[:,np.isnan(ms_v1_to_ms_v2)]
    ms_v2_nonmatches = ms_v2[:,np.isnan(ms_v2_to_ms_v1)]
    
    return ms_v1_matches, ms_v2_matches, ms_v1_nonmatches, ms_v2_nonmatches


def partition_nl_vectors(ms_v1_precusor_mz, ms_v2_precusor_mz, ms_v1, ms_v2, mz_tolerance, resolve_by):
    return partition_ms_vectors(np.array([ms_v1[0] - (ms_v1_precusor_mz - ms_v2_precusor_mz), ms_v1[1]]), ms_v2, mz_tolerance, resolve_by)


def align_vectors(ms_v1_matches, ms_v2_matches, ms_v1_nonmatches, ms_v2_nonmatches):
    ms_v1_aligned = np.concatenate((ms_v1_matches, ms_v1_nonmatches, np.zeros(ms_v2_nonmatches.shape)), axis=1)
    ms_v2_aligned = np.concatenate((ms_v2_matches, np.zeros(ms_v1_nonmatches.shape), ms_v2_nonmatches), axis=1)
    
    return ms_v1_aligned, ms_v2_aligned


def weigh_vector_by_mz_and_intensity(ms_v, mz_power=1, intensity_power=.6):
    """
    Returns (mz^mz_power)*(intensity^intensity_power) vector 
    
    :param ms_v: numpy 2d array, ms_v[0] is m/z values and ms_v[1] is intensity values
    :param mz_power: float, power scales m/z
    :param intensity_power: float, power scales intensity
    
    :return: numpy 1d array
    """
    
    return np.multiply(np.power(ms_v[0], mz_power), np.power(ms_v[1], intensity_power))


def calc_ratio_of_pairs(v1, v2):
    """
    Returns "Ratio of Peak Pairs" as defined in Table 1 of paper below
    
    :param v1: numpy 1d array
    :param v2: numpy 1d array

    :return: float
    
    Stein, S. E., & Scott, D. R. (1994). 
    Optimization and testing of mass spectral library search algorithms for compound identification. 
    Journal of the American Society for Mass Spectrometry, 
    5(9), 859-866. doi:10.1016/1044-0305(94)87009-8
    """
    
    #Find common nonzero indices between v1 and v2
    shared = np.intersect1d(np.nonzero(v1), np.nonzero(v2))
    
    #Select shared indices from v1 and v2
    v1_shared = v1[shared]
    v2_shared = v2[shared]
    
    #Calculate ratio pairs
    ratios = np.multiply(np.divide(v1_shared[1:], v1_shared[:-1]), 
                         np.divide(v2_shared[:-1], v2_shared[1:]))
    
    #Inverse ratio if greater than one
    ratios = np.minimum(ratios, 1/ratios)
    
    return np.nan_to_num(np.sum(ratios)/shared.size)


def score_vectors_dot(v1, v2, normalize=True):
    """
    Returns dot product of two vectors, either normalized to a magnitude of 1 or not.
    
    :param v1: numpy 1d array
    :param v2: numpy 1d array 
    :param normalize: bool, scales to unit vector if true
    
    :return: 
    """
    
    if normalize:
        return np.dot(v1/np.linalg.norm(v1), v2/np.linalg.norm(v2))
    else:
        return np.dot(v1, v2)

    
def score_vectors_composite_dot(data_vector, ref_vector, mass_power=1, intensity_power=.6):
    """
    Returns "Composite" dot product score as defined in Table 1 of paper below
    
    :param data_vector: numpy 2d array, data[0] is m/z and data[1] is intensities
    :param ref_vector: numpy 2d array, ref[0] is m/z and ref[1] is intensities
    
    :return: float
    
    Stein, S. E., & Scott, D. R. (1994). 
    Optimization and testing of mass spectral library search algorithms for compound identification. 
    Journal of the American Society for Mass Spectrometry, 
    5(9), 859-866. doi:10.1016/1044-0305(94)87009-8
    """
    
    #Find number of peaks shared between data and reference
    num_shared = np.intersect1d(np.nonzero(data_vector), np.nonzero(ref_vector)).size
    
    #Weigh vectors according to mass and intensity
    weighted_data_vector = weigh_vector_by_mz_and_intensity(data_vector, mass_power, intensity_power)
    weighted_ref_vector = weigh_vector_by_mz_and_intensity(ref_vector, mass_power, intensity_power)

    #Compute dot product and ratio of pairs
    dot = score_vectors_dot(weighted_data_vector, weighted_ref_vector)
    ratio = calc_ratio_of_pairs(weighted_data_vector, weighted_ref_vector)
    
    return ((data_vector.size*dot) + (num_shared*ratio)) / (num_shared+data_vector.size)
