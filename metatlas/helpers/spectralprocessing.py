import re
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster


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


def sort_by_mz(a):
    """
    Returns array sorted by m/z with associated intensity
    
    :param a: numpy 2d array, v[0] is m/z values and v[1] is intensity values
    
    :return: numpy 2d array
    """
    
    return a.T[a.T[:, 0].argsort()].T


def align_ms_vectors(data, ref, distance, resolve="max"):
    """
    Clusters combined m/z values from data and reference, 
    creates 2xn data and reference vectors where n is the number of clusters with
    index corresponding to a cluster, and assigns m/z values and intensities
    to appropriate cluster index before returning data and reference vectors.
    
    Intensity values assigned to the same cluster index can be resolved in one of two ways:
    max or sum of the intensities. 
    M/z values assigned to the same cluster index use the m/z with the highest intensity.
    
    :param data: numpy 2d array, data[0] is m/z and data[1] is intensities sorted by m/z
    :param ref: numpy 2d array, ref[0] is m/z and ref[1] is intensities sorted by m/z
    :param distance: float, roughly +/- the size of the cluster
    :param resolve: ['max', 'sum'], how to resolve same cluster intensity values
    
    :return data_vector: numpy 2d array, data_vector[0] is m/z values and data_vector[1] is intensity values
    :return ref_vector: numpy 2d array, ref_vector[0] is m/z values and ref_vector[1] is intensity values
    """
    
    assert(resolve == 'max' or resolve == 'sum')
    
    #Trick clustering designed for 2D to use 1D by duplicating row
    combined_mz = np.array([np.append(data[0], ref[0]), np.append(data[0], ref[0])])
    
    #Cluster combined m/z values
    d = pdist(combined_mz.T)
    y = linkage(d,'average')
    clusters  = fcluster(y, distance, 'distance')

    #Create data and ref vectors
    data_vector = np.zeros([2, len(np.unique(clusters))])
    ref_vector = np.zeros([2, len(np.unique(clusters))])
    
    #Create counters to keep track of iteration through data and ref
    r = 0
    d = 0
    
    #Iterate over the clusters assigning m/z values 
    #and intensities to appropriate cluster index
    for i in range(clusters.size):
        if data[0].size > d and data[0][d] == combined_mz[0][i]:
            if data[1][d] > data_vector[1][clusters[i] - 1]:
                data_vector[0][clusters[i] - 1] += data[0][d]
                if resolve == 'max':
                    data_vector[1][clusters[i] - 1] = data[1][d]
            if resolve == 'sum':
                data_vector[1][clusters[i] - 1] += data[1][d]
            d += 1
        if ref[0].size > r and ref[0][r] == combined_mz[0][i]:
            if ref[1][r] > ref_vector[1][clusters[i] - 1]:
                ref_vector[0][clusters[i] - 1] = ref[0][r]
                if resolve == 'max': 
                    ref_vector[1][clusters[i] - 1] = ref[1][r]
            if resolve == 'sum':
                ref_vector[1][clusters[i] - 1] += ref[1][r]
            r += 1

    return data_vector, ref_vector


def weigh_vector_by_mz_and_intensity(v, mz_power=0, intensity_power=.6):
    """
    Returns (mz^mz_power)*(intensity^intensity_power) vector 
    
    :param v: numpy 2d array, v[0] is m/z values and v[1] is intensity values
    :param mz_power: float, power scales m/z
    :param intensity_power: float, power scales intensity
    
    :return: numpy 1d array
    """
    
    return np.multiply(np.power(v[0], mz_power), np.power(v[1], intensity_power))


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

    
def score_vectors_composite_dot(data_vector, ref_vector, mass_power=0, intensity_power=.6):
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
