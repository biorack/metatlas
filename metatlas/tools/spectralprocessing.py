import sys
import re
from copy import deepcopy
from collections import Counter
import os
import numpy as np
from matchms.similarity import CosineHungarian
from matchms import Spectrum
from scipy.optimize import linear_sum_assignment
try:
    from numpy.fft.fftpack import fft,ifft
except:
    from numpy.fft import fft,ifft


from typing import List, Tuple, TypeAlias, Any
import numpy.typing as npt
import pandas as pd
from matplotlib import pyplot as plt
from six.moves import range
from six.moves import zip

# mass differences ranked by commonality in biochemical databases
mass_differences = pd.DataFrame([{"formula":"Same","mass":0.0,"rank":0},{"formula":"+O","mass":15.99491500,"rank":1.00000000},{"formula":"+H 2","mass":2.01565000,"rank":2.00000000},{"formula":"+H 2+O","mass":18.01056500,"rank":3.00000000},{"formula":"+H 2+C","mass":14.01565000,"rank":4.00000000},{"formula":"+H 2+C+O","mass":30.01056500,"rank":5.00000000},{"formula":"+H 8+C 5","mass":68.06260000,"rank":6.00000000},{"formula":"+O-H 2","mass":13.97926500,"rank":7.00000000},{"formula":"+H12+C 9-O 6","mass":24.12441000,"rank":8.00000000},{"formula":"+C","mass":12.00000000,"rank":9.00000000},{"formula":"+H10+C 6+O 5","mass":162.05282500,"rank":10.00000000},{"formula":"+C+O","mass":27.99491500,"rank":11.00000000},{"formula":"+O-H 3","mass":12.97144000,"rank":12.00000000},{"formula":"+O-H 2-C","mass":1.97926500,"rank":13.00000000},{"formula":"+O-H","mass":14.98709000,"rank":14.00000000},{"formula":"+O 2","mass":31.98983000,"rank":15.00000000},{"formula":"+H 4+C 2","mass":28.03130000,"rank":16.00000000},{"formula":"+H 8+C 5-O","mass":52.06768500,"rank":17.00000000},{"formula":"+H 4","mass":4.03130000,"rank":18.00000000},{"formula":"+O 2-H 3","mass":28.96635500,"rank":19.00000000},{"formula":"+O 3+P-H","mass":77.95068300,"rank":20.00000000},{"formula":"+H 6+C 5-O","mass":50.05203500,"rank":21.00000000},{"formula":"+H 2+C 2","mass":26.01565000,"rank":22.00000000},{"formula":"+C+O 2-H","mass":42.98200500,"rank":23.00000000},{"formula":"+H 4+C","mass":16.03130000,"rank":24.00000000},{"formula":"+H 2+C 2+O","mass":42.01056500,"rank":25.00000000},{"formula":"+O 2-H","mass":30.98200500,"rank":26.00000000},{"formula":"+H12+C 9-O 5","mass":40.11932500,"rank":27.00000000},{"formula":"+H 4+C 2+O","mass":44.02621500,"rank":28.00000000},{"formula":"+O-C","mass":3.99491500,"rank":29.00000000},{"formula":"+H 4+C 2-O","mass":12.03638500,"rank":30.00000000},{"formula":"+H26+C15+O","mass":222.19836500,"rank":31.00000000},{"formula":"+C 2","mass":24.00000000,"rank":32.00000000},{"formula":"+O 9+P-H13-C 9","mass":53.82627300,"rank":33.00000000},{"formula":"+O 2-H 2","mass":29.97418000,"rank":34.00000000},{"formula":"+H 4+C+O","mass":32.02621500,"rank":35.00000000},{"formula":"+O11-H 2-C 3","mass":137.92841500,"rank":36.00000000},{"formula":"+H","mass":1.00782500,"rank":37.00000000},{"formula":"+O-H 4","mass":11.96361500,"rank":38.00000000},{"formula":"+H 2+C+O 2","mass":46.00548000,"rank":39.00000000},{"formula":"+H16+C10","mass":136.12520000,"rank":40.00000000},{"formula":"+H 4+O","mass":20.02621500,"rank":41.00000000},{"formula":"+H14+C10-O 5","mass":54.13497500,"rank":42.00000000},{"formula":"+C+O-H","mass":26.98709000,"rank":43.00000000},{"formula":"+C-H 2","mass":9.98435000,"rank":44.00000000},{"formula":"+O 5-H 4-C 4","mass":27.94327500,"rank":45.00000000},{"formula":"+H 2+C 2-O","mass":10.02073500,"rank":46.00000000},{"formula":"+H+O","mass":17.00274000,"rank":47.00000000},{"formula":"+H 8+C 5+O","mass":84.05751500,"rank":48.00000000},{"formula":"+H 3","mass":3.02347500,"rank":49.00000000},{"formula":"+H 2+O 2","mass":34.00548000,"rank":50.00000000}])


# typing
MS2Spectrum: TypeAlias = npt.NDArray[npt.NDArray[np.float64]]
MS2MatchedPeaks: TypeAlias = List[Tuple[Tuple[int, int], float]]
MS2Aligned = Tuple[npt.NDArray[npt.NDArray[Any]], npt.NDArray[npt.NDArray[Any]]]


################################################################################
#Misc Functions
################################################################################

def make_feature_label(row,polarity_attr='polarity',mz_attr='mz',rt_attr='rt_peak'):
    """
    For consistency these are my preferred way to do this:
    * polarity: a "+" or "-" sign
    * mz: float (typically mz_centroid)
    * rt: float (typically rt_peak)

    Returns:
    * string that looks like "+1234.1234@1.23" for
    positive mode, m/z=1234.1234, & rt=1.23
    """
    return '%s%.4f@%.2f'%(row[polarity_attr],row[mz_attr],row[rt_attr])

def make_edges(mz_vec,label_vec,delta_mass,tolerance=5,edge_label=''):
    """
    Create an edge list based on chemical differences.

    Inputs:
    * mz_vec is a list or array of m/z values
    * label_vec is a list or array of strings to label nodes
    * delta mass is a float that is the difference in masses between two features
    * tolerance is in ppm
    * edge_label (str) is usually chemical formula and is an attribute of mass differences.

    Outputs:
    * edge_list: a list of dicts with the keys
       - source and target identified by the label_vec
       - edge labeled by "edge_label"
       - weight: calculated by the 1/num_edges for a particular difference

    See Also: metatlas.spectralprocessing.make_feature_label() for a consistent way to label features.

    Usage:
    # df is a typical table with mz and label values
    # it should have "polarity",""
    df['label'] = df.apply(make_feature_label,axis=1)
    all_edges = []
    for i in range(10): #do the top 10 mass differences
        e = make_edges(df['mz'].values,
                   df['label'].values,
                   mass_differences.loc[i,'mass'],
                  tolerance=3,
                   edge_label=mass_differences.loc[i,'formula'])
        all_edges.extend(e)
    """

    mz_vec = np.asarray(mz_vec)
    label_vec = np.asarray(label_vec)
    delta_mass = float(delta_mass)
    tolerance = float(tolerance)

    m = np.abs(np.abs(mz_vec[:,np.newaxis] - mz_vec) - delta_mass)
    e = (mz_vec[:,np.newaxis] + mz_vec) * tolerance / 2. / 1000000.
    idx = np.argwhere(m<e)
    idx = idx[idx[:,0]<idx[:,1],:] #keep lower left triangle
    zipped_lists = list(zip(label_vec[idx[:,0]],label_vec[idx[:,1]]))
    edge_list = [{'source' : v1, 'target' : v2, 'edge':edge_label} for v1,v2 in zipped_lists]
    for item in edge_list:
        item.update({'weight':1.0/float(len(e))})
    return edge_list

def histogram_edge_frequency(edge_count):
    """
    Input:
    * edge_count: list of dicts generated by metatlas.spectralprocessing.make_edges()

    Output:
    * fig, ax: a figure and axes handles for the frequency of occurnce of various edges in your network

    See Also:
    metatlas.spectralprocessing.make_edges()

    """
    plt.ioff()
    ec_df = pd.DataFrame(edge_count)
    sorted_mass_difference = pd.merge(ec_df,mass_differences,left_on='edge_index',right_index=True).drop(columns=['edge_index']).sort_values('edge_count',ascending=False)
    sorted_mass_difference.reset_index(inplace=True)
    fig,ax = plt.subplots(figsize=(14,7))
    sorted_mass_difference.plot.bar(x='formula',y='edge_count',ax=ax,color='black')
    ax.set_ylabel('edge count')
    plt.ion()
    return fig,ax

def transform_network_type(input_list,output_type=('dataframe','networkx'),do_mst=False):
    """
    For a input:
    * list of dicts with attributes (source,target,edge,weight)
    Other inputs:
    * do_mst: if True and output_type is networkx, returns a minimum spanning tree network

    Return:
    * for output_type='dataframe': a pandas dataframe each row is a connection
    * for output_type='networkx': a networkx graph with corresponding node and edge attributes

    See Also:
    * metatlas.spectralprocessing.make_edges()

    """
    if output_type=='dataframe':
        return pd.DataFrame(input_list)
    elif output_type=='networkx':
        import networkx as nx
        G = nx.Graph()
        for e in input_list:
            G.add_edge(e['source'],e['target'], weight=e['weight'],difference=e['edge'])
        if do_mst==True:
            G = nx.minimum_spanning_tree(G)
            return G

def write_network_for_cytoscape(G,output_file,edge_label='difference',weight_label='weight'):
    """
    Takes in a networkx graph and outputs an edgelist.csv that cytoscape can import
    """
    import networkx as nx
    nx.write_edgelist(G, output_file, data=[edge_label,weight_label],delimiter=',')


def remove_ms_vector_noise(msv, threshold=1e-3):
    """
    Returns ms vector with intensity below threshold removed

    :param msv: numpy 2d array, msv[0] is m/z values and msv[1] is intensity values
    :param threshold: float

    :return: numpy 2d array
    """

    return msv[:, msv[1, :] > threshold * np.max(msv[1])]


def sort_ms_vector_by_mz(msv):
    """
    Returns ms vector sorted by m/z

    :param msv: numpy 2d array, msv[0] is m/z values and msv[1] is intensity values

    :return: numpy 2d array
    """

    return msv[:, msv[0].argsort()]


def filter_frag_refs(metatlas_dataset, frag_refs, compound_idx, file_idx,
                     condition):
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

    :return: pandas dataframe, subset of frag_refs meeting condition
    """

    if 'data' in list(metatlas_dataset[file_idx][compound_idx]['data']['msms'].keys()) and \
                 metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['mz'].size > 0:

        condition_list = list([_f for _f in re.split(r"(\w+|[()])", condition) if _f])

        for i in range(len(condition_list)):
            # Identifiers
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

            # Physical attributes
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

            # Logic
            if condition_list[i] == "and":
                condition_list[i] = "&"
            if condition_list[i] == "or":
                condition_list[i] = "|"
            if condition_list[i] == "not":
                condition_list[i] = "~"

        return frag_refs[eval(''.join(condition_list))]

    else:
        return frag_refs[(frag_refs.index != frag_refs.index)]

def combine_ms_vectors(msv_alignment, combine_by, weights=None):
    """
    Returns ms vector that is the combination of ms vectors in msv_alignment
    with combined m/z being the mean m/z and intensity being either the mean
    or sum of the intensities

    :param msv_alignment: numpy 3d array,
        msv_alignment[i] is ms vector at position i,
        msv_alignment[i,0] is  m/z values at position i,
        msv_alignment[i,1] is intensity values at position i
    :param combine_by: 'mean', 'max', or 'sum': dictates how to combine intensities
    :param weights: list of floats of length msv_alignment[i,j],
        weighs msv_alignment[i] with weights[i] when computing mean

    :return: numpy 2d array
    """

    assert combine_by == 'mean' or combine_by == 'max' or combine_by == 'sum'

    def nanaverage(v, weights, axis):
        v = np.ma.masked_invalid(v)
        v = np.ma.average(v, weights=weights, axis=axis)
        v.mask = np.ma.nomask
        return v

    if combine_by == 'mean':
        return sort_ms_vector_by_mz(np.array([nanaverage(msv_alignment[:,0], weights, 0),
                                              nanaverage(msv_alignment[:,1], weights, 0)]))

    if combine_by == 'max':
        return sort_ms_vector_by_mz(np.array([nanaverage(msv_alignment[:,0], weights, 0),
                                              np.nanmax(msv_alignment[:,1], axis = 0)]))

    if combine_by == 'sum':
        return sort_ms_vector_by_mz(np.array([nanaverage(msv_alignment[:,0], weights, 0),
                                              np.nansum(msv_alignment[:,1], axis = 0)]))

# TODO
# def lower_ms_vector_resolution(msv, mz_tolerance):
#     """
#     Finds indices where m/z values in msv_1 and msv_2 'best' match within +/- mz_tolerance
#     by finding all matches (see find_all_ms_matches) and resolving non-one-to-one cases by
#     assigning matches in order of:
#         'distance': smallest to largest difference in matching m/zs
#         'shape': smallest to largest difference in matching normalized intensities
#         'intensity': largest to smallest product of matching intensities
#     Returns two 1d numpy arrays with indices corresponding to indices of msv_1 and msv_2
#     respectively where each value is the matching index in the other ms vector or nan if
#     there is no matching index.
#
#     :param msv_1: numpy 2d array, msv_1[0] is m/z values and msv_1[1] is intensity values
#     :param msv_2:  numpy 2d array, msv_2[0] is m/z values and msv_2[1] is intensity values
#     :param mz_tolerance: limits matches to be within +/- mz_tolerance
#     :param resolve_by: 'distance', 'shape', or 'intensity',
#         determines what to prioritize when there are multiple matching m/zs within mz_tolerance
#
#     :return: msv_1_to_msv_2, numpy_1d_array of length msv_1
#     :return: msv_2_to_msv_1, numpy_1d_array of length msv_2
#     """
#
#     msv_1 = np.asarray(msv, dtype=float)
#     assert np.any(np.isnan(msv)) == False
#     assert msv.shape[0] == 2
#
#     matches = find_all_ms_matches(msv, msv, mz_tolerance)[0]
#
#     # Create match matrix where row i and column j is value determined by resolve_by if msv_1[i] matches msv_2[j]
#     # and a default value if otherwise
#     match_matrix = np.fromfunction(np.vectorize(lambda i, j: np.absolute(msv[0, int(i)] - msv[0, int(j)]) if int(j) in matches[int(i)] else np.inf), (msv[0].size, msv[0].size)).astype(float)
#
#     # Initialize numpy arrays to return matching indices
#     msv_1_to_msv_2 = np.full_like(msv_1[0], np.nan, dtype=float)
#     msv_2_to_msv_1 = np.full_like(msv_2[0], np.nan, dtype=float)
#
#     # Flatten match matrix
#     flat_match_matrix = match_matrix.ravel()
#
#     # Sort match matrix indices in order determined by resolve_by
#     if resolve_by == 'distance' or resolve_by == 'shape':
#         masked_flat_match_matrix = np.ma.array(flat_match_matrix,
#                                                mask=np.isinf(flat_match_matrix),
#                                                fill_value=np.inf)
#         sorted_match_indices = np.dstack(
#             np.unravel_index(np.ma.argsort(masked_flat_match_matrix),
#                              (msv_1[0].size, msv_2[0].size)))[0]
#
#     # Assign matching indices in order determined by resolve_by
#     for msv_1_idx, msv_2_idx in sorted_match_indices:
#         if np.isinf(match_matrix[msv_1_idx, msv_2_idx]):
#             break
#
#         if msv_1_idx not in msv_2_to_msv_1 and msv_2_idx not in msv_1_to_msv_2:
#             msv_1_to_msv_2[msv_1_idx] = msv_2_idx
#             msv_2_to_msv_1[msv_2_idx] = msv_1_idx
#
#     return msv_1_to_msv_2, msv_2_to_msv_1


################################################################################
#Score Related Functions
################################################################################

def find_all_ms_matches(msv_1, msv_2, mz_tolerance):
    """
    Taken from pactolus and tweaked. Finds all indices where m/z values in
    msv_1 and msv_2 match within +/- mz_tolerance and returns list of lists of
    integers where index i of the outer list corresponds to the m/z value of
    msv_1 at index i matching the m/z values of msv_2 at indices given by the
    integers in the inner list at index i,
    and unique indices of msv_2 whose m/z values match m/z values in msv_1.

    :param msv_1: numpy 2d array, msv_1[0] is m/z values and msv_1[1] is intensity values
    :param msv_2:  numpy 2d array, msv_2[0] is m/z values and msv_2[1] is intensity values
    :param mz_tolerance: m/z tolerance for m/z values to be considering matching

    :return: matches, list of lists
    :return: unique_matches, numpy 1d array
    """

    msv_1 = np.asarray(msv_1, dtype=float)
    msv_2 = np.asarray(msv_2, dtype=float)
    assert not np.any(np.isnan(msv_1)) and not np.any(np.isnan(msv_2))
    assert msv_1.shape[0] == 2 and msv_2.shape[0] == 2

    #start_idx is element for which inserting msv_1 directly ahead of it maintains sort order
    start_idxs = np.searchsorted(msv_2[0], msv_1[0] - mz_tolerance)

    #end_idx is element for which inserting msv_1 directly after it maintains sort order
    #found by searching negative reversed list since np.searchshorted requires increasing order
    end_idxs = len(msv_2[0]) - np.searchsorted(-msv_2[0][::-1], -(msv_1[0] + mz_tolerance))

    #if the start and end idx is the same, the peak is too far away in mass from the data and will be empty
    matches = [list(range(start, end)) for start, end in zip(start_idxs, end_idxs)]

    #flattening the list
    unique_matches = np.unique(np.concatenate(matches))

    return matches, unique_matches


def match_ms_vectors(msv_1, msv_2, mz_tolerance, resolve_by):
    """
    Finds indices where m/z values in msv_1 and msv_2 'best' match within +/- mz_tolerance
    by finding all matches (see find_all_ms_matches) and resolving non-one-to-one cases by
    assigning matches in order of:
        'distance': smallest to largest difference in matching m/zs
        'shape': smallest to largest difference in matching normalized intensities
        'intensity': largest to smallest product of matching intensities
    Returns two 1d numpy arrays with indices corresponding to indices of msv_1 and msv_2
    respectively where each value is the matching index in the other ms vector or nan if
    there is no matching index.

    :param msv_1: numpy 2d array, msv_1[0] is m/z values and msv_1[1] is intensity values
    :param msv_2:  numpy 2d array, msv_2[0] is m/z values and msv_2[1] is intensity values
    :param mz_tolerance: limits matches to be within +/- mz_tolerance
    :param resolve_by: 'distance', 'shape', or 'intensity',
        determines what to prioritize when there are multiple matching m/zs within mz_tolerance

    :return: msv_1_to_msv_2, numpy_1d_array of length msv_1
    :return: msv_2_to_msv_1, numpy_1d_array of length msv_2
    """

    msv_1 = np.asarray(msv_1, dtype=float)
    msv_2 = np.asarray(msv_2, dtype=float)
    assert not np.any(np.isnan(msv_1)) and not np.any(np.isnan(msv_2))
    assert msv_1.shape[0] == 2 and msv_2.shape[0] == 2
    assert resolve_by in ['distance', 'shape', 'intensity']

    matches = find_all_ms_matches(msv_1, msv_2, mz_tolerance)[0]

    # Create match matrix where row i and column j is value determined by resolve_by
    # if msv_1[i] matches msv_2[j] and a default value if otherwise
    def match_distance(i, j):
        return np.absolute(msv_1[0, int(i)] - msv_2[0, int(j)]) if j in matches[int(i)] else np.inf

    def match_shape(i, j):
        ms_scale = np.median([msv_1[1, i]/msv_2[1, j] for i, l in enumerate(matches) for j in l])
        return np.absolute(msv_1[1, int(i)] - (ms_scale*msv_2[1, int(j)])) if j in matches[int(i)] else np.inf

    def match_intensity(i, j):
        return msv_1[1, int(i)] * msv_2[1, int(j)] if j in matches[int(i)] else -np.inf

    match_fun = {"distance": match_distance, "shape": match_shape, "intensity": match_intensity}[resolve_by]
    match_matrix = np.fromfunction(np.vectorize(match_fun, otypes=[float]), (msv_1[0].size, msv_2[0].size))

    # Initialize numpy arrays to return matching indices
    msv_1_to_msv_2 = np.full_like(msv_1[0], np.nan, dtype=float)
    msv_2_to_msv_1 = np.full_like(msv_2[0], np.nan, dtype=float)

    # Flatten match matrix
    flat_match_matrix = match_matrix.ravel()

    # Sort match matrix indices in order determined by resolve_by
    if resolve_by == 'distance' or resolve_by == 'shape':
        masked_flat_match_matrix = np.ma.array(flat_match_matrix,
                                               mask=np.isinf(flat_match_matrix),
                                               fill_value=np.inf)
        sorted_match_indices = np.dstack(
            np.unravel_index(np.ma.argsort(masked_flat_match_matrix),
                             (msv_1[0].size, msv_2[0].size)))[0]
    if resolve_by == 'intensity':
        masked_flat_match_matrix = np.ma.array(flat_match_matrix,
                                               mask=np.isinf(flat_match_matrix),
                                               fill_value=-np.inf)

        sorted_match_indices = np.dstack(np.unravel_index(
            np.ma.argsort(masked_flat_match_matrix, endwith=False)[::-1],
            (msv_1[0].size, msv_2[0].size)))[0]

    # Assign matching indices in order determined by resolve_by
    for msv_1_idx, msv_2_idx in sorted_match_indices:
        if np.isinf(match_matrix[msv_1_idx, msv_2_idx]):
            break

        if msv_1_idx not in msv_2_to_msv_1 and msv_2_idx not in msv_1_to_msv_2:
            msv_1_to_msv_2[msv_1_idx] = msv_2_idx
            msv_2_to_msv_1[msv_2_idx] = msv_1_idx

    return msv_1_to_msv_2, msv_2_to_msv_1


def partition_ms_vectors(msv_1, msv_2, mz_tolerance, resolve_by):
    """
    Finds which m/z values in msv_1 and msv_2 'best' match within +/- mz_tolerance
    (see match_ms_vectors) and returns 4 numpy 2d arrays with msv_1 and msv_2 partitioned
    into matches and nonmatches. Index i of matching partitions match with index i of the
    other and all indices of nonmatching partitions do not match any index of the other.

    :param msv_1: numpy 2d array, msv_1[0] is m/z and msv_1[1] is intensities sorted by m/z
    :param msv_2: numpy 2d array, msv_2[0] is m/z and msv_2[1] is intensities sorted by m/z
    :param mz_tolerance: float, limits matches to be within +/- mz_tolerance
    :param resolve_by: 'distance', 'shape', or 'intensity',
        determines what to prioritize when there are multiple matching m/zs within mz_tolerance
        (see match_ms_vectors)

    :return msv_1_matches: numpy 2d array,
    :return msv_2_matches: numpy 2d array,
    :return msv_1_nonmatches: numpy 2d array,
    :return msv_2_nonmatches: numpy 2d array,
    """

    msv_1 = np.asarray(msv_1, dtype=float)
    msv_2 = np.asarray(msv_2, dtype=float)
    assert (np.any(np.isnan(msv_1)) and np.any(np.isnan(msv_2))) == False
    assert msv_1.shape[0] == 2 and msv_2.shape[0] == 2

    # Find one to one matches within mz_tolerance
    msv_1_to_msv_2, msv_2_to_msv_1 = match_ms_vectors(msv_1, msv_2, mz_tolerance, resolve_by)

    # Assign matching ms vector components to same dimensions
    msv_1_matches = msv_1[:,~np.isnan(msv_1_to_msv_2)]
    msv_2_matches = msv_2[:,msv_1_to_msv_2[~np.isnan(msv_1_to_msv_2)].astype(int)]

    # Assign nonmatching ms vector components to differing dimensions
    msv_1_nonmatches = msv_1[:,np.isnan(msv_1_to_msv_2)]
    msv_2_nonmatches = msv_2[:,np.isnan(msv_2_to_msv_1)]

    return msv_1_matches, msv_2_matches, msv_1_nonmatches, msv_2_nonmatches


def partition_nl_vectors(msv_1_precusor_mz, msv_2_precusor_mz, msv_1, msv_2, mz_tolerance, resolve_by):
    return partition_ms_vectors(
        np.array([msv_1[0] - (msv_1_precusor_mz - msv_2_precusor_mz),msv_1[1]]),
        msv_2,
        mz_tolerance, resolve_by)


def match_peaks(spec1: MS2Spectrum, spec2: MS2Spectrum, frag_mz_tolerance: float) -> MS2MatchedPeaks:
    """
    Match MS2 fragment peaks within m/z tolerance.

    Outputs both the coordinates of the match in addition to the product of their matched intensities.
    The product of matched intensities is maximized when choosing between two matched peaks within tolerance.
    """
    matched_peaks = []
    for i in range(spec1.shape[1]):
        tolerance_filter = np.isclose(spec1[0][i], spec2[0], atol=frag_mz_tolerance)
        if tolerance_filter.any():
            for matching_idx in np.argwhere(tolerance_filter):
                match_coords = (i, matching_idx.item())
                match_value = (spec1[1][i].item() * spec2[1][matching_idx].item())
                matched_peaks.append((match_coords, match_value))

    return matched_peaks


def hungarian_assignment(matched_peaks: MS2MatchedPeaks, matrix_size: int) -> List[Tuple[int, int]]:
    """
    Filter matched peaks by maximizing the matched intensity product with the hungarian algorithm.

    This approach is similarly implemented in MatchMS for reconciling multiple matched peaks within tolerance.
    """

    cost_matrix = np.zeros((matrix_size, matrix_size))
    for match in matched_peaks:
        cost_matrix[match[0][0], match[0][1]] = match[1]

    row_idx, col_idx = linear_sum_assignment(cost_matrix, maximize=True)
    optimized_coords = [(row_idx[i].item(), col_idx[i].item()) for i in range(row_idx.shape[0])]

    filtered_coords = [coord[0] for coord in matched_peaks if coord[0] in optimized_coords]

    return filtered_coords


def link_alligned_spectra(spec1: MS2Spectrum, spec2: MS2Spectrum,
                          filtered_coords: List[Tuple[int, int]]) -> MS2Aligned:
    """
    Create linked and alligned MS2 spectra using filtered matching fragment indices.

    The output objects from this function are used for generating mirror plots elswhere in the codebase.

    Aligned and linked spectrum 1 is modeled as a numpy array with NaN values used as
    placeholders for unmatched peaks from spectrum 2, followed by all peak data from spectrum 1.

    Aligned and linked spectrum 2 is modeled as a numpy array with all unshared spectrum 2
    peak data, followed by matched spectrum 2 m/zs with NaN values interspersed as placeholders.

    Simple example:
        # two shared peaks: 1 and 2
        # two arrays, first for m/z values and second for intensities
        spec1 = [[1, 2, 3, 4], [5, 5, 5, 10]]
        spec2 = [[1, 2, 5, 6, 7], [5, 5, 5, 5, 10]]

        # nan placeholders for 3 unmatched peaks from spec2
        msv_spec1_aligned = [[nan, nan, nan, 1, 2, 3, 4],
                             [nan, nan, nan, 5, 5, 5, 10]]
        # peak values for 3 unmatched peak from spec2
        # peak values for matched m/z, then nan for placeholder
        msv_spec2_aligned = [[5, 6, 7, 1, 2, nan, nan],
                             [5, 5, 10, 5, 5, nan, nan]]
    """

    shared_spec1_idxs = [coord[0] for coord in filtered_coords]
    shared_spec2_idxs = [coord[1] for coord in filtered_coords]

    shared_spec2_mzs = np.array([spec2[0][i] for i in range(spec2.shape[1]) if i in shared_spec2_idxs])
    shared_spec2_is = np.array([spec2[1][i] for i in range(spec2.shape[1]) if i in shared_spec2_idxs])

    unshared_spec2_mzs = np.array([spec2[0][i] for i in range(spec2.shape[1]) if i not in shared_spec2_idxs])
    unshared_spec2_is = np.array([spec2[1][i] for i in range(spec2.shape[1]) if i not in shared_spec2_idxs])

    spec1_alignment_linker = np.full(unshared_spec2_mzs.shape, np.nan)
    msv_spec1_mzs = np.concatenate((spec1_alignment_linker, spec1[0]))
    msv_spec1_is = np.concatenate((spec1_alignment_linker, spec1[1]))
    msv_spec1_aligned = np.asarray([msv_spec1_mzs, msv_spec1_is], dtype=float)

    spec2_alignment_linker = np.full(spec1[0].shape, np.nan)
    np.put(spec2_alignment_linker, shared_spec1_idxs, shared_spec2_mzs)
    msv_spec2_mzs = np.concatenate((unshared_spec2_mzs, spec2_alignment_linker))

    np.put(spec2_alignment_linker, shared_spec1_idxs, shared_spec2_is)
    msv_spec2_is = np.concatenate((unshared_spec2_is, spec2_alignment_linker))
    msv_spec2_aligned = np.asarray([msv_spec2_mzs, msv_spec2_is], dtype=float)

    return msv_spec1_aligned, msv_spec2_aligned


def pairwise_align_ms_vectors(msv_query: MS2Spectrum, msv_ref: MS2Spectrum,
                              frag_mz_tolerance: float) -> MS2Aligned:
    """
    Generate linked and aligned spectra from MS2 spectral vectors.
    """

    matched_peaks = match_peaks(msv_query, msv_ref, frag_mz_tolerance)
    filtered_coords = hungarian_assignment(matched_peaks, max(msv_query.shape[1], msv_ref.shape[1]))
    msv_query_aligned, msv_ref_aligned = link_alligned_spectra(msv_query, msv_ref, filtered_coords)

    return msv_query_aligned, msv_ref_aligned


def align_all_ms_vectors(msv_pairs: Tuple[MS2Spectrum, MS2Spectrum], frag_tolerance: float) -> List[MS2Aligned]:
    """
    Generate list of linked and aligned spectra given pairs of MS2 spectral vectors.
    """
    all_queries = []
    all_refs = []

    for pair in msv_pairs:
        msv_query_aligned, msv_ref_aligned = pairwise_align_ms_vectors(pair[0], pair[1], frag_tolerance)

        all_queries.append(msv_query_aligned)
        all_refs.append(msv_ref_aligned)

    return all_queries, all_refs


def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)

# def multiple_tree_align_ms_vectors(msv_list, mz_tolerance, resolve_by):
#     msv_list = [np.asarray(msv, dtype=float) for msv in msv_list]
#     assert ~np.any([np.any(np.isnan(msv)) for msv in msv_list])
#     assert np.all([msv.shape[0] == 2 for msv in msv_list])
#
#     score_graph = nx.Graph()
#
#     for i in range(len(msv_list)):
#         score_graph.add_node(i)
#
#     for i, j in np.array(np.tril_indices(len(msv_list), -1)).T:
#         msv_i_aligned, msv_j_aligned = pairwise_align_ms_vectors(msv_list[i], msv_list[j], mz_tolerance, resolve_by)
#         score = score_ms_vectors_composite(msv_i_aligned, msv_j_aligned)
#         if score > 0:
#             score_graph.add_edge(i, j, weight=1.-score)
#
#     return nx.minimum_spanning_tree(score_graph)

# def new_multiple_align_ms_vectors(msv_list, mz_tolerance, resolve_by):
#     msv_list = [np.asarray(msv, dtype=float) for msv in msv_list]
#     assert ~np.any([np.any(np.isnan(msv)) for msv in msv_list])
#     assert np.all([msv.shape[0] == 2 for msv in msv_list])
#
#     # Combine msv_list mz-wise and keep track of mz origin
#     msv_combined = np.concatenate(msv_list, axis=1)
#     msv_tracking = np.concatenate([[i]*len(msv_list[i][0])
#                                    for i in range(len(msv_list))])
#
#     sorted_by_mz = msv_combined[0].argsort()
#
#     msv_combined = msv_combined[:, sorted_by_mz]
#
#     print msv_combined.shape
#
#     msv_tracking = msv_tracking[sorted_by_mz]
#
#     if resolve_by == 'distance':
#         vsd = ValueSortedDict()
#     if resolve_by == 'intensity':
#         vsd = ValueSortedDict(lambda v:-v)
#
#     matches = find_all_ms_matches(msv_combined,
#                                   msv_combined,
#                                   mz_tolerance)[0]
#
#     for i, l in enumerate(matches):
#         if len(l) > 1:
#             for j in l[l.index(i) + 1:]:
#                 if resolve_by == 'distance':
#                     vsd[(i, j)] = msv_combined[0, j] - msv_combined[0, i]
#                 if resolve_by == 'intensity':
#                     vsd[(i, j)] = msv_combined[1, i] * msv_combined[1, j]
#
#     to_merge = {i: SortedSet([i]) for i in range(msv_combined.shape[1] + 1)}
#
#     for (i, j) in vsd.iterkeys():
#         if set([msv_tracking[e] for e in to_merge[i]]) & set([msv_tracking[e] for e in to_merge[j]]) == set():
#             if vsd.has_key((to_merge[i][0], to_merge[j][-1])):
#                 merge_set = SortedSet.union(to_merge[i], to_merge[j])
#                 for s in merge_set:
#                     to_merge[s] = merge_set
#
#     for merge_set in set([tuple(s) for s in to_merge.values()]):
#         if len(merge_set) > 1:
#             msv_combined[:, merge_set[0]] = combine_ms_vectors(msv_combined[:, merge_set].T,
#                                                                'mean')
#             msv_combined[:, merge_set[1:]] = np.nan
#
#     msv_combined = msv_combined[:, ~np.isnan(msv_combined[0])]
#     msv_combined = msv_combined[:, msv_combined[0].argsort()]
#
#     print msv_combined.shape
#
#     msv_alignment = np.empty((len(msv_list), msv_combined.shape[0], msv_combined.shape[1]))
#
#     for i, msv in enumerate(msv_list):
#         msv_alignment[i] = pairwise_align_ms_vectors(msv, msv_combined, mz_tolerance, resolve_by)[0]
#
#     if resolve_by == 'distance':
#         vsd = ValueSortedDict()
#     if resolve_by == 'intensity':
#         vsd = ValueSortedDict(lambda v:-v)
#
#     matches = find_all_ms_matches(msv_combined,
#                                   msv_combined,
#                                   mz_tolerance)[0]
#
#     for i, l in enumerate(matches):
#         if len(l) > 1:
#             for j in l[l.index(i) + 1:]:
#                 if resolve_by == 'distance':
#                     vsd[(i, j)] = msv_combined[0, j] - msv_combined[0, i]
#                 if resolve_by == 'intensity':
#                     vsd[(i, j)] = msv_combined[1, i] * msv_combined[1, j]
#
#     to_compress = {i: SortedSet([i]) for i in range(msv_combined.shape[1] + 1)}
#     compressed_mzi = {}
#
#     for (i, j) in vsd.iterkeys():
#         compress_set = SortedSet.union(to_compress[i], to_compress[j])
#         mzi = combine_ms_vectors(np.concatenate(msv_alignment[:, :, compress_set], axis=1).T,
#                                  'mean')
#         print mzi
#         print msv_combined[:, i]
#         print msv_combined[:, j]
#         print
#
#         if np.all(np.sum(~np.isnan(msv_alignment[:, 0, compress_set]), axis=1) <= 1):
#             if (np.all(np.isnan(msv_alignment[:, 0, compress_set])):
#                 for s in compress_set:
#                     to_compress[s] = compress_set
#                     compressed_mzi[s] = mzi
#             elif ((mzi[0] - mz_tolerance <= np.nanmin(msv_alignment[:, 0, compress_set]))
#             and (np.nanmax(msv_alignment[:, 0, compress_set]) <= mzi[0] + mz_tolerance))):
#
#                 for s in compress_set:
#                     to_compress[s] = compress_set
#                     compressed_mzi[s] = mzi
#     n = 0
#     for compress_set in set([tuple(s) for s in to_compress.values()]):
#         if len(compress_set) > 1:
#             n += 1
#             msv_combined[:, compress_set[0]] = compressed_mzi[compress_set[0]]
#             msv_combined[:, compress_set[1:]] = np.nan
#     print n
#
#     msv_combined = msv_combined[:, ~np.isnan(msv_combined[0])]
#     msv_combined = msv_combined[:, msv_combined[0].argsort()]
#
#     msv_alignment = np.empty((len(msv_list), msv_combined.shape[0], msv_combined.shape[1]))
#
#     print msv_combined.shape
#
#     for i, msv in enumerate(msv_list):
#         msv_alignment[i] = pairwise_align_ms_vectors(msv, msv_combined, mz_tolerance, resolve_by)[0]
#
#     return msv_alignment, msv_combined


# def multiple_align_ms_vectors(msv_list, mz_tolerance, resolve_by, combine_by,
#                               mass_power=0, intensity_power=1, weights=None):
#     """EXPERIMENTAL: Pairwise alignments can lead to M/Z-wise misalignments in the multiple alignment.
#
#     Pairwise aligns each combination of pairs of ms vectors in msv_list (see pairwise_align_ms_vectors),
#     scores them (see score_ms_vectors_composite_dot), assigns their scoreto a corresponding position in
#     a score matrix, combines the two ms vectors with highest score (see combine_ms_vectors) according to
#     their weights and weighs it as the sum of its constituent weights, removes the constituent ms vectors and
#     weights from msv_list and weights, inserts the combined ms vector and weight to the
#     front of msv_list and weights, repeats previous steps until msv_list only contains the combination
#     of all ms vectors, and returns a tuple containing a numpy 3d array containing each original
#     ms vector pairwise aligned to the combination of all ms vectors and this combination ms vector.
#
#     The algorithm is analogous to clustal used for sequence alignment by using a distance matrix produced
#     by pairwise ms alignment and a modified neighbor-joining method.
#
#     :param msv_list: list of numpy 2d array, msv_1[0] is m/z and msv_1[1] is intensities sorted by m/z
#     :param mz_tolerance: float, limits matches to be within +/- mz_tolerance
#     :param resolve_by: 'distance', 'shape', or 'intensity',
#         determines what to prioritize when there are multiple matching m/zs within mz_tolerance
#         (see match_ms_vectors)
#
#     :return msv_alignment: numpy 3d array,
#         msv_alignment[i] is ms vector at position i,
#         msv_alignment[i:,0] is  m/z values at position i,
#         msv_alignment[i:,1] is intensity values at position i
#     :return msv_combined: numpy 2d array, total combination of all ms vectors in msv_list
#     """
#
#
#     print 'EXPERIMENTAL: Pairwise alignments can lead to M/Z-wise misalignments in the multiple alignment.'
#
#     msv_list = [np.asarray(msv, dtype=float) for msv in msv_list]
#     assert np.any([np.any(np.isnan(msv)) for msv in msv_list]) == False
#     assert np.all([msv.shape[0] == 2 for msv in msv_list]) == True
#     assert not (resolve_by == 'intensity' and combine_by == 'mean')
#
#     # Make a copy of the original ms vectors
#     msv_original_list = msv_list[:]
#
#     # Set weights of ms vectors uniformly if weights is not set manually
#     if weights is None:
#         weights = [1. / len(msv_list)] * len(msv_list)
#
#     # Initialize score_matrix
#     num_msv = len(msv_list)
#     score_matrix = np.empty((num_msv, num_msv))
#     score_matrix.fill(np.nan)
#
#     # Find the two most similar ms vectors by score, combine them, remove them from msv_list,
#     # add the combination to msv_list, repeat until there is only the combination of all ms vectors
#     while len(msv_list) > 1:
#
#         # Store number of ms vectors
#         num_msv = len(msv_list)
#
#         if np.isnan(score_matrix).all():
#             # Add weighted score of msv_list[i] and msv_list[j] to score_matrix[i,j] for all combinations
#             for i, j in np.array(np.tril_indices(num_msv, -1)).T:
#                 msv_i_aligned, msv_j_aligned = pairwise_align_ms_vectors(msv_list[i], msv_list[j], mz_tolerance, resolve_by)
#                 score_matrix[i, j] = score_ms_vectors_composite(msv_i_aligned, msv_j_aligned, mass_power=mass_power, intensity_power=intensity_power)
#         else:
#             # Add weighted score of msv_list[i] and msv_list[0] to score_matrix[i,0] for new ms vector only
#             for i in range(1, num_msv):
#                 msv_i_aligned, msv_j_aligned = pairwise_align_ms_vectors(msv_list[i], msv_list[0], mz_tolerance, resolve_by)
#                 score_matrix[i, 0] = score_ms_vectors_composite(msv_i_aligned, msv_j_aligned, mass_power=mass_power, intensity_power=intensity_power)
#
#
#         # Flatten the score_matrix
#         flat_score_matrix = score_matrix.ravel()
#
#         # Mask the score_matrix to only see valid scores
#         masked_flat_score_matrix = np.ma.array(flat_score_matrix,
#                                                mask=np.isnan(flat_score_matrix),
#                                                fill_value=-np.inf)
#
#         # Find indices of score_matrix with the maximum score
#         max_score = np.unravel_index(np.ma.argmax(masked_flat_score_matrix),
#                                      (num_msv, num_msv))
#
#         # Pairwise align the ms vectors yielding the best score
#         pair = pairwise_align_ms_vectors(msv_list[max_score[0]],
#                                          msv_list[max_score[1]],
#                                          mz_tolerance,
#                                          resolve_by)
#
#         # Combine the ms vectors yielding the best score into msv_combined
#         msv_combined = combine_ms_vectors(pair, combine_by,
#                                           weights=[weights[max_score[0]],
#                                                    weights[max_score[1]]])
#
#         # Remove the ms vectors yielding the best scores from msv_list and score_matrix
#         msv_list.pop(max_score[0])
#         msv_list.pop(max_score[1])
#         score_matrix = np.delete(score_matrix, max_score[0], 0)
#         score_matrix = np.delete(score_matrix, max_score[0], 1)
#         score_matrix = np.delete(score_matrix, max_score[1], 0)
#         score_matrix = np.delete(score_matrix, max_score[1], 1)
#
#         # Add the msv_combined to msv_list and score_matrix
#         msv_list.insert(0, msv_combined)
#         score_matrix = np.pad(score_matrix, ((1,0),(1,0)), mode='constant', constant_values=np.nan)
#
#         # Remove constituent weights and add combined weight to weights
#         weights.insert(0, weights.pop(max_score[0]) + weights.pop(max_score[1]))
#
#     # Initialize msv_alignment
#     msv_alignment = np.empty((len(msv_original_list), 2, msv_list[0][0].size))
#
#     # Pairwise align all ms vectors from msv_original_list to final msv_combined
#     # and add to msv_alignment
#     print msv_combined
#     for i, msv in enumerate(msv_original_list):
#         msv_alignment[i] = pairwise_align_ms_vectors(msv, msv_combined, mz_tolerance, resolve_by)[0]
#
#     return msv_alignment, msv_combined


def partition_aligned_ms_vectors(msv_1_aligned, msv_2_aligned):
    """
    Finds indices where m/z values in msv_1_aligned and msv_2_aligned are matching (non-nan for both ms vectors)
    and nonmatching (nan for only one ms vector) and returns 4 numpy 2d arrays with msv_1_aligned and msv_2_aligned
    partitioned into matches and nonmatches. Index i of matching partitions match with index i of the
    other and all indices of nonmatching partitions do not match any index of the other.

    :param msv_1_aligned: numpy 2d array, msv_1_aligned[0] is m/z and msv_1_aligned[1] is intensities sorted by m/z,
        each index of msv_1_aligned has its m/z value match that of msv_2_aligned if non-nan
    :param msv_2_aligned: numpy 2d array, msv_2_aligned[0] is m/z and msv_2_aligned[1] is intensities sorted by m/z,
        each index of msv_2_aligned has its m/z value match that of msv_1_aligned if non-nan

    :return msv_1_matches: numpy 2d array,
    :return msv_2_matches: numpy 2d array,
    :return msv_1_nonmatches: numpy 2d array,
    :return msv_2_nonmatches: numpy 2d array,
    """

    msv_1_aligned = np.asarray(msv_1_aligned, dtype=float)
    msv_2_aligned = np.asarray(msv_2_aligned, dtype=float)
    assert msv_1_aligned.shape[0] == 2 and msv_2_aligned.shape[0] == 2
    assert msv_1_aligned.shape == msv_2_aligned.shape

    msv_1_matches = msv_1_aligned[:,~np.isnan(msv_1_aligned[0]) & ~np.isnan(msv_2_aligned[0])]
    msv_2_matches = msv_2_aligned[:,~np.isnan(msv_1_aligned[0]) & ~np.isnan(msv_2_aligned[0])]
    msv_1_nonmatches = msv_1_aligned[:,~np.isnan(msv_1_aligned[0]) & np.isnan(msv_2_aligned[0])]
    msv_2_nonmatches = msv_2_aligned[:, np.isnan(msv_1_aligned[0]) & ~np.isnan(msv_2_aligned[0])]

    return msv_1_matches, msv_2_matches, msv_1_nonmatches, msv_2_nonmatches


def weigh_vector_by_mz_and_intensity(msv, mz_power, intensity_power):
    """
    Returns (mz^mz_power)*(intensity^intensity_power) vector

    :param msv: numpy 2d array, ms_v[0] is m/z values and ms_v[1] is intensity values
    :param mz_power: float, power scales m/z
    :param intensity_power: float, power scales intensity

    :return: numpy 1d array
    """

    weighted_vector = np.multiply(np.power(msv[0], mz_power),
                                  np.power(msv[1], intensity_power))# + intensity_power**((-intensity_power+1)**-1), intensity_power))

    return weighted_vector / np.linalg.norm(weighted_vector)


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

    # Find shared indices between v1 and v2
    shared = (v1 != 0) & ~np.isnan(v1) & (v2 != 0) & ~np.isnan(v2)

    # Select shared indices from v1 and v2
    v1_shared = v1[shared]
    v2_shared = v2[shared]

    if sum(shared) == 1 and len(v1) != 1 and len(v2) != 1:
        return 0

    # Calculate ratio pairs
    ratios = (v1_shared/np.roll(v1_shared, 1)) * (np.roll(v2_shared, 1)/v2_shared)

    # Inverse ratio if greater than one
    ratios = np.minimum(ratios, 1 / ratios)

    return np.nan_to_num(np.sum(ratios) / np.sum(shared))


def score_ms_vectors_dot(msv_sample_aligned, msv_ref_aligned, **kwargs):
    """
    Returns dot product score similar defined in Table 1 of paper below

    :param msv_1_aligned: numpy 2d array, msv_1_aligned[0] is m/z and msv_1_aligned[1] is intensities sorted by m/z,
        each index of msv_1_aligned has its m/z value match that of msv_2_aligned if non-nan
    :param msv_2_aligned: numpy 2d array, msv_2_aligned[0] is m/z and msv_2_aligned[1] is intensities sorted by m/z,
        each index of msv_2_aligned has its m/z value match that of msv_1_aligned if non-nan
    :param kwargs: Ask Daniel or look at the code itself

    :return: float

    Stein, S. E., & Scott, D. R. (1994).
    Optimization and testing of mass spectral library search algorithms for compound identification.
    Journal of the American Society for Mass Spectrometry,
    5(9), 859-866. doi:10.1016/1044-0305(94)87009-8
    """

    msv_sample_aligned = np.asarray(msv_sample_aligned, dtype=float)
    msv_ref_aligned = np.asarray(msv_ref_aligned, dtype=float)

    assert msv_sample_aligned.shape[0] == 2 and msv_ref_aligned.shape[0] == 2
    assert msv_sample_aligned.shape == msv_ref_aligned.shape

    shared = ~np.isnan(msv_sample_aligned[0]) & ~np.isnan(msv_ref_aligned[0])

    mass_power = kwargs.pop('mass_power', 0)
    intensity_power = kwargs.pop('intensity_power', 1 / np.log2(sum(shared) + 1))

    # Weigh vectors according to mass and intensity
    msv_sample_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_sample_aligned),
                                                           mass_power, intensity_power)
    msv_ref_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_ref_aligned),
                                                        mass_power, intensity_power)

    # Compute dot product
    dot = np.dot(msv_sample_weighted, msv_ref_weighted)

    return dot


def score_ms_vectors_composite(msv_sample_aligned, msv_ref_aligned,
                               **kwargs):
    """
    Returns composite score similar to that defined in Table 1 of paper below

    :param msv_1_aligned: numpy 2d array, msv_1_aligned[0] is m/z and msv_1_aligned[1] is intensities sorted by m/z,
        each index of msv_1_aligned has its m/z value match that of msv_2_aligned if non-nan
    :param msv_2_aligned: numpy 2d array, msv_2_aligned[0] is m/z and msv_2_aligned[1] is intensities sorted by m/z,
        each index of msv_2_aligned has its m/z value match that of msv_1_aligned if non-nan
    :param kwargs: Ask Daniel or look at the code

    :return: float

    Stein, S. E., & Scott, D. R. (1994).
    Optimization and testing of mass spectral library search algorithms for compound identification.
    Journal of the American Society for Mass Spectrometry,
    5(9), 859-866. doi:10.1016/1044-0305(94)87009-8
    """

    msv_sample_aligned = np.asarray(msv_sample_aligned, dtype=float)
    msv_ref_aligned = np.asarray(msv_ref_aligned, dtype=float)

    assert msv_sample_aligned.shape[0] == 2 and msv_ref_aligned.shape[0] == 2
    assert msv_sample_aligned.shape == msv_ref_aligned.shape

    shared = ~np.isnan(msv_sample_aligned[0]) & ~np.isnan(msv_ref_aligned[0])

    if sum(shared) == 0:
        return 0

    mean_type = kwargs.pop('mean_type', 'geo')

    mass_power = kwargs.pop('mass_power', 0)
    dot_mass_power = mass_power
    ratio_mass_power = mass_power
    intensity_power = kwargs.pop('intensity_power', 1 / np.log2(sum(shared) + 1))
    dot_intensity_power = intensity_power
    ratio_intensity_power = intensity_power

    dot_mass_power = kwargs.pop('dot_mass_power', 0)
    dot_intensity_power = kwargs.pop('dot_intensity_power', 1 / np.log2(sum(shared) + 1))
    ratio_mass_power = kwargs.pop('ratio_mass_power', 0)
    ratio_intensity_power = kwargs.pop('ratio_intensity_power', 1 / np.log2(sum(shared) + 1))

    # Weigh vectors according to mass and intensity
    msv_sample_dot_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_sample_aligned),
                                                               dot_mass_power, dot_intensity_power)
    msv_ref_dot_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_ref_aligned),
                                                            dot_mass_power, dot_intensity_power)
    msv_sample_ratio_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_sample_aligned),
                                                                 ratio_mass_power, ratio_intensity_power)
    msv_ref_ratio_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_ref_aligned),
                                                              ratio_mass_power, ratio_intensity_power)

    # Compute dot product and ratio of pairs

    dot = np.dot(msv_sample_dot_weighted, msv_ref_dot_weighted)
    ratio = calc_ratio_of_pairs(msv_sample_ratio_weighted, msv_ref_ratio_weighted)

    dot_weight = np.nansum([msv_sample_dot_weighted, msv_ref_dot_weighted])
    ratio_weight = np.sum([msv_sample_ratio_weighted[shared], msv_ref_ratio_weighted[shared]])


    # Compute mean with 1 of 4 methods
    try:
        if mean_type == 'arith':
            comp = ((dot_weight * dot) + (ratio_weight * ratio)) / (dot_weight + ratio_weight)
        if mean_type == 'geo':
            comp = ((dot ** dot_weight) * (ratio ** ratio_weight)) ** (1. / (dot_weight + ratio_weight))
        if mean_type == 'harm':
            comp = (dot_weight + ratio_weight) / ((dot_weight / dot) + (ratio_weight / ratio))
        if mean_type == 'rms':
            comp = (((dot_weight * (dot ** 2.)) + (ratio_weight * (ratio ** 2.))) / (dot_weight + ratio_weight)) ** .5
    except ZeroDivisionError:
        return 0

    return np.nan_to_num(comp)


# def score_ms_vectors_relevance(msv_query_aligned, msv_target_aligned, **kwargs):
#     """
#     Returns score that is a proxy for the probablity that the queried msv
#     is part of the target msv
#
#     :param msv_1_aligned: numpy 2d array, msv_1_aligned[0] is m/z and msv_1_aligned[1] is intensities sorted by m/z,
#         each index of msv_1_aligned has its m/z value match that of msv_2_aligned if non-nan
#     :param msv_2_aligned: numpy 2d array, msv_2_aligned[0] is m/z and msv_2_aligned[1] is intensities sorted by m/z,
#         each index of msv_2_aligned has its m/z value match that of msv_1_aligned if non-nan
#
#     :return: float
#     """
#
#     msv_query_aligned = np.asarray(msv_query_aligned, dtype=float)
#     msv_target_aligned = np.asarray(msv_target_aligned, dtype=float)
#
#     assert msv_query_aligned.shape[0] == 2 and msv_target_aligned.shape[0] == 2
#     assert msv_query_aligned.shape == msv_target_aligned.shape
#
#     shared = ~np.isnan(msv_query_aligned[0]) & ~np.isnan(msv_target_aligned[0])
#
#     mass_power = kwargs.pop('mass_power', 0)
#     intensity_power = kwargs.pop('intensity_power', 1 / np.log2(sum(shared) + 1))
#
#     # Weigh vectors according to mass and intensity
#     msv_query_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_query_aligned),
#                                                           mass_power, intensity_power)
#     msv_target_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_target_aligned),
#                                                            mass_power, intensity_power)
#
#     ratio = calc_ratio_of_pairs(msv_query_weighted, msv_target_weighted)
#
#     if ratio == 0:
#         return 0
#
#     relevance = ratio * ((np.sum(msv_query_aligned[1, shared]) / np.nansum(msv_query_aligned[1])))
#
#     return relevance


################################################################################
#Spectral search related functions
################################################################################
def gaussian_kernel(mzs1, mzs2, sigma=.01):
    return np.exp(-0.5 * np.subtract.outer(mzs1, mzs2)**2 / sigma**2)

def gaussian_dot_product(msv1, msv2, sigma=.01):
    mzs1,intensities1 = msv1
    mzs2,intensities2 = msv2
    K = gaussian_kernel(mzs1, mzs2, sigma=.01)

    return np.sqrt(intensities1).dot(K).dot(np.sqrt(intensities2))

def shifted_gaussian_dot_product(msv1, msv2, pmz1, pmz2, sigma=.01):
    mzs1,intensities1 = msv1
    mzs2,intensities2 = msv2

    K = gaussian_kernel(pmz1-mzs1, pmz2-mzs2, sigma=.01)

    return np.sqrt(intensities1).dot(K).dot(np.sqrt(intensities2))

def hybrid_gaussian_dot_product(msv1, msv2, pmz1, pmz2, sigma=.01):
    mzs1,intensities1 = msv1
    mzs2,intensities2 = msv2

    K =  .5*gaussian_kernel(mzs1, mzs2, sigma=.01)
    K += .5*gaussian_kernel(pmz1-mzs1, pmz2-mzs2, sigma=.01)

    return np.sqrt(intensities1).dot(K).dot(np.sqrt(intensities2))

def normalized_gaussian_dot_product(msv1, msv2, sigma=.01):
    return gaussian_dot_product(msv1,msv2,sigma)/np.sqrt(gaussian_dot_product(msv1,msv1,sigma)*gaussian_dot_product(msv2,msv2,sigma))

def normalized_shifted_gaussian_dot_product(msv1, msv2, pmz1, pmz2, sigma=.01):
    return shifted_gaussian_dot_product(msv1, msv2, pmz1, pmz2, sigma)/np.sqrt(gaussian_dot_product(msv1,msv1,sigma)*gaussian_dot_product(msv2,msv2,sigma))

def normalized_hybrid_gaussian_dot_product(msv1, msv2, pmz1, pmz2, sigma=.01):
    return hybrid_gaussian_dot_product(msv1, msv2, pmz1, pmz2, sigma)/np.sqrt(gaussian_dot_product(msv1,msv1,sigma)*gaussian_dot_product(msv2,msv2,sigma))


def search_ms_refs(msv_query, **kwargs):
    """
    Returns dataframe with search results and msv_query vs msv_ref score
    from ref_df which matches query.

    :param msv_query: numpy 2d array, msv_query[0] is m/z and msv_query[1] is intensities sorted by m/z
    :param kwargs: Ask Daniel or look at the code

    Common query:
    query='(@polarity == polarity) and ((@precursor_mz - @pre_mz_tolerance) <= precursor_mz <= (@precursor_mz + @pre_mz_tolerance))'

    :return: results_df
    """
    kwargs = dict(locals(), **kwargs)

    # Alignment parameters
    resolve_by = kwargs.pop('resolve_by', 'shape')
    frag_mz_tolerance = kwargs.pop('frag_mz_tolerance', .005)

    # Reference parameters
    ref_loc = kwargs.pop('ref_loc', '/global/project/projectdirs/metatlas/projects/spectral_libraries/msms_refs.tab')
    ref_dtypes = kwargs.pop('ref_dtypes', {'database':str, 'id':str, 'name':str,
                                           'spectrum':object,'decimal':float, 'precursor_mz':float,
                                           'polarity':str, 'adduct':str, 'fragmentation_method':str,
                                           'collision_energy':str, 'instrument':str, 'instrument_type':str,
                                           'formula':str, 'exact_mass':float,
                                           'inchi_key':str, 'inchi':str, 'smiles':str})
    ref_index = kwargs.pop('ref_index', ['database', 'id'])
    query = kwargs.pop('query', 'index == index or index == @pd.NaT')
    post_query = kwargs.pop('post_query', 'index == index or index == @pd.NaT')

    if 'do_gdp' in kwargs:
        do_gdp = kwargs.pop('do_gdp')
    else:
        do_gdp = False

    if 'ref_df' in kwargs:
        ref_df = kwargs.pop('ref_df')
    else:
        ref_df = pd.read_csv(ref_loc,
                             sep='\t',
                             dtype=ref_dtypes
                            ).set_index(ref_index)

    ref_df = ref_df.query(query, local_dict=dict(locals(), **kwargs))

    #setup MatchMS settings
    cos = CosineHungarian(tolerance=frag_mz_tolerance)
    mms_query = Spectrum(mz=msv_query[0], intensities=msv_query[1], metadata={'precursor_mz':np.nan})

    if ref_df['spectrum'].apply(type).eq(str).all():
        ref_df.loc[:, 'spectrum'] = ref_df['spectrum'].apply(lambda s: eval(s)).apply(np.array)
        for _, row in ref_df.iterrows():
            row.loc[:, 'spectrum'] = row['spectrum'][:, row['spectrum'][0] < row['precursor_mz'] + 2.5]

    # Define function to score msv_qury against msv's in reference dataframe
    def score_and_num_matches(msv_ref):
        if do_gdp==True:
#             print(msv_query)
            score = normalized_gaussian_dot_product(msv_query,msv_ref,sigma=frag_mz_tolerance)
            num_matches = gaussian_kernel(msv_query[0],msv_ref[0],sigma=frag_mz_tolerance).sum()
            msv_query_aligned = msv_query
            msv_ref_aligned = msv_ref
        else:
            mms_ref = Spectrum(mz=msv_ref[0], intensities=msv_ref[1], metadata={'precursor_mz':np.nan})
            mms_comparison = cos.pair(mms_query, mms_ref)

            score = mms_comparison['score'].item()
            num_matches = mms_comparison['matches'].item()

            msv_query_aligned, msv_ref_aligned = pairwise_align_ms_vectors(msv_query, msv_ref, frag_mz_tolerance)
        # sensitivity = score_ms_vectors_relevance(msv_ref_aligned, msv_query_aligned,
        #                                          **kwargs)
        # precision = score_ms_vectors_relevance(msv_query_aligned, msv_ref_aligned,
        #                                        **kwargs)

        return pd.Series([score, num_matches, # sensitivity, precision,
                          msv_query_aligned, msv_ref_aligned],
                         index=['score', 'num_matches', # 'sensitivity', 'precision',
                                'msv_query_aligned', 'msv_ref_aligned'])

    score_df = ref_df['spectrum'].apply(score_and_num_matches)

    if len(score_df) > 0:
        score_df = score_df.query(post_query)

    return score_df

################################################################################
#Isotope Distribution Related Functions
################################################################################

def parse_formula(formula, elements=[]):
    """
    Takes formula and optional list of elements and returns list of elements found in formula
    (in order of input elements if set) and list of their counts.

    :param formula: string,
    :param elements: list of strings

    :return: list of strings,
    :return: list of ints,
    """

    in_paren = re.compile("\(([A-Z0-9]+?)\)(\d*)")

    expanded = formula

    #Removing parentheses by expanding what they contain until there are none
    while re.search(in_paren, expanded) != None:
        for match in re.finditer(in_paren, expanded):
            expanded = re.sub(in_paren, int(match.group(2)) * match.group(1), expanded, 1)

    #Count elements in expanded formula
    element_count = Counter()
    for e, c in [(e, int(c if c != '' else 1)) for e, c in
                 re.findall(r'([A-Z][a-z]*)(\d*)', expanded)]:
        element_count[e] += c

    if not elements:
        elements = list(element_count.keys())

    return list(zip(*[(e, element_count[e]) for e in elements]))


def make_isotope_matrix(elements,
                        isotope_file=os.path.join(os.path.dirname(__file__), 'isotope_dict.json'),
                        scale_factor=100000):
    """
    Takes list of elements, and optional isotope dictionary and scaling factor.
    Returns an isotope matrix and vector containing masses removed from it to be used in make_isotope_distribution.

    :param elements: list of strings, elements must be in keys of isotope_dict
    :param isotope_dict: dictionary with elements as keys and isotope information as values
    :param scale_factor: integer, scales masses up in isotope_matrix so FFT in make_isotope_distribution
        is faster

    :return isotope_matrix: numpy 3d matrix of shape (max_elements, max_isotopes, 2),
    :return mass_removed_vec: list of ints,
    """

    import json

    with open(isotope_file) as json_file:
        isotope_dict = json.load(json_file)

    max_elements = len(elements) + 1
    max_isotopes = 0

    for e in elements:
        if isotope_dict[e]['num_isotopes'] > max_isotopes:
            max_isotopes = isotope_dict[e]['num_isotopes']

    isotope_matrix = np.zeros((max_elements, max_isotopes, 2))

    isotope_matrix[-1, 0] = [scale_factor, 1]
    mass_removed_vec = []

    for i, e in enumerate(elements):
        mass_removed_vec.append(int(isotope_dict[e]['isotopes'][0]['mass_number'] - 1))
        for j, iso in enumerate(isotope_dict[e]['isotopes']):
            isotope_matrix[i, j] = [scale_factor * (iso['atomic_mass'] - int(isotope_dict[e]['isotopes'][0]['mass_number'] - 1)), iso['abundance']]

    mass_removed_vec.append(-1)

    return np.nan_to_num(isotope_matrix), mass_removed_vec


def add_synthetic_element(synthetic,
                          isotope_file=os.path.join(os.path.dirname(__file__), 'isotope_dict.json')):
    """
    Takes a dictionary with synthetic element labels (like 'D' for deuterium or 'H' if you wish to override H) as keys
    and tuples as values, and optional isotope dictionary. Tuple[1] is a list of floats which sum to 1 and override the
    abundances of the isotopes of the element given by Tuple[0] and isotope_dict.
    Returns an new isotope dictionary with the synthetic element(s) added to it.

    :param synthetic: dictionary of form {synth:(ele,[a,b,...]), ...},
        synth is synthetic element label,
        ele must be in keys of isotope_dict,
        a,b,... are floats that sum to 1,
    :param isotope_dict: dictionary with elements as keys and isotope information as values,

    :return new_isotope_dict: dictionary with elements as keys and isotope information as values,
    """

    import json

    with open(isotope_file) as json_file:
        isotope_dict = json.load(json_file)

    # Make a deep copy of isotope_dict to modify
    new_isotope_dict = deepcopy(isotope_dict)

    # Modify the deep copy
    for e in synthetic.keys():
        new_isotope_dict[e] = deepcopy(isotope_dict[synthetic[e][0]])
        for i in range(new_isotope_dict[e]['num_isotopes']):
            new_isotope_dict[e]['isotopes'][i]['abundance'] = synthetic[e][1][i]

    return new_isotope_dict


def make_isotope_distribution(element_vector, isotope_matrix, mass_removed_vec,
                              scale_factor=100000,
                              cutoff=1e-4, window_size=100, resolution=0.001):
    '''
    %
    % Calculates isotopic distributions including isotopic fine structure
    % of molecules using FFT and various scaling 'tricks'. Easily adopted
    % to molecules of any elemental composition (by altering max_elements
    % and the nuclide matrix A). To simulate spectra, convolute with peak
    % shape using FFT.
    %
    % (C) 1999 by Magnus Palmblad, Division of Ion Physics, Uppsala Univ.
    % Acknowledgements:
    % Lars Larsson-Cohn, Dept. of Mathematical Statistics, Uppsala Univ.,
    % for help on theory of convolutions and FFT.
    % Jan Axelsson, Div. of Ion Physics, Uppsala Univ. for comments and ideas
    %
    % Contact Magnus Palmblad at magnus.palmblad@gmail.com if you should
    % have any questions or comments.
    %
    Converted to Python 1/10/08 by
    Brian H. Clowers bhclowers@gmail.com
    October 31, 2014
    Added Phosphorous and chemical formula parsing
    Added conditional specification of stable isotope composition
    Ben Bowen, ben.bowen@gmail.com
    fVec is a vector representing the chemical formula including deuterium
    # [H, C, N, O, S, P, D]
    DAbund is the amount of deuterium [0-1], 0.05 is typical
    September 25, 2017
    Added arbitrary elemental composition capability
    Added ability to track which nuclides contribute to which peaks
    Daniel Treen, dgct@jps.net
    '''

    def next2pow(x):
        return 2 ** int(np.ceil(np.log(float(x)) / np.log(2.0)))

    element_vector = np.append(element_vector, 0)
    max_elements = isotope_matrix[:, 0, 0].size  # add 1 due to mass correction 'element'
    max_isotopes = isotope_matrix[0, :, 0].size  # maxiumum # of isotopes for one element

    if resolution < 0.00001:  # % minimal mass step allowed
        resolution = 0.00001
    elif resolution > 0.5:  # maximal mass step allowed
        resolution = 0.5

    R = 0.00001 / resolution  # % R is used to scale nuclide masses (see below)

    window_size = window_size / resolution  # convert window size to new mass units
    window_size = next2pow(window_size)  # fast radix-2 fast-Fourier transform algorithm

    if window_size < np.round(496708 * R) + 1:
        window_size = nextpow2(np.round(496708 * R) + 1)  # just to make sure window is big enough

    M = np.array(element_vector)

    monoisotopic = 0.0
    for i, e in enumerate(element_vector):
        monoisotopic = monoisotopic + ((mass_removed_vec[i] * scale_factor + isotope_matrix[i, 0, 0]) * e / scale_factor)

    Mmi = np.append(np.round(isotope_matrix[:-1, 0, 0] * R), 0) * M  # % (Virtual) monoisotopic mass in new units
    Mmi = Mmi.sum()

    folded = np.floor(Mmi / (window_size - 1)) + 1  # % folded folded times (always one folding due to shift below)

    M[max_elements - 1] = np.ceil(((window_size - 1) - np.mod(Mmi, window_size - 1) + np.round(100000 * R)) * resolution)
    mass_removed = np.array(mass_removed_vec) * M  # ';  % correction for 'virtual' elements and mass shift
    mass_removed = mass_removed.sum()

    start = (folded * (window_size - 1) + 1) * resolution + mass_removed
    stop = (folded + 1) * (window_size - 1) * resolution + mass_removed
    step = window_size - 1

    MA = np.linspace(start, stop, step)

    itAs = np.zeros((max_elements, max_isotopes, window_size), dtype=np.complex_)
    tAs = np.zeros((max_elements, window_size), dtype=np.complex_)

    ptA = np.ones(window_size)

    for i in range(max_elements):
        for j in range(max_isotopes):
            if isotope_matrix[i, j, 0] != 0:
                itA = np.zeros(window_size)
                idx = np.round(isotope_matrix[i, j, 0] * R).astype(int)

                itA[idx] = isotope_matrix[i, j, 1]
                itAs[i, j] = fft(itA)

                tAs[i] += itAs[i, j]

    for i in range(max_elements):
        ptA = ptA * (tAs[i] ** M[i])

    fft_ptA = ptA.copy()
    ptA = ifft(ptA).real

    cutoff_idx = np.where(ptA > cutoff)[0]

    contributions = np.zeros((max_elements - 1, max_isotopes, cutoff_idx.size), dtype=int)

    for i in range(max_elements - 1):
        for j in range(max_isotopes):
            c = ~np.isclose(ptA[cutoff_idx], ifft((fft_ptA / (tAs[i])) * ((tAs[i] - itAs[i, j]))).real[cutoff_idx])
            contributions[i, j, c] = np.round((ifft((fft_ptA / (tAs[i])) * (itAs[i, j])).real[cutoff_idx][c] / ptA[cutoff_idx][c]) * element_vector[i])

    return np.array([MA[cutoff_idx], ptA[cutoff_idx]]), contributions


def sort_msms_hits(msms_hits: pd.DataFrame, sorting_method: str) -> tuple[pd.DataFrame, str]:
    """
    Takes an msms hits dataframe and returns a sorted version of it based on the sorting method. Typically
    this function is called while iterating over compounds, so each dataframe input will likely be for a single compound.

    Note: Every time you add a possible scoring method to this function, update the validate function in the class Atlas in
    the script metatlas.tools.config 
    """
    if sorting_method is None:
        sorting_method = 'cosine_score'
    valid_scoring_method_list = ['cosine_score', 'sums', 'weighted', 'numeric_hierarchy', 'quantile_hierarchy']
    if sorting_method not in valid_scoring_method_list:
        raise ValueError(f"Invalid MSMS sorting method: {sorting_method}. Should be one of {valid_scoring_method_list}.")

    if sorting_method == "cosine_score":
        sorted_msms_hits = msms_hits.sort_values('score', ascending=False)

    elif sorting_method == "sums":
        sorted_msms_hits = msms_hits.copy()
        sorted_msms_hits.loc[:, 'summed_ratios_and_score'] = (
                                                              (sorted_msms_hits['num_matches'] / sorted_msms_hits['ref_frags']) +
                                                              (sorted_msms_hits['num_matches'] / sorted_msms_hits['data_frags']) +
                                                              (sorted_msms_hits['score'])
                                                              )
        sorted_msms_hits = sorted_msms_hits.sort_values('summed_ratios_and_score', ascending=False)

    elif sorting_method == "weighted":
        sorted_msms_hits = msms_hits.copy()
        weights = {'score': 0.5,
                    'match_to_data_frag_ratio': 0.25,
                    'match_to_ref_frag_ratio': 0.25
                    }
        sorted_msms_hits.loc[:, 'weighted_score'] = (
                                                    (sorted_msms_hits['score'] * weights['score']) +
                                                    ((sorted_msms_hits['num_matches'] / sorted_msms_hits['data_frags']) * weights['match_to_data_frag_ratio']) +
                                                    ((sorted_msms_hits['num_matches'] / sorted_msms_hits['ref_frags']) * weights['match_to_ref_frag_ratio'])
                                                    )
        sorted_msms_hits = sorted_msms_hits.sort_values('weighted_score', ascending=False)

    elif sorting_method == "numeric_hierarchy":
        sorted_msms_hits = msms_hits.copy()
        bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99, 1]
        labels = ['0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-0.95', '0.95-0.97', '0.97-0.99', '0.99-1']
        sorted_msms_hits = sorted_msms_hits.dropna(subset=['score'])
        sorted_msms_hits['score_bin'] = pd.cut(sorted_msms_hits['score'], bins=bins, labels=labels, right=False)
        sorted_msms_hits = sorted_msms_hits.sort_values(by=['score_bin', 'num_matches', 'score'], ascending=[False, False, False])

    elif sorting_method == "quantile_hierarchy":
        sorted_msms_hits = msms_hits.copy()
        sorted_msms_hits = sorted_msms_hits.dropna(subset=['score'])
        sorted_msms_hits['score_temp'] = sorted_msms_hits['score'] + np.random.normal(0, 1e-8, size=len(sorted_msms_hits))  # Add small noise to handle duplicates
        sorted_msms_hits['score_bin'] = pd.qcut(sorted_msms_hits['score_temp'], duplicates='drop', q=5)
        sorted_msms_hits = sorted_msms_hits.sort_values(by=['score_bin', 'num_matches', 'score_temp'], ascending=[False, False, False])
        sorted_msms_hits = sorted_msms_hits.drop(columns=['score_temp'])

    return sorted_msms_hits, sorting_method