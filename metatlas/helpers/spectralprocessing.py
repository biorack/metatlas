import re
from copy import deepcopy
from collections import Counter
import os
import pickle
import numpy as np
import numpy.fft.fftpack as F


################################################################################
#Misc Functions
################################################################################

def remove_ms_vector_noise(msv, threshold=1e-5):
    """
    Returns ms vector with intensity below threshold removed

    :param msv: numpy 2d array, msv[0] is m/z values and msv[1] is intensity values
    :param threshold: float

    :return: numpy 2d array
    """
    
    return msv[:, msv[1, :] > threshold]


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

    if 'data' in metatlas_dataset[file_idx][compound_idx]['data']['msms'].keys() and \
                 metatlas_dataset[file_idx][compound_idx]['data']['msms']['data']['mz'].size > 0:

        condition_list = list(filter(None, re.split(r"(\w+|[()])", condition)))

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
    assert (np.any(np.isnan(msv_1)) and np.any(np.isnan(msv_2))) == False
    assert msv_1.shape[0] == 2 and msv_2.shape[0] == 2

    #start_idx is element for which inserting msv_1 directly ahead of it maintains sort order
    start_idxs = np.searchsorted(msv_2[0], msv_1[0] - mz_tolerance)

    #end_idx is element for which inserting msv_1 directly after it maintains sort order
    #found by searching negative reversed list since np.searchshorted requires increasing order
    end_idxs = len(msv_2[0]) - np.searchsorted(-msv_2[0][::-1], -(msv_1[0] + mz_tolerance))

    #if the start and end idx is the same, the peak is too far away in mass from the data and will be empty
    matches = [range(start, end) for start, end in zip(start_idxs, end_idxs)]

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
    assert (np.any(np.isnan(msv_1)) and np.any(np.isnan(msv_2))) == False
    assert msv_1.shape[0] == 2 and msv_2.shape[0] == 2
    
    assert resolve_by == 'distance' or resolve_by == 'shape' or resolve_by == 'intensity'
    
    matches = find_all_ms_matches(msv_1, msv_2, mz_tolerance)[0]
    
    # Create match matrix where row i and column j is value determined by resolve_by if msv_1[i] matches msv_2[j]
    # and a default value if otherwise
    if resolve_by == 'distance':
        match_matrix = np.fromfunction(np.vectorize(lambda i, j: np.absolute(msv_1[0, int(i)] - msv_2[0, int(j)]) if int(j) in matches[int(i)] else np.inf), (msv_1[0].size, msv_2[0].size)).astype(float)
    if resolve_by == 'shape':
        ms_scale = np.median([msv_1[1,i]/msv_2[1,j] for i,l in enumerate(matches) for j in l])
        match_matrix = np.fromfunction(np.vectorize(lambda i, j: np.absolute(msv_1[1, int(i)] - (ms_scale*msv_2[1, int(j)])) if int(j) in matches[int(i)] else np.inf), (msv_1[0].size, msv_2[0].size)).astype(float)
    if resolve_by == 'intensity':
        match_matrix = np.fromfunction(np.vectorize(lambda i, j: msv_1[1, int(i)] * msv_2[1, int(j)] if int(j) in matches[int(i)] else -np.inf), (msv_1[0].size, msv_2[0].size)).astype(float)

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


# def align_vectors(msv_1_matches, msv_2_matches, msv_1_nonmatches,
#                   msv_2_nonmatches):
#         (msv_1_matches, msv_1_nonmatches, np.zeros(msv_2_nonmatches.shape)),
#         axis=1)
#     ms_v2_aligned = np.concatenate(
#         (msv_2_matches, np.zeros(msv_1_nonmatches.shape), msv_2_nonmatches),
#         axis=1)
#
#     return ms_v1_aligned, ms_v2_aligned


def pairwise_align_ms_vectors(msv_1, msv_2, mz_tolerance, resolve_by):
    """
    Finds which m/z values in msv_1 and msv_2 'best' match within +/- mz_tolerance
    (see match_ms_vectors) and returns a numpy 3d array containing aligned versions of msv_1
    and msv_2 meaning index i of one matches index i of the other and nonmatches are filled 
    with nan.
    
    :param msv_1: numpy 2d array, msv_1[0] is m/z and msv_1[1] is intensities sorted by m/z
    :param msv_2: numpy 2d array, msv_2[0] is m/z and msv_2[1] is intensities sorted by m/z
    :param mz_tolerance: float, limits matches to be within +/- mz_tolerance
    :param resolve_by: 'distance', 'shape', or 'intensity',
        determines what to prioritize when there are multiple matching m/zs within mz_tolerance
        (see match_ms_vectors)

    :return msv_alignment: numpy 3d array,
        msv_alignment[i] is ms vector at position i,
        msv_alignment[i:,0] is  m/z values at position i,
        msv_alignment[i:,1] is intensity values at position i
    """
    
    msv_1 = np.asarray(msv_1, dtype=float)
    msv_2 = np.asarray(msv_2, dtype=float)
    assert (np.any(np.isnan(msv_1)) and np.any(np.isnan(msv_2))) == False
    assert msv_1.shape[0] == 2 and msv_2.shape[0] == 2

    # Get matches and nonmatches between msv_1 and msv_2
    msv_1_matches, msv_2_matches, msv_1_nonmatches, msv_2_nonmatches = partition_ms_vectors(msv_1, msv_2, mz_tolerance, resolve_by)
    
    # Initialize msv_1_aligned and msv_2_aligned
    msv_1_aligned = np.asarray(msv_1_matches, dtype=float)
    msv_2_aligned = np.asarray(msv_2_matches, dtype=float)
    
    # Insert msv_1_nonmatches into msv_1_aligned and nan shaped like msv_1_nonmatches in msv_2_aligned 
    msv_1_to_insert = np.searchsorted(msv_1_aligned[0], msv_1_nonmatches[0])
    msv_1_aligned = np.insert(msv_1_aligned, msv_1_to_insert, msv_1_nonmatches, 1)
    msv_2_aligned = np.insert(msv_2_aligned, msv_1_to_insert, np.full_like(msv_1_nonmatches, np.nan), 1)

    # Insert msv_2_nonmatches into msv_2_aligned and nan shaped like msv_2_nonmatches in msv_1_aligned 
    msv_2_to_insert = np.searchsorted(msv_2_aligned[0], msv_2_nonmatches[0])
    msv_2_aligned = np.insert(msv_2_aligned, msv_2_to_insert, msv_2_nonmatches, 1)
    msv_1_aligned = np.insert(msv_1_aligned, msv_2_to_insert, np.full_like(msv_2_nonmatches, np.nan), 1)
    
    msv_alignment = np.dstack((np.array([msv_1_aligned.T, msv_2_aligned.T]))).T

    return msv_alignment


def multiple_align_ms_vectors(msv_list, mz_tolerance, resolve_by, combine_by,
                              mass_power=0, intensity_power=1, weights=None):
    """DO NOT USE. DOES NOT WORK AS WHEN THERE ARE CROSSING MATCHES.
    
    Pairwise aligns each combination of pairs of ms vectors in msv_list (see pairwise_align_ms_vectors), 
    scores them (see score_ms_vectors_composite_dot), assigns the product of their score with 
    the sum of their weights to a corresponding position in a score matrix, combines the 
    two ms vectors with highest score (see combine_ms_vectors) according to their weights 
    and weighs it as the sum of its constituent weights, removes the constituent ms vectors and 
    weights from msv_list and weights, inserts the combined ms vector and weight to the 
    front of msv_list and weights, repeats previous steps until msv_list only contains the combination
    of all ms vectors, and returns a tuple containing a numpy 3d array containing each original
    ms vector pairwise aligned to the combination of all ms vectors and this combination ms vector.
    
    The algorithm is analogous to clustal used for sequence alignment by using a distance matrix produced 
    by pairwise ms alignment and a modified neighbor-joining method.
    
    :param msv_list: list of numpy 2d array, msv_1[0] is m/z and msv_1[1] is intensities sorted by m/z
    :param mz_tolerance: float, limits matches to be within +/- mz_tolerance
    :param resolve_by: 'distance', 'shape', or 'intensity',
        determines what to prioritize when there are multiple matching m/zs within mz_tolerance
        (see match_ms_vectors)

    :return msv_alignment: numpy 3d array,
        msv_alignment[i] is ms vector at position i,
        msv_alignment[i:,0] is  m/z values at position i,
        msv_alignment[i:,1] is intensity values at position i
    :return msv_combined: numpy 2d array, total combination of all ms vectors in msv_list
    """
    
    
    try:
        assert False
    except:
        print 'Do not use. Read the doc for details.'
        raise
    
    msv_list = [np.asarray(msv, dtype=float) for msv in msv_list]
    assert np.any([np.any(np.isnan(msv)) for msv in msv_list]) == False
    assert np.all([msv.shape[0] == 2 for msv in msv_list]) == True
    assert not (resolve_by == 'intensity' and combine_by == 'mean')
    
    # Make a copy of the original ms vectors
    msv_original_list = msv_list[:]

    # Set weights of ms vectors uniformly if weights is not set manually
    if weights == None:
        weights = [1. / len(msv_list)] * len(msv_list)

    # Find the two most similar ms vectors by score, combine them, remove them from msv_list, 
    # add the combination to msv_list, repeat until there is only the combination of all ms vectors
    while len(msv_list) > 1:
        
        # Initialize score_matrix
        num_msv = len(msv_list)
        score_matrix = np.empty((num_msv, num_msv))
        score_matrix.fill(np.nan)

        # Add weighted score of msv_list[i] and msv_list[j] to score_matrix[i,j] for all combinations
        for i, j in np.array(np.tril_indices(num_msv, -1)).T:
            msv_i_aligned, msv_j_aligned = pairwise_align_ms_vectors(msv_list[i], msv_list[j], mz_tolerance, resolve_by)
            score_matrix[i, j] = (weights[i] + weights[j]) * score_ms_vectors_composite_dot(msv_i_aligned, msv_j_aligned, mass_power, intensity_power)

        # Flatten the score_matrix
        flat_score_matrix = score_matrix.ravel()

        # Mask the score_matrix to only see valid scores
        masked_flat_score_matrix = np.ma.array(flat_score_matrix,
                                               mask=np.isnan(flat_score_matrix),
                                               fill_value=-np.inf)
        
        # Find indices of score_matrix with the maximum score
        max_score = np.unravel_index(np.ma.argmax(masked_flat_score_matrix),
                                     (num_msv, num_msv))

        # Pairwise align the ms vectors yielding the best score
        pair = pairwise_align_ms_vectors(msv_list[max_score[0]],
                                         msv_list[max_score[1]],
                                         mz_tolerance,
                                         resolve_by)

        # Combine the ms vectors yielding the best score into msv_combined
        msv_combined = combine_ms_vectors(pair, combine_by,
                                          weights=[weights[max_score[0]],
                                                   weights[max_score[1]]])

        # Remove the ms vectors yielding the best scores from msv_list
        msv_list.pop(max_score[0])
        msv_list.pop(max_score[1])

        # Add the msv_combined to msv_list
        msv_list.insert(0, msv_combined)
        # Remove constituent weights and add combined weight to weights
        weights.insert(0, weights.pop(max_score[0]) + weights.pop(max_score[1]))

    # Initialize msv_alignment
    msv_alignment = np.empty((len(msv_original_list), 2, msv_list[0][0].size))

    # Pairwise align all ms vectors from msv_original_list to final msv_combined
    # and add to msv_alignment
    for i, msv in enumerate(msv_original_list):
        msv_alignment[i] =  pairwise_align_ms_vectors(msv, msv_combined, mz_tolerance, resolve_by)[0]

    return msv_alignment, msv_combined


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


def weigh_vector_by_mz_and_intensity(msv, mz_power=1, intensity_power=.6):
    """
    Returns (mz^mz_power)*(intensity^intensity_power) vector

    :param msv: numpy 2d array, ms_v[0] is m/z values and ms_v[1] is intensity values
    :param mz_power: float, power scales m/z
    :param intensity_power: float, power scales intensity

    :return: numpy 1d array
    """

    return np.multiply(np.power(msv[0], mz_power),
                       np.power(msv[1], intensity_power))


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
    shared = ~np.isnan(v1) & ~np.isnan(v2)

    # Select shared indices from v1 and v2
    v1_shared = v1[shared]
    v2_shared = v2[shared]

    # Calculate ratio pairs
    ratios = (v1_shared/np.roll(v1_shared, 1)) * (np.roll(v2_shared, 1)/v2_shared)

    # Inverse ratio if greater than one
    ratios = np.minimum(ratios, 1 / ratios)

    return np.nan_to_num(np.sum(ratios) / np.sum(shared))


def score_vectors_dot(v1, v2, normalize=True):
    """
    Returns dot product of two vectors, either normalized to a magnitude of 1 or not.

    :param v1: numpy 1d array
    :param v2: numpy 1d array
    :param normalize: bool, scales to unit vector if true

    :return:
    """

    if normalize:
        return np.dot(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
    else:
        return np.dot(v1, v2)


def score_ms_vectors_composite_dot(msv_sample_aligned, msv_ref_aligned, mass_power=1,
                                intensity_power=.6):
    """
    Returns "Composite" dot product score as defined in Table 1 of paper below

    :param msv_1_aligned: numpy 2d array, msv_1_aligned[0] is m/z and msv_1_aligned[1] is intensities sorted by m/z,
        each index of msv_1_aligned has its m/z value match that of msv_2_aligned if non-nan
    :param msv_2_aligned: numpy 2d array, msv_2_aligned[0] is m/z and msv_2_aligned[1] is intensities sorted by m/z,
        each index of msv_2_aligned has its m/z value match that of msv_1_aligned if non-nan

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
    
    # Find number of peaks shared between sample and reference and in sample
    num_shared = np.sum(~np.isnan(msv_sample_aligned) & ~np.isnan(msv_ref_aligned))
    num_sample = msv_sample_aligned[0,~np.isnan(msv_sample_aligned[0])].size

    # Weigh vectors according to mass and intensity
    msv_sample_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_sample_aligned),
                                                           mass_power,
                                                           intensity_power)
    msv_ref_weighted = weigh_vector_by_mz_and_intensity(np.nan_to_num(msv_ref_aligned),
                                                        mass_power,
                                                        intensity_power)

    # Compute dot product and ratio of pairs
    dot = score_vectors_dot(msv_sample_weighted, msv_ref_weighted)
    ratio = calc_ratio_of_pairs(msv_sample_aligned[1], msv_ref_aligned[1])

    return ((num_sample * dot) + (num_shared * ratio)) / (num_sample + num_shared)

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
        elements = element_count.keys()

    return zip(*[(e, element_count[e]) for e in elements])


def make_isotope_matrix(elements, 
                        isotope_dict=pickle.load(open(os.path.join(os.path.dirname(__file__), 'isotope_dict.pkl'), 'rb')),
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
                          isotope_dict=pickle.load(open(os.path.join(os.path.dirname(__file__), 'isotope_dict.pkl'), 'rb'))):
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
    Added ability to track which isotopes of which elements contribute to which peaks
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

    for i in xrange(max_elements):
        for j in xrange(max_isotopes):
            if isotope_matrix[i, j, 0] != 0:
                itA = np.zeros(window_size)
                idx = np.round(isotope_matrix[i, j, 0] * R).astype(int)

                itA[idx] = isotope_matrix[i, j, 1]
                itAs[i, j] = F.fft(itA)

                tAs[i] += itAs[i, j]

    for i in xrange(max_elements):
        ptA = ptA * (tAs[i] ** M[i])

    fft_ptA = ptA.copy()
    ptA = F.ifft(ptA).real

    cutoff_idx = np.where(ptA > cutoff)[0]

    contributions = np.zeros((max_elements - 1, max_isotopes, cutoff_idx.size), dtype=int)

    for i in xrange(max_elements - 1):
        for j in xrange(max_isotopes):
            c = ~np.isclose(ptA[cutoff_idx], F.ifft((fft_ptA / (tAs[i])) * ((tAs[i] - itAs[i, j]))).real[cutoff_idx])
            contributions[i, j, c] = np.round((F.ifft((fft_ptA / (tAs[i])) * (itAs[i, j])).real[cutoff_idx][c] / ptA[cutoff_idx][c]) * element_vector[i])

    return np.array([MA[cutoff_idx], ptA[cutoff_idx]]), contributions