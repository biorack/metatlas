#!python
"""
Compute the fuzzy distance between centroided spectra

"""
__authors__ = 'Curt R. Fischer, Oliver Ruebel, Benjamin P. Bowen'
__copyright__ = 'Lawrence Berkeley National Laboratory and Authors, 2015.  All rights currently reserved.'


import numpy as np

from score_frag_dag import normalize_intensities
from scipy.stats import norm


def calc_lambda(mz_1, mz_2, mass_tol):
    """
    Finds the mass-uncertainty-adjusted product of intensities between two peaks close in mass.

    :param mz_1: float, m/z value in Da for peak 1
    :param mz_2: float, m/z value in Da for peak 2
    :param mass_tol:, float, maximum allowable mass difference in Da for two peaks to be matched.
    :return: lambda value, float, mass-difference-adjusted equality between two mass values.
    The name lambda refers to http://pubs.acs.org/doi/abs/10.1021/ac5014783, not to Python's "lambda function" concept.
    """
    epsilon = abs(mz_1 - mz_2)
    return 0 if epsilon > mass_tol else 2 * (1-norm.cdf(abs(epsilon), scale=mass_tol/2))


def sparse_uncertain_dot(mz_intensity_arr_1, mz_intensity_arr_2, mass_tol):
    """
    Computes dot product of two (sparse) mass spectra accounting for mass uncertainty.

    :param mz_intensity_arr_1: numpy ndarray of shape (n_peaks_1, 2) with columns m/z and intensitiy
    :param mz_intensity_arr_2: numpy ndarray of shape (n_peaks_2, 2) with columns m/z and intensitiy
    :param mass_tol: float, maximum allowable mass difference in Da for two peaks to be matched.
    The name lambda refers to http://pubs.acs.org/doi/abs/10.1021/ac5014783, not to Python's "lambda function" concept.
    """
    mzs_1, intensities_1 = mz_intensity_arr_1[:, 0], mz_intensity_arr_1[:, 1]
    mzs_2, intensities_2 = mz_intensity_arr_2[:, 0], mz_intensity_arr_2[:, 1]
    unc_peak_prod_ufunc = np.frompyfunc(lambda x, y: calc_lambda(x, y, mass_tol), 2, 1)
    lambda_mat = np.ufunc.outer(unc_peak_prod_ufunc, mzs_1, mzs_2)
    return np.dot(intensities_1, lambda_mat).dot(intensities_2)


def calc_dot_matrix(scan_list, params):
    """
    Computes the sparse uncertain dot product for every unique pair of scans in a list of scans and
    returns the proximity matrix, i.e. matrix of dot products.

    :param scan_list: list of numpy ndarrays, each of shape (n_rows, 2) with columns m/z and intensity
    :param params:               dict, has keys:  'mass_tol': a float with allowed mass discrepancy in Da
                                                  'neutral_losses': list of floats [or dict with float values],
                                                                    allowed neutral losses.  When peaks in
                                                                    two spectra differ by any of these amounts,
                                                                    they are still counted as matching
                                                  'noise_cutoff':          (optional) float, value below which
                                                                                             intensities will be
                                                                                             set to zero
                                                  'normalize': int or None, whether or not to normalize intensities,
                                                                            if int indicates order of norm to use.
                                                                            For computing dot products L2 norm should
                                                                            be used.  For now, only values supported
                                                                            are 1 and 2.
                                                  'want_match_matrix': boolean, whether or not to return a matrix
                                                                                of dimension n by m, where n = n1 + n2 =
                                                                                number of unique masses in both spectra,
                                                                                and m is number of neutralizations
                                                  'metric': function to use for comparing spectra
    :return: proximity matrix, a numpy ndarray of shape (n_scans, n_scans) of floats,
                                non-zero only in upper triangular part
    """
    n_scans = len(scan_list)
    matrix = np.zeros(shape=(n_scans, n_scans), dtype=float)
    for i in xrange(n_scans):
        for j in xrange(i, n_scans):
            matrix[i, j] = score_spectra_pairwise(scan_list[i], scan_list[j], params)
    return matrix


def threshold_intensities(mz_intensity_arr, cutoff_intensity):
    """
    Remove noise from a mass spectrum by removing all peaks below a provided intensity value

    :param mz_intensity_arr: numpy nd_array with shape (n_peaks, 2) with columns m/z (in Da) and intensity in counts
    :param cutoff_intensity: float, value (in counts i.e. same units as intensity) below which all peaks are removed
    :return: thresholded_mz_intensity_arr
    """

    return mz_intensity_arr[mz_intensity_arr[:, 1] >= cutoff_intensity]


def score_spectra_pairwise(mz_intensity_arr_1, mz_intensity_arr_2, params=None):
    """
    Uses mass-tolerance-aware

    :param mz_intensity_arr_1:   numpy nd_array with shape (n_peaks_1, 2) with columns m/z (in Da) and intensity
    :param mz_intensity_arr_2:   numpy nd_array with shape (n_peaks_2, 2) with columns m/z (in Da) and intensity
    :param params:               dict, has keys:  'mass_tol': a float with allowed mass discrepancy in Da
                                                  'neutral_losses': list of floats [or dict with float values],
                                                                    allowed neutral losses.  When peaks in
                                                                    two spectra differ by any of these amounts,
                                                                    they are still counted as matching
                                                  'noise_cutoff':          (optional) float, value below which
                                                                                             intensities will be
                                                                                             set to zero
                                                  'normalize': int or None, whether or not to normalize intensities,
                                                                            if int indicates order of norm to use.
                                                                            For computing dot products L2 norm should
                                                                            be used.  For now, only values supported
                                                                            are 1 and 2.
                                                  'want_match_matrix': boolean, whether or not to return a matrix
                                                                                of dimension n by m, where n = n1 + n2 =
                                                                                number of unique masses in both spectra,
                                                                                and m is number of neutralizations
                                                  'metric': function to use for comparing spectra
    :return:  score or (score, match_matrix)
    """
    # default parameters:
    default_params = {
                      'mass_tol': 0.05,
                      'neutral_losses': {None: 0, '13C': 1.0034},
                      'noise_cutoff': None,
                      'normalize': 2,
                      'want_match_matrix': False,
                      'metric': None
                      }

    if params is None:
        params = default_params

    # check integrity of params
    if set(params).difference(set(default_params)):
        raise ValueError('Unknown parameter supplied in params dictionary')

    # define any unspecified parameters
    for undefined_parameter in set(default_params).difference(set(params)):
        params[undefined_parameter] = default_params[undefined_parameter]

    # unpack parameters
    scan_1, scan_2 = mz_intensity_arr_1, mz_intensity_arr_2
    mass_tol = float(params['mass_tol'])
    neutral_losses = params['neutral_losses']
    noise_cutoff = params['noise_cutoff']
    normalize = params['normalize']
    want_match_matrix = params['want_match_matrix']
    metric = params['metric']

    # TODO
    if metric is not None:
        raise NotImplementedError

    # do noise cutoff if requested:
    if noise_cutoff:
        scan_1 = threshold_intensities(scan_1, noise_cutoff)
        scan_2 = threshold_intensities(scan_2, noise_cutoff)

    # do intensity normalization if requested
    if normalize:
        scan_1 = normalize_intensities(scan_1, order=normalize)
        scan_2 = normalize_intensities(scan_2, order=normalize)

    # duplicate peaks in scan 2 by adding/subtracting every neutral loss to it
    if isinstance(neutral_losses, dict):
        neutral_losses = neutral_losses.values()

    old_scan_2 = scan_2.copy()
    for loss in neutral_losses:
        if loss != 0:
            shifted_scan_2 = old_scan_2
            shifted_mzs = old_scan_2[:, 0] + loss
            shifted_scan_2[:, 0] = shifted_mzs
            scan_2 = np.vstack((scan_2, shifted_scan_2))

    return sparse_uncertain_dot(scan_1, scan_2, mass_tol)

