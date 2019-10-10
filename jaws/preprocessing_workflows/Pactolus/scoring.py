"""
Score spectra/scans against a collection of molecular fragmentation directed acyclic graphs (trees).

"""

__authors__ = 'Curt R. Fischer, Oliver Ruebel, Benjamin P. Bowen'
__copyright__ = 'Lawrence Berkeley National Laboratory and Authors, 2015.  All rights currently reserved.'


# standard libraries
import re
import os
import sys
import time
import datetime

# numpy and scipy
import numpy as np
from scipy.stats import norm
from scipy.sparse import csc_matrix


# rdkit for finding masses and inchi keys
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

# global variables
PATH_TO_TEST_DATA = './test_data/'

HIT_TABLE_DTYPE = np.dtype({'names': ['score', 'id', 'name',  'mass', 'n_peaks', 'n_match'],
                            'formats': ['f4', 'a100', 'a100', 'f4',   'i4',       'i4']})
"""
Numpy data type (dtype) used for hit tables
"""

FILE_LOOKUP_TABLE_DTYPE = np.dtype([('filename', 'a400'),
                                    ('ms1_mass', 'f8'), ])
"""
Numpy data dtype (dtype) used for the tree file lookup table
"""

METACYC_DTYPE = np.dtype([('metacyc_id', 'a30'),
                          ('name', 'a100'),
                          ('inchi', 'a1000'),
                          ('lins', 'a400'),
                          ('inchi_key', 'a40'),
                          ('mass', float)])
"""
Numpy data type (dtype) used to store METACYC molecular databases
"""


def calculate_lambda(fragment_masses, data_masses, mass_tol):
    """
    Calculates lambda values as described in the MIDAS paper [dx.doi.org/10.1021/ac5014783]

    :param fragment_masses:    np.ndarray of floats w/ shape (n,), theoretical masses from frag tree, in Da
    :param data_masses:        np.ndarray of floats w/ shape (n,), (neutralized) masses in Da detected in MS2 scan
    :param mass_tol:           float, maximum allowed mass difference in Da between theoretical and detected peaks
    :return lambda:            np.ndarray of floats containing the lambda scores
    """
    epsilons = fragment_masses - data_masses
    return np.piecewise(epsilons,
                        condlist=np.abs(epsilons) <= mass_tol,
                        funclist=[lambda x: 2 * (1-norm.cdf(abs(x), scale=mass_tol/2)),
                                  0,
                                  ]
                        )


def bfs_plausibility_score(plaus_score, tree, matches, nodes=None,):
    """
    Modifies in place a numpy vector of plausibility scores for each fragment in a tree through breadth-first search.

    :param plaus_score:     numpy 1D vector with len same as tree, value is plausibility score
    :param tree:            numpy structured array output by FragTree.numpy_tree
    :param matches:         numpy 1D vector of ints; the rows of tree that match any data peaks
    :param nodes:           numpy array with 1 dimension; row indices of tree currently being scored
    """
    # if nodes is None, start at the root:
    if nodes is None:
        nodes = np.where(tree['parent_vec'] == -1)[0]
        plaus_score[nodes] = 1

    # find direct children of current nodes and if there are any, score them
    children = np.where(np.in1d(tree['parent_vec'], nodes))[0]
    if children.size:
        parents = tree['parent_vec'][children]
        # depth of frag i = num of bonds broken i.e. num of Trues in bond_bool_arr[i, :]
        depths = (tree[children]['bond_bool_arr']).sum(axis=1)
        parent_depths = (tree[parents]['bond_bool_arr']).sum(axis=1)
        b = depths - parent_depths

        base = np.select([np.in1d(parents, matches)], [0.5], default=0.1,)
        plaus_score[children] = np.power(base, b) * plaus_score[parents]

        # continue recursive bfs by looking for children of these children
        bfs_plausibility_score(plaus_score, tree, matches, nodes=children)

    else:
        return


def find_matching_fragments(data_masses, tree, mass_tol):
    """
    Find node sets in a tree whose mass is within mass_tol of a data_mz value

    :param data_masses: numpy 1D array, float, *neutralized* m/z values of data from an MS2 or MSn scan
    :param tree: numpy structured array as output by FragDag
    :param mass_tol: precursore m/z mass tolerance

    :return: matching_frag_sets, list of lists; len is same as data_mzs; sublists are idxs to rows of tree that match
    :return: unique_matching_frags, numpy 1d array, a flattened numpy version of unique idxs in matching_frag_sets
    """
    # start_idx is element for which inserting data_mz directly ahead of it maintains sort order
    start_idxs = np.searchsorted(tree['mass_vec'], data_masses-mass_tol)

    # end_idx is element for which inserting data_mz directly after it maintains sort order
    #  found by searching negative reversed list since np.searchshorted requires increasing order
    length = len(tree)
    end_idxs = length - np.searchsorted(-tree['mass_vec'][::-1], -(data_masses+mass_tol))

    # if the start and end idx is the same, the peak is too far away in mass from the data and will be empty
    matching_frag_sets = [range(start, end) for start, end in zip(start_idxs, end_idxs)]  # if range(start, end)]

    # flattening the list
    unique_matching_frags = np.unique(np.concatenate(matching_frag_sets))

    # Returning both the flat index array and the sets of arrays is good:
    #       matching_frag_sets makes maximizing the MIDAS score easy
    #       unique_matching_frags makes calculating the plausibility score easy
    return matching_frag_sets, unique_matching_frags


def normalize_intensities(mz_intensity_arr, order=1):
    """
    Normalizes intensities in a 2D numpy array of m/z values & intensities.
    Designed to be used non-neutralized data.

    :param mz_intensity_arr: numpy float with shape (num_peaks, 2)
    :param order: int, if 1 then L1 norm is returned; otherwise L2 norm is computed
    :return: out, a normalized version of mz_intensity_arr with intensities that sum to 1

    Note: this function assumes intensities can never be negative.
    """
    out = mz_intensity_arr
    norm_ = out[:, 1].sum() if order == 1 else np.linalg.norm(out[:, 1])
    out[:, 1] = out[:, 1] / norm_
    return out


def neutralize_mzs(data_mz_intensity_arr, neutralizations):
    """
    Converts data m/z values to possible data neutral mass values (in Da) by combinatorial subtraction from ionizations.

    :param data_mz_intensity_arr: numpy float with shape (num_peaks, 2), assumed normalized i.e. x[:, 1].sum() = 1
    :param neutralizations: masses: in Da of singly charged ionic species whose removal results in a neutral mass.
    :return: neutralized_mass_intensity_arr
    """
    num_peaks, _ = data_mz_intensity_arr.shape
    mzs, intensities = data_mz_intensity_arr[:, 0], data_mz_intensity_arr[:, 1]
    shifted_arrays = tuple(np.array([mzs+shift, intensities]).T for shift in neutralizations)
    return np.vstack(shifted_arrays)


def calculate_MIDAS_score(mz_intensity_arr, tree, mass_tol, neutralizations, max_depth=None, want_match_matrix=False):
    """
    Score the the plausibility that a given MS2 (or MSn) scan arose from a given compound.

    :param mz_intensity_arr:    ndarray w/ shape (n_peaks, 2).  m/z values in col 0, intensities in col 1
    :param tree:                numpy structured array for a frag. tree of a molecule, usually from FragTree.numpy_tree
    :param mass_tol:            maximum mass in Da by which two MS2 (or MSn) peaks can differ
    :param neutralizations:     list of floats, adjustments (in Da) added to data peaks in order to neutralize them
    :param max_depth:           int, truncates tree, keeping only nodes <= max_depth bond breakages away from root
    :param want_match_matrix:   bool, if True then tuple of (score, match_matrix) is returned, else return score only
    :return:                    score, float, a MIDAS score
    :return:                    match_matrix, a bool matrix with a shape of (n_tree_nodes, n_peaks).
                                              elements are True if given peak matches given node of frag_dag
                                              Returned only if want_match_matrix is True, othewise None is returned.
    """
    # NOTE: in the following we refer to axis=1 as rows and axis=0 as columns
    # attempt to truncate tree if max_depth is supplied
    if max_depth:
        node_depths = tree['bond_bool_arr'].sum(axis=1)
        nodes_to_keep = np.where(node_depths <= max_depth)[0]
        tree = tree[nodes_to_keep]

    # normalize and neutralize data
    mz_rel_intensity_arr = normalize_intensities(mz_intensity_arr)
    mass_rel_intensity_arr = neutralize_mzs(mz_rel_intensity_arr, neutralizations)
    data_masses, data_rel_intensities = mass_rel_intensity_arr[:, 0], mass_rel_intensity_arr[:, 1]

    # find matching fragments
    matching_frag_sets, unique_frag_arr = find_matching_fragments(data_masses, tree, mass_tol)

    # if there are no matching fragments, the score is 0
    if unique_frag_arr.size == 0:
        score = 0
        match_matrix = None
    # Compute the score and match matrix
    else:
        # initialize and modify in place the plaus_score array:
        plaus_score = np.zeros(len(tree))
        bfs_plausibility_score(plaus_score, tree, matches=unique_frag_arr)

        # construct match_matrix to matrix (csr format) where columns are data peaks, rows are tree nodes.
        # matrix is constructed as sparse but converted to dense.
        # TODO:  slight speedup probably possible for full sparse implementation.
        inds = np.concatenate(matching_frag_sets)  # row indices
        indptr = np.append(0, np.cumsum(([len(el) for el in matching_frag_sets])))
        data = np.ones(shape=inds.shape, dtype=bool)
        match_matrix = csc_matrix((data, inds, indptr), dtype=bool).toarray()

        # loop over rows of match_matrix to calculate score
        score = 0
        rows_with_matches = np.where(match_matrix.any(axis=1))[0]
        for row_idx in rows_with_matches:
            col_idx = np.where(match_matrix[row_idx, :])[0]
            lambdas = calculate_lambda(tree[row_idx]['mass_vec'], data_masses[col_idx], mass_tol)
            subscore = (plaus_score[row_idx] * data_rel_intensities[col_idx] * lambdas).max()
            score += subscore

    # Return the reuqested results
    if want_match_matrix:
        #if match_matrix is not None:
        #    print (match_matrix.shape, mz_intensity_arr.shape, tree.shape)
        return score, match_matrix
    else:
        return score, None



def load_file_lookup_table(path):
    """
    Load or generate a file lookup table based on the given path.

    :param path: This may be either the path to: 1) the .npy file with the file lookup table prepared via
                 make_file_lookup_table_by_MS1_mass, 2) a path to an HDF5 file ending in '.h5'
                 and dataset defined via <filename>:<datasetpath
                 3) a textfile where each line a path to a tree
                 4) the path to a directory with all trees that should be used. In the latter 2 cases the
                 make_file_lookup_table_by_MS1_mass function is called to compile the lookup table, whereas
                 in the first case the table is simply restored from file.

    :return: tree_files_by_ms1_mass, a numpy structured array with columns (i) filename and (ii) MS1 mass.
             The numpy dtype is type = [('filename', 'a400'),('ms1_mass', 'f8'), ]

    """
    if os.path.isfile(path):
        if path.endswith('.npy'):
            file_lookup_table = np.load(path)
        else:
            in_treefile = open(path, 'r')
            tree_files = [line.rstrip('\n') for line in in_treefile]
            in_treefile.close()
            file_lookup_table = make_file_lookup_table_by_MS1_mass(tree_files=tree_files)
            in_treefile.close()
    elif os.path.isdir(path):
        file_lookup_table = make_file_lookup_table_by_MS1_mass(path=path)
    else:
        split_path = path.split(':')
        if len(split_path) == 2 and split_path[0].endswith('.h5'):
            hdf_file = h5py.File(split_path[0], 'r')
            file_lookup_table = hdf_file[split_path[1]][:]
        else:
            file_lookup_table = None

    return file_lookup_table


def compound_metadata_from_trees(table_file,
                                 mpi_comm=None,
                                 mpi_root=0):
    """
    Compile metadata for the trees in the given table_file

    :param table_file:  Full path to .npy file containing tree file names and parent masses or the numpy array directly.
    :param mpi_comm: The MPI communicator to be used
    :param mpi_root: The MPI root rank to be used for writing

    :return: Dictionary of 1D arrays with metadata information about the compounds. This includes the:
            inchi, num_atoms, num_bonds, inchi_key, mass
    """
    # load .npy file & extract inchi keys from the filename
    if isinstance(table_file, basestring):
        file_lookup_table = np.load(table_file)
    elif isinstance(table_file, np.ndarray):
        file_lookup_table = table_file
    else:
        raise ValueError('Invalid table-file given')

    # Initalize the compound metadata return values
    num_trees = len(file_lookup_table)
    compound_metadata = {'inchi': np.zeros(num_trees, dtype='a1000'),
                         'num_atoms': np.zeros(num_trees, dtype='int'),
                         'num_bonds': np.zeros(num_trees, dtype='int'),
                         'inchi_key': np.zeros(num_trees, dtype='a100'),
                         'mass': np.zeros(num_trees, dtype='float'),}
                        # TODO Need the name
                        # TODO Need the id (or metacyc_id)
                        # TODO need the lins (???_

    # Compile the metadata based on the metadata information stored in the trees
    for treeindex, treefile in enumerate(file_lookup_table):

        filename = treefile['filename']
        file_reader = h5py.File(filename)
        group_key = file_reader.keys()[0]
        data_key = file_reader[group_key].keys()[0]
        compound_metadata['inchi_key'][treeindex] = group_key
        compound_metadata['mass'][treeindex] = file_reader[group_key][data_key][-1]['mass_vec']
        try:
            compound_metadata['inchi'][treeindex] = file_reader[group_key].attrs['inchi']
            compound_metadata['num_atoms'][treeindex] = file_reader[group_key].attrs['num_atoms']
            compound_metadata['num_bonds'][treeindex] = file_reader[group_key].attrs['num_bonds']
        except:
            log_helper.warning(__name__, "Could not determine metadata for tree: " + str(treeindex),
                               root=mpi_root, comm=mpi_comm)
        file_reader.close()

    # return the compound metadata results
    return compound_metadata




