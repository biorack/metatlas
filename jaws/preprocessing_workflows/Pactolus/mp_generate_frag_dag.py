#!python
"""
Generate fragmentation directed acyclic graph (ie: trees) for molecules.

"""
__authors__ = 'Curt R. Fischer, Oliver Ruebel, Benjamin P. Bowen'
__copyright__ = 'Lawrence Berkeley National Laboratory and Authors, 2015.  All rights currently reserved.'



# standard library
import sys
import os
from itertools import chain, combinations, repeat
from collections import OrderedDict
from random import choice
from string import lowercase
from time import time

# command line args
import argparse
# from pactolus.third_party.argparse_helper import RawDescriptionDefaultHelpArgParseFormatter

# other common libraries
import csv

# numpy
import numpy as np

# rdkit
from rdkit import Chem
from rdkit.Chem.rdmolops import GetFormalCharge
from rdkit.Chem import rdqueries

# h5py
import h5py

import multiprocessing as mp
"""
----------------
HELPER FUNCTIONS
----------------
"""


def nonprefix_combinations(iterable, max_len):
    """
    Create list of tuples of sorted integers of length 1 to max_len while
    excluding tuples that are a prefix of another tuple.

    :param iterable:
    :param max_len: Maximum tuple length
    """
    # http://stackoverflow.com/a/31282125/4480692
    iterable = list(iterable)
    if not (iterable and max_len):
        return
    for comb in combinations(iterable, max_len):
        yield comb
    for length in xrange(max_len-2, -1, -1):
        for comb in combinations(iterable[:-1], length):
            yield comb + (iterable[-1],)


def remove_bonds(rw_mol, mol, bond_indices, undo=False):
    """
    Helper function that modifies in place a rdkit.Chem.rdchem.RWMol object by removing bonds.

    Useful because RWMol.RemoveBond() method annoyingly takes beginning and ending atom indices,
    not bond indices, and also because there is no RWMol.RemoveBonds() vectorized version.
    """
    bonds = [mol.GetBondWithIdx(idx) for idx in bond_indices]
    bond_atom_indices = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in bonds]
    if not undo:
        [rw_mol.RemoveBond(atoms[0], atoms[1]) for atoms in bond_atom_indices]
    else:
        [rw_mol.AddBond(atoms[0], atoms[1]) for atoms in bond_atom_indices]
    return


def randomword(length):
    """
    Generates a random string of lowercase alphabetic characters.
    """
    return ''.join(choice(lowercase) for _ in xrange(length))


def get_isotope_dict(isostope_file='../data/max_abundance_isotopes.csv'):
    """
    Returns isotope_dict (element symbols (strings) are keys and mass numbers (ints) are values) by reading a file.
    """
    try:
        isotope_dict = {}
        with open(isostope_file, 'r') as csvfile:
            file_contents = csv.reader(csvfile)
            for line in file_contents:
                el, mass_number = tuple(line)
                isotope_dict[el] = int(mass_number)
        assert isotope_dict
        return isotope_dict
    except IOError:
        raise TypeError("""
                        A 'max_abundance_isotopes.csv'  file must be present in this directory;
                        (or you can pass in your own isotope_dict manually from within python).
                        """)
    except AssertionError:
        raise RuntimeError("""
                           Could not properly construct isotope_dict, check contents of
                           'max_abundance_isotopes.csv'
                           """)


"""
-----------------------------
PRIMARY FUNCTIONS AND CLASSES
-----------------------------
"""


class FragTree(object):
    """
    A class to represent fragmentation trees of rdkit.Chem.rdchem.Mol objects.

    Note that bonds broken come from self.mol, while calculation of correct monoisotopic
    masses requires self.molH. Thus, functionality of this code requires that Chem.AddHs()
    maintains bond and atom indices between mol and molH.  For all cases checked by the authors,
    this is true.

    :param mol: rdkit.Chem.rdchem.Mol object, the parent molecule of the fragmentation tree
                this molecule must (i) be fully connected by covalent bonds and
                (ii) have only _implicit_ hydrogens
    :param molH: rdkit.Chem.rdchem.Mol object, a version of mol with but with explicit hydrogens
    :param atom_masses: list of floats, the mass of each atom in molecule, stored outside rdkit for speed
    :param breakable_bonds: list of ints, indices of bonds which should be broken during tree formation
    :param isotope_dict: dictionary with element symbols as keys and mass numbers as values

    """
    def __init__(self, mol, isotope_dict, max_depth):
        """
        Creates a fragmentation tree for a given rdkit.Chem.rdchem.Mol object.

        :param mol: rdkit.Chem.rdchem.Mol object, the parent molecule of the fragmentation tree
                    this molecule must (i) be fully connected by covalent bonds and
                    (ii) have only _implicit_ hydrogens
        :param isotope_dict: dictionary with element symbols as keys and mass numbers (ints) as values
        """
        # perform checks on integrity of input molecule
        if len(Chem.GetMolFrags(mol, sanitizeFrags=False, asMols=False)) != 1:
            raise TypeError('Molecule must be fully connected by covalent bonds.')

        if Chem.RemoveHs(mol).GetNumAtoms() != mol.GetNumAtoms():
            raise TypeError('Molecule must have only implicit H atoms.')

#        #TODO This expression should allow for molecules with a permanent charge, but forbid molecules that can be neutralized
#        if GetFormalCharge(mol) != 0:
#            raise TypeError('Molecule must have overall neutral charge.')

        # begin initialization
        self.mol = mol
        self.isotope_dict = isotope_dict

        # set isotopes according to the supplied isotope dictionary
        self.set_isotopes(self.mol)

        # fill hydrogens and set isotopes in H'ed molecule
        self.molH = Chem.AddHs(self.mol)
        self.set_isotopes(self.molH)

        # to speed calculation of monoisotopic masses, store each atom's mass outside of rdkit
        root_atoms = [atom.GetIdx() for atom in self.molH.GetAtoms()]
        root_bonds = [bond.GetIdx() for bond in self.mol.GetBonds()]
        self.atom_masses = [self.molH.GetAtomWithIdx(i).GetMass() for i in root_atoms]

        # store number of bonds and atoms
        self.num_atoms = len(root_atoms)
        self.num_bonds = len(root_bonds)

        # make the fragmentation tree (as a dictionary)
        self.max_depth = max_depth
        self.fragment_dict = self.make_fragment_dict()
        self.annotate_fragments()

        # store length of tree
        self.num_frags = len(self.fragment_dict)

        # convert to numpy / binary for archival and rapid search
        self.numpy_tree = self.convert_dict_to_numpy()

    def set_isotopes(self, mol):
        """
        Modifies in place the isotope property of all atoms in a molecule, according to the supplied dictionary.
        """
        [atom.SetIsotope(self.isotope_dict[atom.GetSymbol()]) for atom in mol.GetAtoms()]
        return

    def calculate_monoisotopic_mass(self, atom_indices):
        """
        Calculate the monoisotopic mass of a subset of atoms in self.molH from a list of atom indices.

        :param atom_indices: a list of unique integers, indices of the atoms that are in the fragment.
        """
        return sum(self.atom_masses[idx] for idx in atom_indices)

    def make_fragment_dict(self):
        """
        Uses an iterable to compute all nodes of a fragmentation tree, but without retaining parent/child info.

        This function is used modify self.fragment_dict in place in the ``__init__`` constructor

        :returns:  fragment_dict ; Python dict with

            **keys**:   (sorted) tuple of atom indices in the fragment

            **values**: Python dict with

                    * **keys**: 'path'
                    * **values**: (sorted) tuple of unique shortest bond breakage path giving rise to fragment
        """
        bond_sets = chain.from_iterable(combinations(xrange(self.num_bonds), depth)
                                        for depth in range(self.max_depth+1))
        fragment_dict = {}

        wm = Chem.RWMol(self.molH)
        for bond_set in bond_sets:
            # Remove bonds from H'ed molecule
            remove_bonds(wm, self.mol, bond_set)

            # Get a list of atom sets: tuples with each tuple containing indices for atoms in each fragment
            frag_atoms_list = Chem.GetMolFrags(wm, sanitizeFrags=False, asMols=False)

            # store a fragment and bond breakage path in the dictionary if it is (a) not there already or
            # (b) the bond breakage path that generates the fragment is shorter than the one previously added
            fragment_dict.update({atoms: {'path': bonds}
                                  for (atoms, bonds) in zip(frag_atoms_list, repeat(bond_set))
                                  if atoms not in fragment_dict
                                  or len(fragment_dict[atoms]) > len(bonds)})

            # Restore broken bonds
            remove_bonds(wm, self.mol, bond_set, undo=True)
        return fragment_dict

    def annotate_fragments(self):
        """
        Modify self.fragment_dict in place, adding info on  the parent nodes and monoisotopic masses of fragments.

        :return: `fragment_dict`, which is a Python `dict` with:

            * `keys`:   (sorted) tuple of atom indices in the fragment
            * `values`: dict with the following keys: values:

                * ``path``:   (sorted) tuple of unique shortest bond breakage path giving rise to fragment
                * ``parent``: the 'path' of the parent node
                * ``mass`` :  the monoisotopic mass of the fragment
        """
        for fragment, frag_data in self.fragment_dict.iteritems():
            # Add mass
            frag_data['mass'] = self.calculate_monoisotopic_mass(fragment)

            # Add parent nodes
            parent_path = frag_data['path']
            # The parent is a prefix subset of the node's bond breakage path.
            # We test prefix subsets in decreasing length order.
            # For linear bonds, the first prefix subset always results in a hit.
            # For ring bonds, the parent may arise from a more truncated subset
            #    b/c multiple bonds must break to release the fragment.
            while parent_path:
                candidate_parents = [key for key, data in self.fragment_dict.items()
                                     if data['path'] == parent_path[:-1]
                                     and set(fragment).issubset(set(key))]
                if candidate_parents:
                    frag_data['parent'] = candidate_parents[0]
                    break
                else:
                    parent_path = parent_path[:-1]
                    continue

        # set the parent of the root node
        all_atoms = tuple(range(self.num_atoms))
        self.fragment_dict[all_atoms]['parent'] = None
        return

    def convert_dict_to_numpy(self):
        """
        Convert a fragment dict to mass-sorted numpy format for fast and easy searching, scoring, and storage.
        """
        # Sort dictionary by mass (this simplifies parent reference by index vs. sorting later)
        ordered_frag_dict = OrderedDict(sorted(self.fragment_dict.items(), key=lambda t: t[1]['mass']))

        # Define output dtype
        dtype = [('atom_bool_arr', 'b1', self.num_atoms),
                 ('bond_bool_arr', 'b1', self.num_bonds),
                 ('mass_vec', 'f8'),
                 ('parent_vec', 'i8'), ]

        atom_idxs, bond_idxs = xrange(self.num_atoms), xrange(self.num_bonds)

        # atom_bool_arr must be iterated over separately and first b/c filling parent_vec requires it
        # atom_bool_arr[i, j] is True if atom j is in fragment i
        atom_bool_arr = np.array([[i in key for i in atom_idxs] for key in ordered_frag_dict], dtype=bool)

        # initialize and fill other arrays / vectors
        # bond_bool_arr[i, j] is True if bond j must be broken to generate fragment i
        bond_bool_arr = np.zeros((self.num_frags, self.num_bonds), dtype=bool)
        mass_vec = np.zeros(self.num_frags, dtype=float)
        parent_vec = np.zeros(self.num_frags, dtype=int)

        for idx, (k, val) in enumerate(ordered_frag_dict.iteritems()):
            bond_bool_arr[idx, :] = [i in val['path'] for i in bond_idxs]
            mass_vec[idx] = val['mass']
            try:
                parent_vec[idx] = [val['parent'] == key for key in ordered_frag_dict].index(True)
            except ValueError:
                # The root will have None as a parent which is not a valid dict key
                parent_vec[idx] = -1

        # fill numpy structured array by column
        numpy_tree = np.zeros(shape=self.num_frags, dtype=dtype)
        numpy_tree['atom_bool_arr'] = atom_bool_arr
        numpy_tree['bond_bool_arr'] = bond_bool_arr
        numpy_tree['mass_vec'] = mass_vec
        numpy_tree['parent_vec'] = parent_vec
        return numpy_tree

        

def grow_tree_from_inchi(inchi,max_depth=None, isotope_dict=None, file_params=None):
    """
    Makes a FragTree from an InChI string; this is desiged to be parallizable over molecules by mpi.

    :param: inchi, string, an InChI string for a molecule
    """
    # make fragmentation tree for each inchi
    mol = Chem.MolFromInchi(inchi)

    file_name = file_params['output_hdf5_file_base'] + '_' + str(max_depth) + '_' + Chem.InchiToInchiKey(inchi) + '.h5'
    file_path = os.path.join(file_params['output_directory'], file_name)

    # check for existence of file and ignore if complete
    # TODO: check whether max_depth used >= current max_depth before skipping (instead of just an implicit !=)
    if os.path.exists(file_path):
        print 'SKIPPING FragTree for Molecule %s' % (inchi,)
        return

    # Check for success of InChI conversion
    if mol is None:

        with open(file_params['output_error_log'], 'a') as log_file:
            log_file.write('%s was not parsed by RDKit.  Is it valid InChI?\n' % inchi)
            print '%s was not parsed by RDKit.  Is it valid InChI?\n' % inchi
        return

    # Try to make the tree and if successful, record time it takes
    try:
        start = time()
        mol_tree = FragTree(mol, max_depth=max_depth, isotope_dict=isotope_dict)
        stop = time()
        print 'FragTree for Molecule %s made in %s seconds' % (inchi, stop-start)

        # Output results to hdf5 file
        with h5py.File(file_path, 'a') as f:
            inchi_key = Chem.InchiToInchiKey(inchi)
            hdf5_group = f.create_group(inchi_key)
            dset_name = 'FragTree_at_max_depth=%s' % max_depth
            hdf5_group.create_dataset(name=dset_name, data=mol_tree.numpy_tree)

            # set attributes
            hdf5_group.attrs['inchi'] = inchi
            hdf5_group.attrs['permanent_charge'] = 0
            #TODO: calculate and store the permanent charge automatically
            hdf5_group.attrs['num_atoms'] = mol_tree.num_atoms
            hdf5_group.attrs['num_bonds'] = mol_tree.num_bonds
            hdf5_group.attrs['num_fragments'] = mol_tree.num_frags
            hdf5_group.attrs['max_depth'] = mol_tree.max_depth
            hdf5_group.attrs['time_to_build'] = stop - start

    except TypeError as bad_mol_err:
        print bad_mol_err
        with open(file_params['output_error_log'], 'a') as log_file:
            log_file.write('InChI %s gave an error %s: \n' % (inchi, bad_mol_err.message))
    return


"""
----------------
COMMAND LINE USAGE
----------------
"""

def command_line_params():
    """
    Setup and parse the command line parameters

    :return: dict with command line parameters
    """
    parser = argparse.ArgumentParser(add_help=False)
                                     # formatter_class=RawDescriptionDefaultHelpArgParseFormatter)
    parser.add_argument("--inchi_file", '-inchi',
                        action='store',
                        default='test_FragTreeLibrary_inchis.txt',
                        type=str,
                        required=False,
                        help='Text file with the inchi string of the molecules to be processed.')
    parser.add_argument("--output_base_name", '-of',
                        action='store',
                        default='FragTreeLibrary_test_hdf5',
                        type=str,
                        required=False,
                        help='Base name of the output file(s)')
    parser.add_argument("--output_dir", '-od',
                        action='store',
                        default='results',
                        type=str,
                        required=False,
                        help='Directory where the results should be stored.')
    parser.add_argument("--error_log", '-el',
                        action='store',
                        default='FragTreeLibrary_test_log.txt',
                        type=str,
                        required=False,
                        help="Name of the output error log file")
    parser.add_argument("--isotope_file", '-iso',
                        action='store',
                        default='../data/max_abundance_isotopes.csv',
                        type=str,
                        required=False,
                        help="Text file with the maximum abundant isotopes. ")
    parser.add_argument("--max_depth", "-md",
                        action="store",
                        default=2,
                        type=int,
                        required=False,
                        help="Integer indicting the maximum fragmentation depth.")
    # Add the command line help argument
    parser.add_argument('-h', '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='show this help message and exit')
    # Parse the command line arguments and return the dict with the results
    return vars(parser.parse_args())


def print_arguments(cl_params):
    """
    Print the settings from the dict with command line arguments

    :param cl_params: Dict with the command line arguments
    """
    print ""
    print "Settings:"
    print "---------"
    if mpi_helper.get_rank() == 0:
        for param_key, param_val in cl_params.iteritems():
            print str(param_key) + " = " + str(param_val)
    print ""

def grow_tree_mp(arg):
    inchi = arg[0]
    md = arg[1]
    isotope_dict = arg[2]
    file_params = arg[3]
    print inchi
    print mp.current_process()

    grow_tree_from_inchi(inchi,max_depth=md, isotope_dict=isotope_dict, file_params=file_params)

def main():
    """
    Main function
    """
    sys.path.append(os.path.dirname(__file__))

    # Parse the command line arguments
    cl_params = command_line_params()

    # Create the file parameter for the FragTreeLibrary
    file_params = {'input_inchi_file': cl_params['inchi_file'],
                   'output_directory': cl_params['output_dir'],
                   'output_hdf5_file_base': cl_params['output_base_name'],
                   'output_error_log': cl_params['error_log']}

    isotope_dict = get_isotope_dict(isostope_file=cl_params['isotope_file'])

    # Make output directory if it does not exist
    if not os.path.isdir(file_params['output_directory']):
        try:
            os.mkdir(file_params['output_directory'])
        except OSError:
            # When executed in parallel it is possible that another rank already created the dir
            # in the meantime. We can safely ignore this error.
            if os.path.isdir(file_params['output_directory']):
                pass
            else:
                raise

    # Get isotope dictionary (if none was provided)
    if isotope_dict is None:
        isotope_dict = get_isotope_dict()
    else:
        isotope_dict = isotope_dict
       
    # make list of inchis
    inchi_list = []
    with open(file_params['input_inchi_file'], 'r') as inchi_file:
        for line in inchi_file:
            inchi_list.append(line.strip())
   
    with open(file_params['output_error_log'], 'w') as _:
        pass

#    for inchi in inchi_list:
#        grow_tree_from_inchi(inchi,max_depth=cl_params['max_depth'], isotope_dict=isotope_dict, file_params=file_params)
    mp_params = []
    for inchi in inchi_list:
        mol = Chem.MolFromInchi(inchi)
        q = rdqueries.AtomNumEqualsQueryAtom(6)
        try:
            # In the code above, you have to have at least one bond to break or it will crash.  
            if len(mol.GetAtomsMatchingQuery(q)) > 1:
                #grow_tree_from_inchi(inchi,max_depth=cl_params['max_depth'], isotope_dict=isotope_dict, file_params=file_params)
                mp_params.append((inchi,cl_params['max_depth'], isotope_dict, file_params))
        except:
            print(inchi)
            
    pool = mp.Pool(processes=10)
    pool.map(grow_tree_mp, mp_params)
    pool.close()
    # for p in mp_params:
        # grow_tree_mp(p)

    return
#test
if __name__ == "__main__":
    main()
