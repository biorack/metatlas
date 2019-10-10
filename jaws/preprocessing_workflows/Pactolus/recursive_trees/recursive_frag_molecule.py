import numpy as np

# rdkit
from rdkit import Chem
#from rdkit.Chem.rdmolops import GetFormalCharge
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

def remove_bonds(rw_mol, mol, bond_indices, undo=False):
    """
    Helper function that modifies in place a rdkit.Chem.rdchem.RWMol object by removing bonds.

    Useful because RWMol.RemoveBond() method annoyingly takes beginning and ending atom indices,
    not bond indices, and also because there is no RWMol.RemoveBonds() vectorized version.
    """
    if hasattr(bond_indices,'__iter__'):
        bonds = [mol.GetBondWithIdx(idx) for idx in bond_indices]
    else:
        bonds = [mol.GetBondWithIdx(bond_indices)]
        
    bond_atom_indices = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in bonds]
    if not undo:
        [rw_mol.RemoveBond(atoms[0], atoms[1]) for atoms in bond_atom_indices]
    else:
        [rw_mol.AddBond(atoms[0], atoms[1]) for atoms in bond_atom_indices]
    return


###MolFragmentToSmiles.
# This might be a way to find fragments from regular pactolus trees
# It uses the atom coordinates and the original mol to produce a smiles for the subset of atoms

###Also, you could use EditableMol and Remove Atom

###Do substructure match: http://www.rdkit.org/docs/RDKit_Book.html#atom-atom-matching-in-substructure-queries
def recursive_tree(all_frags,all_frag_masses,relationships,mol):
    f_tree = FragTree(mol)
    #parent = mol
    parent = Chem.MolToSmiles(Chem.AddHs(mol),False)
    if not parent in all_frags:
        all_frags.append(parent)
        all_frag_masses.append(CalcExactMolWt(mol))
    #print len(all_frags),len(relationships)
    for i,f in enumerate(f_tree.fragment_list):
        if f['frag_mol'].GetNumBonds() > 6:
            fragment = f['frag_smiles'] #the identifier to look up the molecule with
            #fragment = Chem.MolToSmiles(Chem.AddHs(f['frag_mol']),False)
            if not fragment in all_frags:
                all_frags.append(fragment)
                all_frag_masses.append(f['frag_mass'])
            myrel = (all_frags.index(parent),all_frags.index(fragment))
            if not myrel in relationships:
                relationships.append(myrel)
                recursive_tree(all_frags,all_frag_masses,relationships,f['frag_mol'])
    return all_frags,all_frag_masses,relationships
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

    def __init__(self, mol, max_depth=1):
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
        # set isotopes according to the supplied isotope dictionary

        # fill hydrogens and set isotopes in H'ed molecule
        self.molH = Chem.AddHs(self.mol)

        root_bonds = [bond.GetIdx() for bond in self.mol.GetBonds()]
        self.num_bonds = len(root_bonds)

        # make the fragmentation tree (as a dictionary)
        self.max_depth = max_depth
        self.fragment_list = self.make_fragment_list()


    def make_fragment_list(self):
        """
        Uses an iterable to compute all nodes of a fragmentation tree, but without retaining parent/child info.

        This function is used modify self.fragment_dict in place in the ``__init__`` constructor

        :returns:  fragment_dict ; Python dict with

            **keys**:   (sorted) tuple of atom indices in the fragment

            **values**: Python dict with

                    * **keys**: 'path'
                    * **values**: (sorted) tuple of unique shortest bond breakage path giving rise to fragment
        """
        fragment_list = []
        wm = Chem.RWMol(self.molH)
        for bond in range(self.num_bonds):
            # Remove bonds from H'ed molecule
            remove_bonds(wm, self.mol, bond)
            for mol in Chem.GetMolFrags(wm,sanitizeFrags=False,asMols=True):
                fragment_list.append({'frag_mol_h': mol,
                                      'frag_mol': Chem.RemoveHs(mol),
                                      'frag_mass': CalcExactMolWt(mol),
                                      'frag_smiles': Chem.MolToSmiles(mol,True),
                                      'fragment_mass':CalcExactMolWt(mol)})

            # Restore broken bonds
            remove_bonds(wm, self.mol, bond, undo=True)
        return fragment_list
    