""" cheminformatics related functions """

import logging

import matchms
import numpy as np

from rdkit import Chem

from metatlas.interfaces.compounds import structure_cleaning as cleaning

logger = logging.getLogger(__name__)


def get_parent_mass(precursor_mz: float, adduct: str) -> float:
    """Returns the mass of the input molecule that would result in the supplied precursor_mz and adduct"""
    dummy = matchms.Spectrum(
        mz=np.array([]), intensities=np.array([]), metadata={"precursor_mz": precursor_mz, "adduct": adduct}
    )
    updated = matchms.filtering.add_parent_mass(dummy)
    return updated.metadata["parent_mass"]


def get_precursor_mz(parent_mass: float, adduct: str) -> float:
    """For an input molecule with parent_mass that generates adduct, return the resutling precursor_mz"""
    adducts = matchms.importing.load_adducts_dict()
    if adduct not in adducts:
        raise KeyError("Adduct '%s' is not supported")
    multiplier = adducts[adduct]["mass_multiplier"]
    correction_mass = adducts[adduct]["correction_mass"]
    return (parent_mass + correction_mass) / multiplier


def is_positive_mode(adduct: str) -> bool:
    """Returns True if the MS mode for an adduct is positive"""
    adducts = matchms.importing.load_adducts_dict()
    if adduct not in adducts:
        raise KeyError("Adduct '%s' is not supported")
    return adducts[adduct]["ionmode"] == "positive"


def is_valid_inchi_pair(test_inchi: str, test_inchi_key: str) -> bool:
    """True if if test_inchi has the inchi key test_inchi_key"""
    if not matchms.utils.is_valid_inchi(test_inchi):
        return False
    return test_inchi_key == Chem.inchi.InchiToInchiKey(test_inchi)


def is_valid_inchi_smiles_pair(test_inchi: str, test_smiles: str) -> bool:
    """
    True if test_inchi and test_smiles have the same structure.
    """
    mol_from_inchi = Chem.inchi.MolFromInchi(test_inchi)
    if mol_from_inchi is None:
        return False
    mol_from_smiles = Chem.MolFromSmiles(test_smiles)
    if mol_from_smiles is None:
        return False
    return are_equal(mol_from_inchi, mol_from_smiles)


def inchi_to_smiles(inchi: str) -> str:
    """Convert Inchi to smiles"""
    out = Chem.MolToSmiles(Chem.inchi.MolFromInchi(inchi))
    if out is None:
        raise ValueError(f"'{inchi}' is not a valid Inchi")
    return out


def smiles_to_inchi(smiles: str) -> str:
    """Convert smiles to Inchi"""
    out = Chem.inchi.MolFromInchi(Chem.MolFromSmiles(smiles))
    if out is None:
        raise ValueError(f"'{smiles}' is not a valid smiles")
    return out


def inchi_or_smiles_to_molecule(molecule_id: str) -> Chem.rdchem.Mol:
    """Convert Inchi or smiles to rdkit Mol"""
    out = Chem.inchi.MolFromInchi(molecule_id) or Chem.MolFromSmiles(molecule_id)
    if out is None:
        raise ValueError(f"'{molecule_id}' is not a valid Inchi or smiles")
    return out


def inchi_or_smiles_to_inchi(molecule_id: str) -> str:
    """Inchi or smiles string to smiles string"""
    out = Chem.inchi.MolToInchi(inchi_or_smiles_to_molecule(molecule_id))
    if out is None:
        raise ValueError(f"'{molecule_id}' is not a valid Inchi or smiles")
    return out


def inchi_or_smiles_to_smiles(molecule_id: str) -> str:
    """Inchi or smiles string to smiles string"""
    out = Chem.MolToSmiles(inchi_or_smiles_to_molecule(molecule_id))
    if out is None:
        raise ValueError(f"'{molecule_id}' is not a valid Inchi or smiles")
    return out


def normalize_molecule(mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """Removes salt and neutralizes charges"""
    desalted, _ = cleaning.desalt(mol)
    return cleaning.NeutraliseCharges(desalted)[0]


def are_equal(molecule1: Chem.rdchem.Mol, molecule2: Chem.rdchem.Mol) -> bool:
    """True if both molecules are substructures of each other"""
    return molecule1.HasSubstructMatch(molecule2) and molecule2.HasSubstructMatch(molecule1)


def is_synonym(name: str, synonym_string: str) -> bool:
    """
    Inputs:
        name: string to check for within synonym_string
        synonym_string: string with /// between names
    Returns True if case insensitive match of name to full name in synonym_string
    """
    return name.lower() in [x.lower() for x in synonym_string.split("///")]


def valid_adduct(value: str) -> bool:
    """
    True if the value is an adduct listed supported by the matchms package
    This is not a comprehensive list, so it will return False for some uncommon adducts
    """
    adducts = matchms.importing.load_adducts_dict()
    return value in adducts
