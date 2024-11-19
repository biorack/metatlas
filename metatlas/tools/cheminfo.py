""" cheminformatics related functions """

import functools
import logging

import ipywidgets as widgets
from matchms import Spectrum, filtering
import numpy as np
from rdkit import Chem

from metatlas.interfaces.compounds import structure_cleaning as cleaning

logger = logging.getLogger(__name__)


@functools.lru_cache
def get_parent_mass(precursor_mz: float, adduct: str) -> float:
    """Returns the mass of the input molecule that would result in the supplied precursor_mz and adduct"""
    dummy = Spectrum(
        mz=np.array([]), intensities=np.array([]), metadata={"precursor_mz": precursor_mz, "adduct": adduct}
    )
    updated = filtering.metadata_processing.derive_ionmode.derive_ionmode(dummy)
    #updated = filtering.metadata_processing.correct_charge.correct_charge(updated) # Can be added to remove warning
    #updated = filtering.metadata_processing.add_parent_mass.add_parent_mass(updated) # Can be added to remove warning
    return updated.metadata["parent_mass"]


@functools.lru_cache
def get_precursor_mz(parent_mass: float, adduct: str) -> float:
    """For an input molecule with parent_mass that generates adduct, return the resutling precursor_mz"""
    adducts = filtering.filter_utils.load_known_adducts.load_known_adducts()
    if adduct not in adducts['adduct'].tolist():
        raise KeyError("Adduct '%s' is not supported")
    multiplier = adducts.loc[adducts['adduct'] == adduct, 'mass_multiplier'].values[0]
    correction_mass = adducts.loc[adducts['adduct'] == adduct, 'correction_mass'].values[0]
    return (parent_mass * multiplier) + correction_mass


@functools.lru_cache
def is_positive_mode(adduct: str) -> bool:
    """Returns True if the MS mode for an adduct is positive"""
    adducts = filtering.filter_utils.load_known_adducts.load_known_adducts()
    if adduct not in adducts['adduct'].values:
        raise KeyError("Adduct '%s' is not supported")
    return adducts.loc[adduct, "ionmode"] == "positive"


@functools.lru_cache
def is_valid_inchi_pair(test_inchi: str, test_inchi_key: str) -> bool:
    """True if if test_inchi has the inchi key test_inchi_key"""
    if not filtering.filter_utils.smile_inchi_inchikey_conversions.is_valid_inchi(test_inchi):
        return False
    return test_inchi_key == Chem.inchi.InchiToInchiKey(test_inchi)


@functools.lru_cache
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


@functools.lru_cache
def inchi_to_smiles(inchi: str) -> str:
    """Convert Inchi to smiles"""
    out = Chem.MolToSmiles(Chem.inchi.MolFromInchi(inchi))
    if out is None:
        raise ValueError(f"'{inchi}' is not a valid Inchi")
    return out


@functools.lru_cache
def smiles_to_inchi(smiles: str) -> str:
    """Convert smiles to Inchi"""
    out = Chem.inchi.MolFromInchi(Chem.MolFromSmiles(smiles))
    if out is None:
        raise ValueError(f"'{smiles}' is not a valid smiles")
    return out


@functools.lru_cache
def inchi_or_smiles_to_molecule(molecule_id: str) -> Chem.rdchem.Mol:
    """Convert Inchi or smiles to rdkit Mol"""
    out = Chem.inchi.MolFromInchi(molecule_id) or Chem.MolFromSmiles(molecule_id)
    if out is None:
        raise ValueError(f"'{molecule_id}' is not a valid Inchi or smiles")
    return out


@functools.lru_cache
def inchi_or_smiles_to_inchi(molecule_id: str) -> str:
    """Inchi or smiles string to smiles string"""
    out = Chem.inchi.MolToInchi(inchi_or_smiles_to_molecule(molecule_id))
    if out is None:
        raise ValueError(f"'{molecule_id}' is not a valid Inchi or smiles")
    return out


@functools.lru_cache
def inchi_or_smiles_to_smiles(molecule_id: str) -> str:
    """Inchi or smiles string to smiles string"""
    out = Chem.MolToSmiles(inchi_or_smiles_to_molecule(molecule_id))
    if out is None:
        raise ValueError(f"'{molecule_id}' is not a valid Inchi or smiles")
    return out


@functools.lru_cache
def normalize_molecule(mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """Removes salt and neutralizes charges"""
    desalted, _ = cleaning.desalt(mol)
    if desalted is not None:
        neutralized, _ = cleaning.NeutraliseCharges(desalted)
    else:
        return None
    if neutralized is not None:
        return neutralized
    else:
        return None


@functools.lru_cache
def are_equal(molecule1: Chem.rdchem.Mol, molecule2: Chem.rdchem.Mol) -> bool:
    """True if both molecules are substructures of each other"""
    return molecule1.HasSubstructMatch(molecule2) and molecule2.HasSubstructMatch(molecule1)


@functools.lru_cache
def is_synonym(name: str, synonym_string: str) -> bool:
    """
    Inputs:
        name: string to check for within synonym_string
        synonym_string: string with /// between names
    Returns True if case insensitive match of name to full name in synonym_string
    """
    return name.lower() in [x.lower() for x in synonym_string.split("///")]


@functools.lru_cache
def valid_adduct(value: str) -> bool:
    """
    True if the value is an adduct listed supported by the matchms package
    This is not a comprehensive list, so it will return False for some uncommon adducts
    """
    adducts = filtering.filter_utils.load_known_adducts.load_known_adducts()
    return value in adducts['adduct'].tolist()


def mol_to_image(mol: Chem.rdchem.Mol, **kwargs) -> widgets.Image:
    """Generated a ipywidget.Image of a molecule from a rkit Mol"""
    d2d = Chem.Draw.MolDraw2DSVG(300, 300)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    text = d2d.GetDrawingText()
    return widgets.Image(value=text.encode("utf-8"), format="svg+xml", **kwargs)
