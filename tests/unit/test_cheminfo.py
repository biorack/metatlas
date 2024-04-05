# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long

from matchms.filtering.filter_utils.load_known_adducts import load_known_adducts
from rdkit import Chem

from metatlas.tools import cheminfo

INCHI = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"
INCHI_KEY = "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
SMILES = "C1[C@@H]([C@@H](CO)O[C@H]1n1cnc2c(N)ncnc12)O"


def test_get_parent_mass01():
    adducts = load_known_adducts()
    original_parent = 100
    for name in adducts:
        pre = cheminfo.get_precursor_mz(original_parent, name)
        parent = cheminfo.get_parent_mass(pre, name)
        assert abs(original_parent - parent) < 1e-7


def test_is_valid_inchi_pair():
    assert cheminfo.is_valid_inchi_pair(INCHI, INCHI_KEY)
    assert not cheminfo.is_valid_inchi_pair("", INCHI_KEY)
    assert not cheminfo.is_valid_inchi_pair(INCHI, "")
    assert not cheminfo.is_valid_inchi_pair(f"{INCHI}foobar!", INCHI_KEY)
    assert not cheminfo.is_valid_inchi_pair(INCHI, f"{INCHI_KEY}foobar!")


def test_is_valid_inchi_smiles_pair():
    assert cheminfo.is_valid_inchi_smiles_pair(INCHI, SMILES)
    assert not cheminfo.is_valid_inchi_smiles_pair("", SMILES)
    assert not cheminfo.is_valid_inchi_smiles_pair(INCHI, "")
    assert not cheminfo.is_valid_inchi_smiles_pair(f"{INCHI}foobar!", SMILES)
    assert not cheminfo.is_valid_inchi_smiles_pair(INCHI, f"{SMILES}foobar!")


def test_are_equal():
    mol1 = Chem.inchi.MolFromInchi(INCHI)
    mol2 = Chem.inchi.MolFromInchi("InChI=1S/H2O/h1H2")
    assert cheminfo.are_equal(mol1, mol1)
    assert cheminfo.are_equal(mol2, mol2)
    assert not cheminfo.are_equal(mol1, mol2)
    assert not cheminfo.are_equal(mol2, mol1)


def test_is_synonym():
    assert cheminfo.is_synonym("foobar", "FOO///bar///FooZoo///FooBar")
    assert cheminfo.is_synonym("foobar", "FOOBAR")
    assert cheminfo.is_synonym("FooBar", "foobar///bar///FooZoo///FooBeeear")
    assert not cheminfo.is_synonym("foobar", "")
    assert not cheminfo.is_synonym("FooBarz", "foobar///bar///FooZoo///FooBeeear")
