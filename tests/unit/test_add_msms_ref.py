""" unit tests for add_msms_refs module """
# pylint: disable=missing-function-docstring,line-too-long

import json
import pytest
import traitlets

from rdkit import Chem

from metatlas.tools import add_msms_ref

INCHI = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"
INCHI_KEY = "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
SMILES = "C1[C@@H]([C@@H](CO)O[C@H]1n1cnc2c(N)ncnc12)O"


def tests_msms_ref01(mocker, compound):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    add_msms_ref.MsmsRef(
        database="my_db",
        name="2'-deoxyadenosine",
        spectrum=add_msms_ref.Spectrum(intensities=[1, 1.4, 2], mzs=[100, 101, 555]),
        decimal=4,
        precursor_mz=251.101839276,
        polarity="negative",
        adduct="[M-H]+",
        fragmentation_method="cid",
        collision_energy="60eV",
        instrument="ThermoTOF-3000",
        instrument_type="LC-ESI-QTOF",
        formula="C10H13N5O3",
        exact_mass=251.101839276,
        inchi_key=INCHI_KEY,
        inchi=INCHI,
    )


def tests_msms_ref02(mocker, compound):
    mocker.patch("metatlas.datastructures.metatlas_objects.retrieve", return_value=[compound])
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.MsmsRef(
            database="my_db",
            name="2'-deoxyadenosine",
            spectrum=add_msms_ref.Spectrum(intensities=[1, 1.4, 2], mzs=[100, 101, 555]),
            decimal=4,
            precursor_mz=251.101839276,
            polarity="negative",
            adduct="[M-H]+",
            fragmentation_method="cid",
            collision_energy="60eV",
            instrument="ThermoTOF-3000",
            instrument_type="LC-ESI-QTOF",
            formula="C10H13N5O3",
            exact_mass=251.101839276,
            inchi_key="xxx",
            inchi=INCHI,
        )


def test_is_inchi01():
    assert add_msms_ref.is_inchi(INCHI)
    assert not add_msms_ref.is_inchi("f{INCHI}BLAH")
    assert not add_msms_ref.is_inchi("")
    assert not add_msms_ref.is_inchi("InChI=")


def test_is_valid_inchi_pair():
    assert add_msms_ref.is_valid_inchi_pair(INCHI, INCHI_KEY)
    assert not add_msms_ref.is_valid_inchi_pair("", INCHI_KEY)
    assert not add_msms_ref.is_valid_inchi_pair(INCHI, "")
    assert not add_msms_ref.is_valid_inchi_pair(f"{INCHI}foobar!", INCHI_KEY)
    assert not add_msms_ref.is_valid_inchi_pair(INCHI, f"{INCHI_KEY}foobar!")


def test_is_valid_inchi_smiles_pair():
    assert add_msms_ref.is_valid_inchi_smiles_pair(INCHI, SMILES)
    assert not add_msms_ref.is_valid_inchi_smiles_pair("", SMILES)
    assert not add_msms_ref.is_valid_inchi_smiles_pair(INCHI, "")
    assert not add_msms_ref.is_valid_inchi_smiles_pair(f"{INCHI}foobar!", SMILES)
    assert not add_msms_ref.is_valid_inchi_smiles_pair(INCHI, f"{SMILES}foobar!")


def test_are_equal():
    mol1 = Chem.inchi.MolFromInchi(INCHI)
    mol2 = Chem.inchi.MolFromInchi("InChI=1S/H2O/h1H2")
    assert add_msms_ref.are_equal(mol1, mol1)
    assert add_msms_ref.are_equal(mol2, mol2)
    assert not add_msms_ref.are_equal(mol1, mol2)
    assert not add_msms_ref.are_equal(mol2, mol1)


def test_is_synonym():
    assert add_msms_ref.is_synonym("foobar", "FOO///bar///FooZoo///FooBar")
    assert add_msms_ref.is_synonym("foobar", "FOOBAR")
    assert add_msms_ref.is_synonym("FooBar", "foobar///bar///FooZoo///FooBeeear")
    assert not add_msms_ref.is_synonym("foobar", "")
    assert not add_msms_ref.is_synonym("FooBarz", "foobar///bar///FooZoo///FooBeeear")


def test_spectrum01():
    add_msms_ref.Spectrum(intensities=[1.2, 1, 4], mzs=[123, 145, 256.04])


def test_spectrum02():
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.Spectrum(intensities=[1.2, 1], mzs=[123, 145, 256.04])
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.Spectrum(intensities=[1.2, 1, 4], mzs=[123, 145])
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.Spectrum(intensities=[1.2, 1, 4], mzs=[])
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.Spectrum(intensities=[], mzs=[123])
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.Spectrum(intensities=[1], mzs=[-123])
    with pytest.raises(traitlets.TraitError):
        add_msms_ref.Spectrum(intensities=[1, 1], mzs=[123, 22])


def test_str_to_spectrum():
    spectrum1 = add_msms_ref.str_to_spectrum("[[123.456,145.789],[1.0,2.2]]")
    assert spectrum1.mzs == [123.456, 145.789]
    assert spectrum1.intensities == [1.0, 2.2]
    spectrum2 = add_msms_ref.str_to_spectrum("[ [123.456, 145.789], [1.0, 2.2] ]")
    assert spectrum2.mzs == [123.456, 145.789]
    assert spectrum2.intensities == [1.0, 2.2]
    with pytest.raises(json.JSONDecodeError):
        add_msms_ref.str_to_spectrum("foobar")
