""" unit tests for add_msms_refs module """
# pylint: disable=missing-function-docstring,line-too-long

import pytest
import traitlets

from metatlas.tools import add_msms_ref


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
        inchi_key="OLXZPDWKRNYJJZ-RRKCRQDMSA-N",
        inchi="InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1",
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
            inchi="InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1",
        )
