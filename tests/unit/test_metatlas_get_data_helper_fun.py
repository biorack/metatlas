""" unit testing of helper functions """
# pylint: disable=missing-function-docstring
import pytest
from metatlas.io import metatlas_get_data_helper_fun as gdhf


def test_make_atlas_df(atlas_two_compounds):
    # pylint: disable=line-too-long
    expected = """{"inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"compound_name":{"0":"2\'-deoxyadenosine"},"rt_max":{"0":2.6964640054},"rt_min":{"0":1.6964640054},"rt_peak":{"0":2.1964640054},"rt_units":{"0":"min"},"detected_polarity":{"0":"positive"},"mz":{"0":252.1091393},"mz_tolerance":{"0":20.0},"mz_tolerance_units":{"0":"ppm"},"mono_isotopic_molecular_weight":{"0":251.101839276},"pubchem_compound_id":{"0":"13730"},"synonyms":{"0":"2\'-deoxyadenosine"},"inchi":{"0":"InChI=1S\\/C10H13N5O3\\/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7\\/h3-7,16-17H,1-2H2,(H2,11,12,13)\\/t5-,6+,7+\\/m0\\/s1"},"adduct":{"0":"[M+H]+"},"label":{"0":"2\'-deoxyadenosine"},"ms1_notes":{"0":"keep"},"ms2_notes":{"0":"bad match to ref"},"identification_notes":{"0":"my id note"}}"""  # noqa: E501
    expected = """{"inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","1":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"compound_name":{"0":"2\'-deoxyadenosine","1":"2\'-deoxyadenosine"},"rt_max":{"0":2.6964640054,"1":2.6964640054},"rt_min":{"0":1.6964640054,"1":1.6964640054},"rt_peak":{"0":2.1964640054,"1":2.1964640054},"rt_units":{"0":"min","1":"min"},"detected_polarity":{"0":"positive","1":"positive"},"mz":{"0":252.1091393,"1":252.1091393},"mz_tolerance":{"0":20.0,"1":20.0},"mz_tolerance_units":{"0":"ppm","1":"ppm"},"mono_isotopic_molecular_weight":{"0":251.101839276,"1":251.101839276},"pubchem_compound_id":{"0":"13730","1":"13730"},"synonyms":{"0":"2\'-deoxyadenosine","1":"2\'-deoxyadenosine"},"inchi":{"0":"InChI=1S\\/C10H13N5O3\\/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7\\/h3-7,16-17H,1-2H2,(H2,11,12,13)\\/t5-,6+,7+\\/m0\\/s1","1":"InChI=1S\\/C10H13N5O3\\/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7\\/h3-7,16-17H,1-2H2,(H2,11,12,13)\\/t5-,6+,7+\\/m0\\/s1"},"adduct":{"0":"[M+H]+","1":"[M+H]+"},"label":{"0":"2\'-deoxyadenosine","1":"2\'-deoxyadenosine"},"ms1_notes":{"0":"keep","1":"keep"},"ms2_notes":{"0":"bad match to ref","1":"bad match to ref"},"identification_notes":{"0":"my id note","1":"my id note"}}"""  # noqa: E501
    assert expected == gdhf.make_atlas_df(atlas_two_compounds).to_json()


def test_transfer_identification_data_to_atlas(metatlas_dataset, atlas):
    mod_atlas = atlas.clone(recursive=True)
    mod_atlas.compound_identifications[0].ms1_notes = "ms1_note to overwrite"
    mod_atlas.compound_identifications[0].ms2_notes = "ms2_note to overwrite"
    mod_atlas.compound_identifications[
        0
    ].identification_notes = "identification_note to overwrite"
    out = gdhf.transfer_identification_data_to_atlas(metatlas_dataset, atlas)
    updated = atlas.compound_identifications[0]
    assert updated.ms1_notes == out.compound_identifications[0].ms1_notes
    assert updated.ms2_notes == out.compound_identifications[0].ms2_notes
    assert (
        updated.identification_notes
        == out.compound_identifications[0].identification_notes
    )


def test_set_nested_term_attr(metatlas_dataset):
    gdhf.set_nested(
        metatlas_dataset,
        [0, 0, "identification", "mz_references", 0, "adduct"],
        "[M+NH4]+",
    )
    assert (
        metatlas_dataset[0][0]["identification"].mz_references[0].adduct == "[M+NH4]+"
    )


def test_set_nested_term_attr_tuple(metatlas_dataset):
    gdhf.set_nested(
        metatlas_dataset,
        [0, 0, "identification", "mz_references", 0, ("adduct",)],
        "[M+NH4]+",
    )
    assert (
        metatlas_dataset[0][0]["identification"].mz_references[0].adduct == "[M+NH4]+"
    )


def test_set_nested_term_list(metatlas_dataset):
    gdhf.set_nested(metatlas_dataset, [0, 0], "foobar")
    assert metatlas_dataset[0][0] == "foobar"


def test_set_nested_term_dict(metatlas_dataset):
    gdhf.set_nested(metatlas_dataset, [0, 0, "identification"], "foobar")
    assert metatlas_dataset[0][0]["identification"] == "foobar"


def test_set_nested_term_no_ids(metatlas_dataset):
    with pytest.raises(ValueError):
        gdhf.set_nested(metatlas_dataset, [], "foobar")
