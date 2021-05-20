""" unit testing of helper functions """
# pylint: disable=missing-function-docstring
import pytest
from metatlas.io import metatlas_get_data_helper_fun as gdhf


def test_make_atlas_df(atlas_with_2_cids):
    # pylint: disable=line-too-long
    expected = """{"inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","1":"OIRDTQYFTABQOQ-KQYNXXCUSA-N"},"compound_name":{"0":"2\'-deoxyadenosine","1":"adenosine"},"rt_max":{"0":2.6964640054,"1":3.523318408},"rt_min":{"0":1.6964640054,"1":2.523318408},"rt_peak":{"0":2.1964640054,"1":3.023318408},"rt_units":{"0":"min","1":"min"},"detected_polarity":{"0":"positive","1":"positive"},"mz":{"0":252.1091393,"1":268.1040539},"mz_tolerance":{"0":20.0,"1":20.0},"mz_tolerance_units":{"0":"ppm","1":"ppm"},"mono_isotopic_molecular_weight":{"0":251.101839276,"1":267.096753896},"pubchem_compound_id":{"0":"13730","1":"60961"},"synonyms":{"0":"2\'-deoxyadenosine","1":"adenosine\\/\\/\\/58-61-7\\/\\/\\/Adenocard\\/\\/\\/Adenoscan"},"inchi":{"0":"InChI=1S\\/C10H13N5O3\\/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7\\/h3-7,16-17H,1-2H2,(H2,11,12,13)\\/t5-,6+,7+\\/m0\\/s1","1":"InChI=1S\\/C10H13N5O4\\/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10\\/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)\\/t4-,6-,7-,10-\\/m1\\/s1"},"adduct":{"0":"[M+H]+","1":"[M+H]+"},"label":{"0":"2\'-deoxyadenosine","1":"adenosine"},"ms1_notes":{"0":"keep","1":""},"ms2_notes":{"0":"bad match to ref","1":""},"identification_notes":{"0":"my id note","1":""}}"""  # noqa: E501
    assert expected == gdhf.make_atlas_df(atlas_with_2_cids).to_json()


def test_transfer_identification_data_to_atlas(metatlas_dataset, atlas):
    mod_atlas = atlas.clone(recursive=True)
    mod_atlas.compound_identifications[0].ms1_notes = "ms1_note to overwrite"
    mod_atlas.compound_identifications[0].ms2_notes = "ms2_note to overwrite"
    mod_atlas.compound_identifications[0].identification_notes = "identification_note to overwrite"
    out = gdhf.transfer_identification_data_to_atlas(metatlas_dataset, atlas)
    updated = atlas.compound_identifications[0]
    assert updated.ms1_notes == out.compound_identifications[0].ms1_notes
    assert updated.ms2_notes == out.compound_identifications[0].ms2_notes
    assert updated.identification_notes == out.compound_identifications[0].identification_notes


def test_set_nested_term_attr(metatlas_dataset):
    gdhf.set_nested(
        metatlas_dataset,
        [0, 0, "identification", "mz_references", 0, "adduct"],
        "[M+NH4]+",
    )
    assert metatlas_dataset[0][0]["identification"].mz_references[0].adduct == "[M+NH4]+"


def test_set_nested_term_attr_tuple(metatlas_dataset):
    gdhf.set_nested(
        metatlas_dataset,
        [0, 0, "identification", "mz_references", 0, ("adduct",)],
        "[M+NH4]+",
    )
    assert metatlas_dataset[0][0]["identification"].mz_references[0].adduct == "[M+NH4]+"


def test_set_nested_term_list(metatlas_dataset):
    gdhf.set_nested(metatlas_dataset, [0, 0, "identification", "mz_references"], None)
    assert metatlas_dataset[0][0]["identification"].mz_references is None


def test_set_nested_term_dict(metatlas_dataset):
    gdhf.set_nested(metatlas_dataset, [0, 0, "identification"], "foobar")
    assert metatlas_dataset[0][0]["identification"] == "foobar"


def test_set_nested_term_no_ids(metatlas_dataset):
    with pytest.raises(ValueError):
        gdhf.set_nested(metatlas_dataset, [], "foobar")
