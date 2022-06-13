""" unit testing of helper functions """
# pylint: disable=missing-function-docstring
import pytest
import numpy as np

from metatlas.io import metatlas_get_data_helper_fun as gdhf


def test_make_atlas_df(atlas_with_2_cids):
    # pylint: disable=line-too-long
    expected = """{"inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","1":"OIRDTQYFTABQOQ-KQYNXXCUSA-N"},"compound_name":{"0":"2\'-deoxyadenosine","1":"adenosine"},"rt_max":{"0":2.6964640054,"1":3.523318408},"rt_min":{"0":1.6964640054,"1":2.523318408},"rt_peak":{"0":2.1964640054,"1":3.023318408},"rt_units":{"0":"min","1":"min"},"detected_polarity":{"0":"positive","1":"positive"},"mz":{"0":252.1091393,"1":268.1040539},"mz_tolerance":{"0":20.0,"1":20.0},"mz_tolerance_units":{"0":"ppm","1":"ppm"},"mono_isotopic_molecular_weight":{"0":251.101839276,"1":267.096753896},"pubchem_compound_id":{"0":"13730","1":"60961"},"synonyms":{"0":"2\'-deoxyadenosine","1":"adenosine\\/\\/\\/58-61-7\\/\\/\\/Adenocard\\/\\/\\/Adenoscan"},"inchi":{"0":"InChI=1S\\/C10H13N5O3\\/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7\\/h3-7,16-17H,1-2H2,(H2,11,12,13)\\/t5-,6+,7+\\/m0\\/s1","1":"InChI=1S\\/C10H13N5O4\\/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10\\/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)\\/t4-,6-,7-,10-\\/m1\\/s1"},"adduct":{"0":"[M+H]+","1":"[M+H]+"},"label":{"0":"2\'-deoxyadenosine","1":"adenosine"},"ms1_notes":{"0":"keep","1":""},"ms2_notes":{"0":"-1,bad match to ref","1":""},"identification_notes":{"0":"my id note","1":""}}"""  # noqa: E501
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


def test_set_nested_term_attr01(metatlas_dataset):
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


def test_extract_list01():
    assert gdhf.extract([1, 2], [0]) == 1


def test_extract_list02():
    assert gdhf.extract([1, 2], [5]) is None


def test_extract_list03():
    assert gdhf.extract([1, 2], [5], 99) == 99


def test_extract_list04():
    assert gdhf.extract([], [0], 99) == 99


def test_extract_list05():
    assert gdhf.extract("okay", [], 99) == "okay"


def test_extract_tuple01():
    assert gdhf.extract(("foo", "bar"), [1]) == "bar"


def test_extract_tuple02():
    assert gdhf.extract(("foo", "bar"), [-1]) == "bar"


def test_extract_tuple03():
    assert gdhf.extract(("foo", "bar"), [9], "zoop") == "zoop"


def test_extract_dict01():
    assert gdhf.extract({"foo": 1, "bar": 2}, ["bar"]) == 2


def test_extract_dict02():
    assert gdhf.extract({"foo": 1, "bar": 2}, ["blah"], "zoop") == "zoop"


def test_extract_default01():
    assert gdhf.extract({"foo": 1, "bar": 2}, ["foo", "bar"], "zoop") == "zoop"


def test_extract_metatlas_dataset01(metatlas_dataset):
    ids = [0, 0, "identification", ("compound",), 0, ("inchi_key",)]
    assert gdhf.extract(metatlas_dataset, ids, "zoop") == "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"


def test_extract_metatlas_dataset02(metatlas_dataset):
    ids = [0, 0, "identification", "compound", 0, "inchi_key"]
    assert gdhf.extract(metatlas_dataset, ids, "zoop") == "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"


def test_extract_metatlas_dataset03(metatlas_dataset):
    assert gdhf.extract(metatlas_dataset, ["foo"], "zoop") == "zoop"


def test_get_data_for_atlas_and_lcmsrun(atlas_df, df_container):
    result = gdhf.get_data_for_atlas_and_lcmsrun(atlas_df, df_container, 0.75, 0)
    assert result[0] == [
        {
            "num_ms1_datapoints": 74.0,
            "mz_peak": 252.1090087891,
            "rt_peak": 2.2922415733,
            "mz_centroid": 252.10896296693303,
            "rt_centroid": 2.2720730579808084,
            "peak_height": 2359861.25,
            "peak_area": 57016800.755859375,
        }
    ]
    assert result[1] == [
        {
            "mz": [
                252.1089324951,
                252.1090087891,
                252.1088104248,
                252.1090087891,
                252.10887146,
                252.1089324951,
                252.1089324951,
                252.1088256836,
                252.1088867188,
                252.1090393066,
                252.1089782715,
                252.1089630127,
                252.1089630127,
                252.1089782715,
                252.1090240479,
                252.1089782715,
                252.1090240479,
                252.1089324951,
                252.1090393066,
                252.1088867188,
                252.10887146,
                252.1089324951,
                252.1089630127,
                252.1089935303,
                252.1089172363,
                252.1089477539,
                252.1090545654,
                252.1089630127,
                252.1090240479,
                252.1090087891,
                252.1090393066,
                252.1090240479,
                252.1089935303,
                252.1090240479,
                252.1089630127,
                252.1090087891,
                252.1090240479,
                252.1089172363,
                252.1089019775,
                252.1089477539,
                252.1089324951,
                252.1089477539,
                252.1089477539,
                252.1089477539,
                252.1089782715,
                252.1088867188,
                252.1089172363,
                252.1089324951,
                252.1089782715,
                252.1089477539,
                252.1089172363,
                252.1089324951,
                252.1089630127,
                252.1088867188,
                252.1089630127,
                252.1085205078,
                252.1090545654,
                252.1089935303,
                252.1088104248,
                252.1086578369,
                252.1089935303,
                252.1085510254,
                252.1082763672,
                252.1082458496,
                252.1084136963,
                252.1092224121,
                252.1091766357,
                252.1092834473,
                252.1087493896,
                252.1112518311,
                252.1088409424,
                252.1086425781,
                252.1091766357,
                252.1094055176,
            ],
            "rt": [
                2.1030805111,
                2.1084616184,
                2.1139531136,
                2.1193552017,
                2.1248509884,
                2.1302509308,
                2.135682106,
                2.1411821842,
                2.1459801197,
                2.1513926983,
                2.1568279266,
                2.1622362137,
                2.1676549911,
                2.1730883121,
                2.179015398,
                2.1845297813,
                2.1900422573,
                2.1949694157,
                2.20002985,
                2.2055358887,
                2.2110378742,
                2.2165191174,
                2.2219588757,
                2.2273921967,
                2.2328462601,
                2.2382712364,
                2.2437169552,
                2.2492566109,
                2.2547125816,
                2.2601687908,
                2.2656960487,
                2.2704958916,
                2.2758042812,
                2.2813498974,
                2.2868082523,
                2.2922415733,
                2.2976748943,
                2.3031060696,
                2.308131218,
                2.313628912,
                2.3185498714,
                2.3239560127,
                2.3293914795,
                2.3349123001,
                2.3403663635,
                2.346799612,
                2.3522267342,
                2.3576600552,
                2.3631224632,
                2.3685662746,
                2.3740911484,
                2.3794057369,
                2.3848536015,
                2.3903660774,
                2.3953785896,
                2.4006638527,
                2.4062638283,
                2.411709547,
                2.4171659946,
                2.4226117134,
                2.4302260876,
                2.4357616901,
                2.4407405853,
                2.4461927414,
                2.451615572,
                2.4571509361,
                2.4627010822,
                2.4681572914,
                2.4735822678,
                2.4735822678,
                2.4787945747,
                2.4842174053,
                2.4896612167,
                2.495146513,
            ],
            "intensity": [
                312203.5,
                387914.59375,
                308308.5,
                334653.59375,
                339521.625,
                345527.21875,
                292437.34375,
                413614.53125,
                300285.28125,
                383848.71875,
                404313.21875,
                377231.34375,
                453965.5625,
                431327.0,
                523180.0625,
                510239.8125,
                631459.1875,
                807419.5,
                842647.5625,
                1053031.625,
                1082361.625,
                1198966.625,
                1109162.375,
                1126347.125,
                1373071.5,
                1589018.375,
                1281309.875,
                1660166.75,
                1492912.25,
                2029801.5,
                2029874.125,
                2035966.625,
                2010867.875,
                2036981.375,
                2148879.25,
                2359861.25,
                2054066.125,
                1691976.0,
                1778159.125,
                1776166.125,
                1752154.125,
                1575676.875,
                1199910.625,
                1259708.25,
                1087384.375,
                826077.125,
                802296.875,
                547785.125,
                545340.0625,
                584624.4375,
                468524.8125,
                305931.1875,
                330310.34375,
                309740.625,
                289212.71875,
                230440.9375,
                210549.390625,
                169972.390625,
                140521.234375,
                116637.953125,
                117197.625,
                84652.1171875,
                117615.578125,
                103500.921875,
                89320.9453125,
                76313.9296875,
                55575.00390625,
                76784.6796875,
                28829.162109375,
                26051.6171875,
                42957.18359375,
                50342.6953125,
                37611.33984375,
                38202.83203125,
            ],
        }
    ]
    assert result[2] == [
        {
            "mz": [252.1087036133, 252.1572875977, 252.1090698242, 252.1557617188],
            "i": [93112.0859375, 7624.11328125, 76976.7265625, 6090.6440429688],
            "rt": [2.2203779221, 2.2203779221, 2.3452186584, 2.3452186584],
            "polarity": [1.0, 1.0, 1.0, 1.0],
            "precursor_MZ": [252.10887146, 252.10887146, 252.1089477539, 252.1089477539],
            "precursor_intensity": [2872807.5, 2872807.5, 3046732.75, 3046732.75],
            "collision_energy": [23.3333339691, 23.3333339691, 23.3333339691, 23.3333339691],
        }
    ]


def test_get_data_for_atlas_df_and_file(lcmsrun, group, atlas_df, atlas, username):
    # pylint: disable=line-too-long
    result = gdhf.get_data_for_atlas_df_and_file((lcmsrun.hdf5_file, group, atlas_df, atlas))
    expected = (
        {
            "atlas_name": f"HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_{username}_0_0",
            "atlas_unique_id": "749354f7ad974b288624dad533dcbeec",
            "lcmsrun": "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5",
            "group": {
                "creation_time": "2021-05-04T09:41:17",
                "description": "No description",
                "head_id": "61041d07b5a24ca5b88efbda8f319654",
                "items": [
                    {
                        "acquisition_time": 1604770080,
                        "creation_time": "2020-11-13T15:58:43",
                        "description": "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 "
                        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML",
                        "experiment": "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583",
                        "hdf5_file": "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5",
                        "head_id": "7ce51039cfca4426b4e51999ac45d018",
                        "injection_volume": 0.0,
                        "injection_volume_units": "uL",
                        "last_modified": "2021-08-16T12:04:52",
                        "method": None,
                        "mzml_file": "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML",
                        "name": "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML",
                        "pass_qc": False,
                        "prev_uid": "beec4ed4d4b94190aa0776718b5f4fcb",
                        "sample": None,
                        "unique_id": "7ce51039cfca4426b4e51999ac45d018",
                        "username": username,
                    }
                ],
                "last_modified": "2021-05-04T09:41:17",
                "name": f"20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_{username}_0_0_Cone-S1",
                "prev_uid": "origin",
                "short_name": "POS_Cone-S1",
                "unique_id": "61041d07b5a24ca5b88efbda8f319654",
                "username": username,
            },
            "identification": {
                "compound": [
                    {
                        "chebi_id": "CHEBI:17256",
                        "chebi_url": "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17256",
                        "creation_time": "2016-06-17T18:13:15",
                        "description": "A purine 2'-deoxyribonucleoside having adenine as the " "nucleobase.",
                        "formula": "C10H13N5O3",
                        "head_id": "60cd6743e56545c6a6cb066ec3553450",
                        "hmdb_id": "HMDB00101",
                        "hmdb_url": "http://www.hmdb.ca/metabolites/HMDB00101",
                        "img_abc_id": "",
                        "inchi": "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1",
                        "inchi_key": "OLXZPDWKRNYJJZ-RRKCRQDMSA-N",
                        "iupac_name": "",
                        "kegg_id": "C00559",
                        "kegg_url": "http://www.genome.jp/dbget-bin/www_bget?C00559",
                        "last_modified": "2021-08-16T12:04:52",
                        "lipidmaps_id": "",
                        "lipidmaps_url": "",
                        "metacyc_id": "DEOXYADENOSINE",
                        "mono_isotopic_molecular_weight": 251.101839276,
                        "name": "2'-deoxyadenosine",
                        "neutralized_2d_inchi": "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)",
                        "neutralized_2d_inchi_key": "OLXZPDWKRNYJJZ-UHFFFAOYSA-N",
                        "neutralized_inchi": "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1",
                        "neutralized_inchi_key": "OLXZPDWKRNYJJZ-RRKCRQDMSA-N",
                        "num_free_radicals": 0,
                        "number_components": 1,
                        "permanent_charge": 0,
                        "prev_uid": "35b9603309da430da318a0ab8c180457",
                        "pubchem_compound_id": "13730",
                        "pubchem_url": "http://pubchem.ncbi.nlm.nih.gov/compound/13730",
                        "source": "gnps///chebi///metacyc///hmdb",
                        "synonyms": "2'-deoxyadenosine",
                        "unique_id": "60cd6743e56545c6a6cb066ec3553450",
                        "username": username,
                        "wikipedia_url": "",
                    }
                ],
                "creation_time": "2021-02-10T16:20:49",
                "description": "No description",
                "do_normalization": False,
                "frag_references": [],
                "head_id": "18737c7141cc4efaa4545bead13ac751",
                "identification_grade": None,
                "identification_notes": "my id note",
                "intensity_references": [],
                "internal_standard_id": "",
                "internal_standard_to_use": "",
                "last_modified": "2021-08-16T12:04:52",
                "ms1_notes": "keep",
                "ms2_notes": "-1,bad match to ref",
                "mz_references": [
                    {
                        "adduct": "[M+H]+",
                        "creation_time": "2021-02-10T16:20:50",
                        "description": "No description",
                        "detected_polarity": "positive",
                        "enabled": True,
                        "head_id": "eb6d03c9ef574051b92dad7b2fc259a2",
                        "last_modified": "2021-08-16T12:04:52",
                        "lcms_run": None,
                        "modification": "",
                        "mz": 252.1091393,
                        "mz_tolerance": 20.0,
                        "mz_tolerance_units": "ppm",
                        "name": "Untitled",
                        "observed_formula": "",
                        "prev_uid": "18c3020ee0a8422cac1664e3b7d564c6",
                        "ref_type": "",
                        "unique_id": "eb6d03c9ef574051b92dad7b2fc259a2",
                        "username": username,
                    }
                ],
                "name": "2'-deoxyadenosine",
                "prev_uid": "6d0d6618769b4efea1d1457fae2b3358",
                "rt_references": [
                    {
                        "creation_time": "2021-02-10T16:20:50",
                        "description": "No description",
                        "enabled": True,
                        "head_id": "a845ddfdf8ef4713bcef3bdb84999030",
                        "last_modified": "2021-08-16T12:04:52",
                        "lcms_run": None,
                        "name": "Untitled",
                        "prev_uid": "194498f6d0ac475c943130d51c3c036a",
                        "ref_type": "",
                        "rt_max": 2.6964640053707174,
                        "rt_min": 1.6964640053707174,
                        "rt_peak": 2.1964640053707174,
                        "rt_units": "min",
                        "unique_id": "a845ddfdf8ef4713bcef3bdb84999030",
                        "username": username,
                    }
                ],
                "unique_id": "18737c7141cc4efaa4545bead13ac751",
                "username": username,
            },
            "data": {
                "msms": {
                    "data": {
                        "mz": np.array([252.10870361, 252.1572876, 252.10906982, 252.15576172]),
                        "i": np.array([93112.0859375, 7624.11328125, 76976.7265625, 6090.64404297]),
                        "rt": np.array([2.22037792, 2.22037792, 2.34521866, 2.34521866]),
                        "polarity": np.array([1.0, 1.0, 1.0, 1.0]),
                        "precursor_MZ": np.array([252.10887146, 252.10887146, 252.10894775, 252.10894775]),
                        "precursor_intensity": np.array([2872807.5, 2872807.5, 3046732.75, 3046732.75]),
                        "collision_energy": np.array([23.33333397, 23.33333397, 23.33333397, 23.33333397]),
                    }
                },
                "eic": {
                    "mz": [
                        252.1089324951,
                        252.1090087891,
                        252.1088104248,
                        252.1090087891,
                        252.10887146,
                        252.1089324951,
                        252.1089324951,
                        252.1088256836,
                        252.1088867188,
                        252.1090393066,
                        252.1089782715,
                        252.1089630127,
                        252.1089630127,
                        252.1089782715,
                        252.1090240479,
                        252.1089782715,
                        252.1090240479,
                        252.1089324951,
                        252.1090393066,
                        252.1088867188,
                        252.10887146,
                        252.1089324951,
                        252.1089630127,
                        252.1089935303,
                        252.1089172363,
                        252.1089477539,
                        252.1090545654,
                        252.1089630127,
                        252.1090240479,
                        252.1090087891,
                        252.1090393066,
                        252.1090240479,
                        252.1089935303,
                        252.1090240479,
                        252.1089630127,
                        252.1090087891,
                        252.1090240479,
                        252.1089172363,
                        252.1089019775,
                        252.1089477539,
                        252.1089324951,
                        252.1089477539,
                        252.1089477539,
                        252.1089477539,
                        252.1089782715,
                        252.1088867188,
                        252.1089172363,
                        252.1089324951,
                        252.1089782715,
                        252.1089477539,
                        252.1089172363,
                        252.1089324951,
                        252.1089630127,
                        252.1088867188,
                        252.1089630127,
                        252.1085205078,
                        252.1090545654,
                        252.1089935303,
                        252.1088104248,
                        252.1086578369,
                        252.1089935303,
                        252.1085510254,
                        252.1082763672,
                        252.1082458496,
                        252.1084136963,
                        252.1092224121,
                        252.1091766357,
                        252.1092834473,
                        252.1087493896,
                        252.1112518311,
                        252.1088409424,
                        252.1086425781,
                        252.1091766357,
                        252.1094055176,
                    ],
                    "rt": [
                        2.1030805111,
                        2.1084616184,
                        2.1139531136,
                        2.1193552017,
                        2.1248509884,
                        2.1302509308,
                        2.135682106,
                        2.1411821842,
                        2.1459801197,
                        2.1513926983,
                        2.1568279266,
                        2.1622362137,
                        2.1676549911,
                        2.1730883121,
                        2.179015398,
                        2.1845297813,
                        2.1900422573,
                        2.1949694157,
                        2.20002985,
                        2.2055358887,
                        2.2110378742,
                        2.2165191174,
                        2.2219588757,
                        2.2273921967,
                        2.2328462601,
                        2.2382712364,
                        2.2437169552,
                        2.2492566109,
                        2.2547125816,
                        2.2601687908,
                        2.2656960487,
                        2.2704958916,
                        2.2758042812,
                        2.2813498974,
                        2.2868082523,
                        2.2922415733,
                        2.2976748943,
                        2.3031060696,
                        2.308131218,
                        2.313628912,
                        2.3185498714,
                        2.3239560127,
                        2.3293914795,
                        2.3349123001,
                        2.3403663635,
                        2.346799612,
                        2.3522267342,
                        2.3576600552,
                        2.3631224632,
                        2.3685662746,
                        2.3740911484,
                        2.3794057369,
                        2.3848536015,
                        2.3903660774,
                        2.3953785896,
                        2.4006638527,
                        2.4062638283,
                        2.411709547,
                        2.4171659946,
                        2.4226117134,
                        2.4302260876,
                        2.4357616901,
                        2.4407405853,
                        2.4461927414,
                        2.451615572,
                        2.4571509361,
                        2.4627010822,
                        2.4681572914,
                        2.4735822678,
                        2.4735822678,
                        2.4787945747,
                        2.4842174053,
                        2.4896612167,
                        2.495146513,
                    ],
                    "intensity": [
                        312203.5,
                        387914.59375,
                        308308.5,
                        334653.59375,
                        339521.625,
                        345527.21875,
                        292437.34375,
                        413614.53125,
                        300285.28125,
                        383848.71875,
                        404313.21875,
                        377231.34375,
                        453965.5625,
                        431327.0,
                        523180.0625,
                        510239.8125,
                        631459.1875,
                        807419.5,
                        842647.5625,
                        1053031.625,
                        1082361.625,
                        1198966.625,
                        1109162.375,
                        1126347.125,
                        1373071.5,
                        1589018.375,
                        1281309.875,
                        1660166.75,
                        1492912.25,
                        2029801.5,
                        2029874.125,
                        2035966.625,
                        2010867.875,
                        2036981.375,
                        2148879.25,
                        2359861.25,
                        2054066.125,
                        1691976.0,
                        1778159.125,
                        1776166.125,
                        1752154.125,
                        1575676.875,
                        1199910.625,
                        1259708.25,
                        1087384.375,
                        826077.125,
                        802296.875,
                        547785.125,
                        545340.0625,
                        584624.4375,
                        468524.8125,
                        305931.1875,
                        330310.34375,
                        309740.625,
                        289212.71875,
                        230440.9375,
                        210549.390625,
                        169972.390625,
                        140521.234375,
                        116637.953125,
                        117197.625,
                        84652.1171875,
                        117615.578125,
                        103500.921875,
                        89320.9453125,
                        76313.9296875,
                        55575.00390625,
                        76784.6796875,
                        28829.162109375,
                        26051.6171875,
                        42957.18359375,
                        50342.6953125,
                        37611.33984375,
                        38202.83203125,
                    ],
                },
                "ms1_summary": {
                    "num_ms1_datapoints": 74.0,
                    "mz_peak": 252.1090087891,
                    "rt_peak": 2.2922415733,
                    "mz_centroid": 252.10896296693303,
                    "rt_centroid": 2.2720730579808084,
                    "peak_height": 2359861.25,
                    "peak_area": 57016800.755859375,
                },
            },
        },
    )
    assert len(result) == len(expected)
    assert result[0].keys() == expected[0].keys()
    for key in ["atlas_name", "lcmsrun", "data"]:
        assert str(result[0][key]) == str(expected[0][key])
