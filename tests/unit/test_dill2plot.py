# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long

import pandas as pd

from metatlas.plots import dill2plots


def test_rt_range_overlaps_none():
    i = type("", (), {})()
    i.rt_min = 10
    i.rt_max = 20
    j = type("", (), {})()
    j.rt_min = 30
    j.rt_max = 40
    assert not dill2plots.rt_range_overlaps(i, j)
    assert not dill2plots.rt_range_overlaps(j, i)


def test_rt_range_overlaps_partial():
    i = type("", (), {})()
    i.rt_min = 10
    i.rt_max = 20
    j = type("", (), {})()
    j.rt_min = 15
    j.rt_max = 25
    assert dill2plots.rt_range_overlaps(i, j)
    assert dill2plots.rt_range_overlaps(j, i)


def test_rt_range_overlaps_inside():
    i = type("", (), {})()
    i.rt_min = 10
    i.rt_max = 20
    j = type("", (), {})()
    j.rt_min = 12
    j.rt_max = 18
    assert dill2plots.rt_range_overlaps(i, j)
    assert dill2plots.rt_range_overlaps(j, i)


def test_within_tolerance_yes():
    assert dill2plots.within_tolerance(99, 100, 0.02)


def test_within_tolerance_no():
    assert not dill2plots.within_tolerance(99, 100, 0.002)


def test_filter_metatlas_objects_by_list_remove_all():
    i = type("", (), {})()
    i.myattr = [1, 2]
    j = type("", (), {})()
    j.myattr = [2, 3]
    assert [] == dill2plots.filter_metatlas_objects_by_list([i, j], "myattr", [4, 5])


def test_filter_metatlas_objects_by_list_remove_none():
    i = type("", (), {})()
    i.myattr = [1, 2]
    j = type("", (), {})()
    j.myattr = [2, 3]
    assert [i, j] == dill2plots.filter_metatlas_objects_by_list([i, j], "myattr", [2])


def test_remove_metatlas_objects_by_list_remove_none():
    i = type("", (), {})()
    i.myattr = [1, 2]
    j = type("", (), {})()
    j.myattr = [2, 3]
    assert [i, j] == dill2plots.remove_metatlas_objects_by_list([i, j], "myattr", [4, 5])


def test_remove_metatlas_objects_by_list_remove_all():
    i = type("", (), {})()
    i.myattr = [1, 2]
    j = type("", (), {})()
    j.myattr = [2, 3]
    assert [] == dill2plots.remove_metatlas_objects_by_list([i, j], "myattr", [0, 2, 5])


def test_export_atlas_to_spreadsheet(atlas, username):
    # pylint: disable=line-too-long
    expected = (
        """{"chebi_id":{"0":"CHEBI:17256"},"chebi_url":{"0":"http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17256"},"creation_time":{"0":1466212395.0},"description":{"0":"A purine 2'-deoxyribonucleoside having adenine as the nucleobase."},"formula":{"0":"C10H13N5O3"},"head_id":{"0":"60cd6743e56545c6a6cb066ec3553450"},"hmdb_id":{"0":"HMDB00101"},"hmdb_url":{"0":"http://www.hmdb.ca/metabolites/HMDB00101"},"img_abc_id":{"0":""},"inchi":{"0":"InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"},"inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"iupac_name":{"0":""},"kegg_id":{"0":"C00559"},"kegg_url":{"0":"http://www.genome.jp/dbget-bin/www_bget?C00559"},"last_modified":{"0":1612996604.0},"lipidmaps_id":{"0":""},"lipidmaps_url":{"0":""},"metacyc_id":{"0":"DEOXYADENOSINE"},"mono_isotopic_molecular_weight":{"0":251.101839276},"name":{"0":"2'-deoxyadenosine"},"neutralized_2d_inchi":{"0":"InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)"},"neutralized_2d_inchi_key":{"0":"OLXZPDWKRNYJJZ-UHFFFAOYSA-N"},"neutralized_inchi":{"0":"InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"},"neutralized_inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"num_free_radicals":{"0":0.0},"number_components":{"0":1.0},"permanent_charge":{"0":0.0},"prev_uid":{"0":"origin"},"pubchem_compound_id":{"0":"13730"},"pubchem_url":{"0":"http://pubchem.ncbi.nlm.nih.gov/compound/13730"},"source":{"0":"gnps///chebi///metacyc///hmdb"},"synonyms":{"0":"2'-deoxyadenosine"},"unique_id":{"0":"60cd6743e56545c6a6cb066ec3553450"},"username":{"0":"""  # noqa:  E501
        f'"{username}"'
        """},"wikipedia_url":{"0":""},"label":{"0":"2'-deoxyadenosine"},"id_notes":{"0":"No description"},"ms1_notes":{"0":"keep"},"ms2_notes":{"0":"-1,bad match to ref"},"identification_notes":{"0":"my id note"},"rt_min":{"0":1.6964640054},"rt_max":{"0":2.6964640054},"rt_peak":{"0":2.1964640054},"mz":{"0":252.1091393},"mz_tolerance":{"0":20.0},"adduct":{"0":"[M+H]+"},"polarity":{"0":"positive"}}"""
    )  # noqa:  E501
    assert expected == dill2plots.export_atlas_to_spreadsheet(atlas).to_json().replace(r"\/", "/")


def test_filter_atlas01(metatlas_dataset):
    assert len(dill2plots.filter_atlas(metatlas_dataset.atlas_df, metatlas_dataset, 1, 2.30e6)) == 1
    assert len(dill2plots.filter_atlas(metatlas_dataset.atlas_df, metatlas_dataset, 1, 2.36e6)) == 0
    assert len(dill2plots.filter_atlas(metatlas_dataset.atlas_df, metatlas_dataset, 73, 1e4)) == 1
    assert len(dill2plots.filter_atlas(metatlas_dataset.atlas_df, metatlas_dataset, 74, 1e4)) == 0


def test_strong_signal_compound_idxs(metatlas_dataset):
    assert dill2plots.strong_signal_compound_idxs(metatlas_dataset, 1, 2.30e6) == [0]
    assert dill2plots.strong_signal_compound_idxs(metatlas_dataset, 1, 2.36e6) == []
    assert dill2plots.strong_signal_compound_idxs(metatlas_dataset, 73, 1e4) == [0]
    assert dill2plots.strong_signal_compound_idxs(metatlas_dataset, 74, 1e4) == []


def test_get_msms_hits01(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs)
    hits = dill2plots.get_msms_hits(metatlas_dataset)
    expected = (
        """{"score":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":0.0,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":0.0},"num_matches":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":1,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":1},"msv_query_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1087036133,null,252.1572875977],[null,null,null,null,null,null,null,null,null,null,null,null,null,93112.0859375,null,7624.11328125]],"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1090698242,null,252.1557617188],[null,null,null,null,null,null,null,null,null,null,null,null,null,76976.7265625,null,6090.6440429688]]},"msv_ref_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]],"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]]},"name":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', """
        """2.2203779221)":"2\'-deoxyadenosine","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"2\'-deoxyadenosine"},"adduct":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"[M+H]+","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"[M+H]+"},"inchi_key":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.1091393,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1091393},"measured_precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.10887146,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1089477539},"measured_precursor_intensity":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":2872807.5,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":3046732.75}}"""
    )
    assert len(hits) == 2
    assert expected == hits.to_json()


def test_get_msms_hits02(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs)
    metatlas_dataset.set_rt(0, "rt_max", 2.3)  # reduce the RT bounds to so that only one hit falls within
    hits_small = dill2plots.get_msms_hits(metatlas_dataset)
    expected_small = """{"score":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":0.0},"num_matches":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":1},"msv_query_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1087036133,null,252.1572875977],[null,null,null,null,null,null,null,null,null,null,null,null,null,93112.0859375,null,7624.11328125]]},"msv_ref_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]]},"name":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"2\'-deoxyadenosine"},"adduct":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"[M+H]+"},"inchi_key":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.1091393},"measured_precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.10887146},"measured_precursor_intensity":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":2872807.5}}"""
    assert len(hits_small) == 1
    assert hits_small.to_json() == expected_small


def test_get_msms_hits03(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs)
    metatlas_dataset.set_rt(0, "rt_max", 2.3)  # reduce the RT bounds to so that only one hit falls within
    hits_large = dill2plots.get_msms_hits(
        metatlas_dataset, extra_time=0.75
    )  # with extra_time the second hits is also included
    expected_large = (
        """{"score":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":0.0,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":0.0},"num_matches":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":1,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":1},"msv_query_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1087036133,null,252.1572875977],[null,null,null,null,null,null,null,null,null,null,null,null,null,93112.0859375,null,7624.11328125]],"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1090698242,null,252.1557617188],[null,null,null,null,null,null,null,null,null,null,null,null,null,76976.7265625,null,6090.6440429688]]},"msv_ref_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]],"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]]},"name":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', """
        """2.2203779221)":"2\'-deoxyadenosine","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"2\'-deoxyadenosine"},"adduct":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"[M+H]+","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"[M+H]+"},"inchi_key":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.1091393,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1091393},"measured_precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.10887146,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1089477539},"measured_precursor_intensity":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":2872807.5,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":3046732.75}}"""
    )
    assert len(hits_large) == 2
    assert hits_large.to_json() == expected_large


def test_get_msms_hits04(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs)
    metatlas_dataset.set_rt(0, "rt_max", 2.20)  # reduce the RT bounds to so that no hits fall within
    hits = dill2plots.get_msms_hits(metatlas_dataset)
    assert len(hits) == 0


def test_get_msms_hits05(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs)
    expected = (
        """{"score":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":0.0,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":0.0},"num_matches":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":1,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":1},"msv_query_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1087036133,null,252.1572875977],[null,null,null,null,null,null,null,null,null,null,null,null,null,93112.0859375,null,7624.11328125]],"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[null,null,null,null,null,null,null,null,null,null,null,null,null,252.1090698242,null,252.1557617188],[null,null,null,null,null,null,null,null,null,null,null,null,null,76976.7265625,null,6090.6440429688]]},"msv_ref_aligned":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]],"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[57.0345,63.3177,63.3205,69.0344,71.0499,73.0292,84.9778,99.0447,117.055,118.059,136.062,137.066,236.709,252.109,253.112,null],[176328.0,328818.0,274432.0,197637.0,896360.0,1192020.0,378547.0,3921880.0,15737700.0,266131.0,144220000.0,3455270.0,185227.0,20960800.0,1284450.0,null]]},"name":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', """
        """2.2203779221)":"2\'-deoxyadenosine","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"2\'-deoxyadenosine"},"adduct":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"[M+H]+","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"[M+H]+"},"inchi_key":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.1091393,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1091393},"measured_precursor_mz":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.10887146,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1089477539},"measured_precursor_intensity":{"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":2872807.5,"(\'metatlas\', \'c7dddd297e104ca79caea72a90150532\', \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":3046732.75}}"""
    )
    hits = dill2plots.get_msms_hits(metatlas_dataset, keep_nonmatches=True)
    assert len(hits) == 2
    assert hits.to_json() == expected


def test_get_msms_hits06(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs.iloc[0:0])
    expected = (
        """{"score":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":2872807.5,"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":3046732.75},"num_matches":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":0,"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":0},"msv_query_aligned":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[252.1087036133,252.1572875977],[93112.0859375,7624.11328125]],"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[252.1090698242,252.1557617188],[76976.7265625,6090.6440429688]]},"msv_ref_aligned":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":[[null,null],[null,null]],"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":[[null,null],[null,null]]},"name":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"2\'-deoxyadenosine","(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"2\'-deoxyadenosine"},"adduct":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"[M+H]+","(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"[M+H]+"},"inchi_key":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N","(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"precursor_mz":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221"""
        """)":252.1091393,"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1091393},"measured_precursor_mz":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":252.10887146,"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":252.1089477539},"measured_precursor_intensity":{"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.2203779221)":2872807.5,"(nan, nan, \'20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\', 2.3452186584)":3046732.75}}"""
    )
    hits = dill2plots.get_msms_hits(metatlas_dataset, keep_nonmatches=True)
    assert len(hits) == 2
    assert hits.to_json() == expected


def test_get_msms_hits07(metatlas_dataset, msms_refs, mocker):
    mocker.patch("metatlas.plots.dill2plots.get_refs", return_value=msms_refs.iloc[0:0])
    hits = dill2plots.get_msms_hits(metatlas_dataset)
    assert len(hits) == 0


def test_instructions01(instructions, mocker):
    mocker.patch("pandas.read_csv", return_value=instructions)
    inst = dill2plots.InstructionSet("fake_path")
    assert inst.query("FAKE_INCHI_KEY", "", "", "") == ["No instructions for this data"]
    out1 = inst.query("OIRDTQYFTABQOQ-KQYNXXCUSA-N", "", "", "")
    assert len(out1) == 1
    assert out1 == ["Note 4 is column and polarity independent"]
    assert inst.query("OLXZPDWKRNYJJZ-RRKCRQDMSA-N", "", "HILICZ", "") == ["No instructions for this data"]
    out2 = inst.query("OLXZPDWKRNYJJZ-RRKCRQDMSA-N", "", "", "")
    assert out2 == ["Note 5 has a fake column and should not match"]
    out3 = inst.query("HXACOUQIXZGNBF-UHFFFAOYSA-N", "[M+H]+", "HILICZ", "positive")
    assert out3 == ["Note 2 contain a comma, right?"]
    out4 = inst.query("OIRDTQYFTABQOQ-KQYNXXCUSA-N", "", "C18", "negative")
    assert out4 == ["Note 4 is column and polarity independent"]


def test_make_atlas_from_spreadsheet(mocker):
    csv_data = pd.DataFrame({"rt_min": [1.1], "rt_peak": [1.3], "rt_max": [1.5], "mz": [234.6578]})
    mocker.patch("metatlas.plots.dill2plots._get_dataframe", return_value=csv_data)
    atlas = dill2plots.make_atlas_from_spreadsheet(
        "foo.csv", "test_atlas_99", "csv", polarity="positive", mz_tolerance=5
    )
    assert len(atlas.compound_identifications) == 1
    assert atlas.compound_identifications[0].rt_references[0].rt_peak == 1.3
