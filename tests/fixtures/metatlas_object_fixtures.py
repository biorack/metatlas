# /pylint: disable=line-too-long, missing-function-docstring, missing-module-docstring

import pandas as pd
import pytest
from metatlas.datastructures import metatlas_objects as metob


@pytest.fixture(name="compound")
def fixture_compound():
    compound = metob.Compound()
    compound.unique_id = "60cd6743e56545c6a6cb066ec3553450"
    compound.mono_isotopic_molecular_weight = 251.101839276
    compound.creation_time = 1466212395
    compound.synonyms = "2'-deoxyadenosine"  # value was pruned down
    compound.inchi_key = "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
    compound.chebi_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17256"
    compound.permanent_charge = 0
    compound.img_abc_id = ""
    compound.neutralized_2d_inchi = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)"  # noqa: E501
    compound.lipidmaps_url = ""
    compound.source = "gnps///chebi///metacyc///hmdb"
    compound.kegg_url = "http://www.genome.jp/dbget-bin/www_bget?C00559"
    compound.hmdb_url = "http://www.hmdb.ca/metabolites/HMDB00101"
    compound.wikipedia_url = ""
    compound.head_id = "60cd6743e56545c6a6cb066ec3553450"
    compound.formula = "C10H13N5O3"
    compound.number_components = 1
    compound.iupac_name = ""
    compound.username = "wjholtz"
    compound.pubchem_compound_id = "13730"
    compound.description = "A purine 2'-deoxyribonucleoside having adenine as the nucleobase."
    compound.metacyc_id = "DEOXYADENOSINE"
    compound.kegg_id = "C00559"
    compound.hmdb_id = "HMDB00101"
    compound.chebi_id = "CHEBI:17256"
    compound.inchi = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"  # noqa: E501
    compound.neutralized_inchi_key = "OLXZPDWKRNYJJZ-RRKCRQDMSA-N"
    compound.prev_uid = "origin"
    compound.neutralized_inchi = "InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"  # noqa: E501
    compound.name = "2'-deoxyadenosine"
    compound.neutralized_2d_inchi_key = "OLXZPDWKRNYJJZ-UHFFFAOYSA-N"
    compound.num_free_radicals = 0
    compound.lipidmaps_id = ""
    compound.last_modified = 1612996604
    compound.pubchem_url = "http://pubchem.ncbi.nlm.nih.gov/compound/13730"
    return compound


@pytest.fixture(name="rt_reference")
def fixture_rt_reference():
    rt_ref = metob.RtReference()
    rt_ref.unique_id = "a845ddfdf8ef4713bcef3bdb84999030"
    rt_ref.username = "wjholtz"
    rt_ref.rt_units = "min"
    rt_ref.description = "No description"
    rt_ref.rt_peak = "2.1964640053707174"
    rt_ref.enabled = True
    rt_ref.creation_time = 1613002850
    rt_ref.lcms_run = None
    rt_ref.rt_min = 1.6964640053707174
    rt_ref.last_modified = 1613002979
    rt_ref.ref_type = ""
    rt_ref.prev_uid = "origin"
    rt_ref.rt_max = 2.6964640053707174
    rt_ref.name = "Untitled"
    rt_ref.head_id = "a845ddfdf8ef4713bcef3bdb84999030"
    return rt_ref


@pytest.fixture(name="mz_reference")
def fixture_mz_reference():
    mz_ref = metob.MzReference()
    mz_ref.unique_id = "eb6d03c9ef574051b92dad7b2fc259a2"
    mz_ref.username = "wjholtz"
    mz_ref.adduct = "[M+H]+"
    mz_ref.description = "No description"
    mz_ref.mz_tolerance_units = "ppm"
    mz_ref.enabled = True
    mz_ref.mz = 252.1091393
    mz_ref.creation_time = 1613002850
    mz_ref.lcms_run = None
    mz_ref.mz_tolerance = 20.0
    mz_ref.last_modified = 1613002979
    mz_ref.detected_polarity = "positive"
    mz_ref.modification = ""
    mz_ref.ref_type = ""
    mz_ref.observed_formula = ""
    mz_ref.prev_uid = "origin"
    mz_ref.name = "Untitled"
    mz_ref.head_id = "eb6d03c9ef574051b92dad7b2fc259a2"
    return mz_ref


@pytest.fixture(name="compound_identification")
def fixture_compound_identification(compound, rt_reference, mz_reference):
    ident = metob.CompoundIdentification()
    ident.unique_id = "18737c7141cc4efaa4545bead13ac751"
    ident.username = "wjholtz"
    ident.description = "No description"
    ident.creation_time = 1613002849
    ident.last_modified = 1613002979
    ident.identification_grade = None
    ident.prev_uid = "origin"
    ident.name = "2'-deoxyadenosine"
    ident.head_id = "18737c7141cc4efaa4545bead13ac751"
    ident.internal_standard_to_use = ""
    ident.internal_standard_id = ""
    ident.do_normalization = False
    ident.identification_notes = "my id note"
    ident.ms2_notes = "bad match to ref"
    ident.ms1_notes = "keep"
    ident.frag_references = []
    ident.intensity_references = []
    ident.compound = [compound]
    ident.mz_references = [mz_reference]
    ident.rt_references = [rt_reference]
    return ident


@pytest.fixture(name="atlas")
def fixture_atlas(compound_identification):
    small_atlas = metob.Atlas()
    small_atlas.compound_identifications = [compound_identification]
    return small_atlas


@pytest.fixture(name="compound_2")
def fixture_compound_2():
    compound = metob.Compound()
    compound.chebi_id = "CHEBI:16335"
    compound.chebi_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16335"
    compound.creation_time = 1466212384
    compound.description = "A ribonucleoside composed of a molecule of adenine attached to a ribofuranose moiety via a beta1N9-glycosidic bond."
    compound.formula = "C10H13N5O4"
    compound.head_id = "1ad02275f47b4033a451e99874f4764f"
    compound.hmdb_id = "HMDB00050"
    compound.hmdb_url = "http://www.hmdb.ca/metabolites/HMDB00050"
    compound.img_abc_id = ""
    compound.inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)/t4-,6-,7-,10-/m1/s1"
    compound.inchi_key = "OIRDTQYFTABQOQ-KQYNXXCUSA-N"
    compound.iupac_name = ""
    compound.kegg_id = "C00212"
    compound.kegg_url = "http://www.genome.jp/dbget-bin/www_bget?C00212"
    compound.last_modified = 1612996604
    compound.lipidmaps_id = ""
    compound.lipidmaps_url = ""
    compound.metacyc_id = "ADENOSINE"
    compound.mono_isotopic_molecular_weight = 267.096753896
    compound.name = "adenosine"
    compound.neutralized_2d_inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)"
    compound.neutralized_2d_inchi_key = "OIRDTQYFTABQOQ-UHFFFAOYSA-N"
    compound.neutralized_inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)/t4-,6-,7-,10-/m1/s1"
    compound.neutralized_inchi_key = "OIRDTQYFTABQOQ-KQYNXXCUSA-N"
    compound.num_free_radicals = 0
    compound.number_components = 1
    compound.permanent_charge = 0
    compound.prev_uid = "origin"
    compound.pubchem_compound_id = "60961"
    compound.pubchem_url = "http://pubchem.ncbi.nlm.nih.gov/compound/60961"
    compound.source = "chebi///wikidata///metacyc///gnps///hmdb"
    compound.synonyms = "adenosine///58-61-7///Adenocard///Adenoscan"  # this value was pruned down
    compound.unique_id = "1ad02275f47b4033a451e99874f4764f"
    compound.username = "wjholtz"
    compound.wikipedia_url = ""
    return compound


@pytest.fixture(name="rt_reference_2")
def fixture_rt_reference_2():
    rt_ref = metob.RtReference()
    rt_ref.creation_time = 1613002857
    rt_ref.description = "No description"
    rt_ref.enabled = True
    rt_ref.head_id = "f74622bcef924f5390ba6e127633e731"
    rt_ref.last_modified = 1613002980
    rt_ref.lcms_run = None
    rt_ref.name = "Untitled"
    rt_ref.prev_uid = "origin"
    rt_ref.ref_type = ""
    rt_ref.rt_max = 3.5233184079926665
    rt_ref.rt_min = 2.5233184079926665
    rt_ref.rt_peak = 3.0233184079926665
    rt_ref.rt_units = "min"
    rt_ref.unique_id = "f74622bcef924f5390ba6e127633e731"
    rt_ref.username = "wjholtz"
    return rt_ref


@pytest.fixture(name="mz_reference_2")
def fixture_mz_reference_2():
    mz_ref = metob.MzReference()
    mz_ref.adduct = "[M+H]+"
    mz_ref.creation_time = 1613002857
    mz_ref.description = "No description"
    mz_ref.detected_polarity = "positive"
    mz_ref.enabled = True
    mz_ref.head_id = "b0e3cf0df44a4079be7908c6b525d3ac"
    mz_ref.last_modified = 1613002980
    mz_ref.lcms_run = None
    mz_ref.modification = ""
    mz_ref.mz = 268.1040539
    mz_ref.mz_tolerance = 20.0
    mz_ref.mz_tolerance_units = "ppm"
    mz_ref.name = "Untitled"
    mz_ref.observed_formula = ""
    mz_ref.prev_uid = "origin"
    mz_ref.ref_type = ""
    mz_ref.unique_id = "b0e3cf0df44a4079be7908c6b525d3ac"
    mz_ref.username = "wjholtz"
    return mz_ref


@pytest.fixture(name="compound_identification_2")
def fixture_compound_identification_2(compound_2, rt_reference_2, mz_reference_2):
    ident = metob.CompoundIdentification()
    ident.creation_time = 1613002856
    ident.description = "No description"
    ident.do_normalization = False
    ident.frag_references = []
    ident.head_id = "6cca7aa44c0e4a109f695ba980d69472"
    ident.identification_grade = None
    ident.identification_notes = ""
    ident.intensity_references = []
    ident.internal_standard_id = ""
    ident.internal_standard_to_use = ""
    ident.last_modified = 1613002980
    ident.ms1_notes = ""
    ident.ms2_notes = ""
    ident.name = "adenosine"
    ident.prev_uid = "origin"
    ident.unique_id = "6cca7aa44c0e4a109f695ba980d69472"
    ident.username = "wjholtz"
    ident.frag_references = []
    ident.intensity_references = []
    ident.compound = [compound_2]
    ident.mz_references = [mz_reference_2]
    ident.rt_references = [rt_reference_2]
    return ident


@pytest.fixture(name="atlas_with_2_cids")
def fixture_atlas_with_2_cids(compound_identification, compound_identification_2):
    small_atlas = metob.Atlas()
    small_atlas.compound_identifications = [
        compound_identification,
        compound_identification_2,
    ]
    return small_atlas


@pytest.fixture(name="lcmsrun")
def fixture_lcmsrun():
    run = metob.LcmsRun()
    run.unique_id = "7ce51039cfca4426b4e51999ac45d018"
    run.username = "root"
    run.hdf5_file = "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5"  # noqa: E501
    run.description = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML"  # noqa: E501
    run.creation_time = 1605311923
    run.sample = None
    run.last_modified = 1620101765
    run.mzml_file = "/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML"  # noqa: E501
    run.prev_uid = "28323058b6e84a9db0f9e802544764e3"
    run.method = None
    run.name = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.mzML"  # noqa: E501
    run.head_id = "7ce51039cfca4426b4e51999ac45d018"
    run.experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    run.injection_volume = 0.0
    run.injection_volume_units = "uL"
    run.acquisition_time = 1604770080
    run.pass_qc = False
    return run


@pytest.fixture(name="group")
def fixture_group(lcmsrun):
    grp = metob.Group()
    grp.items = [lcmsrun]
    grp.unique_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.username = "root"
    grp.description = "No description"
    grp.creation_time = 1620146477
    grp.last_modified = 1620146477
    grp.prev_uid = "origin"
    grp.name = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root0_Cone-S1"
    grp.head_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.short_name = "POS_Cone-S1"
    return grp


@pytest.fixture(name="group_with_2_lcmsruns")
def fixture_group_with_2_lcmsruns(lcmsrun):
    grp = metob.Group()
    grp.items = [lcmsrun, lcmsrun]
    grp.unique_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.username = "root"
    grp.description = "No description"
    grp.creation_time = 1620146477
    grp.last_modified = 1620146477
    grp.prev_uid = "origin"
    grp.name = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root0_Cone-S1"
    grp.head_id = "61041d07b5a24ca5b88efbda8f319654"
    grp.short_name = "POS_Cone-S1"
    return grp


@pytest.fixture(name="hits")
def fixture_hits():
    hits_plus = pd.DataFrame(
        data={
            "score": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": 0.7253785748,
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": 0.8688691781,
            },
            "num_matches": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": 6,
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": 7,
            },
            "msv_query_aligned": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": [
                    [
                        None,
                        None,
                        None,
                        None,
                        56.7212257385,
                        59.0436058044,
                        71.0422821045,
                        73.0214157104,
                        None,
                        89.1910018921,
                        99.0413742065,
                        104.3592529297,
                        104.3681869507,
                        117.0548171997,
                        None,
                        118.9432754517,
                        136.0619506836,
                        None,
                        None,
                        None,
                        145.9665527344,
                        163.9772491455,
                        169.9678497314,
                        177.1133270264,
                        187.9771575928,
                        205.9878387451,
                        210.9933166504,
                        229.0038452148,
                        252.0215606689,
                        252.1087036133,
                        252.1572875977,
                        252.2064666748,
                    ],
                    [
                        None,
                        None,
                        None,
                        None,
                        3361.7712402344,
                        6589.943359375,
                        6501.9853515625,
                        4987.177734375,
                        None,
                        3257.0708007812,
                        13393.138671875,
                        3280.0544433594,
                        4276.0112304688,
                        57809.1875,
                        None,
                        4965.7436523438,
                        648640.5625,
                        None,
                        None,
                        None,
                        11511.76171875,
                        10362.68359375,
                        5714.70703125,
                        9354.2353515625,
                        73409.0078125,
                        257685.234375,
                        53554.28125,
                        193491.515625,
                        5038.1469726562,
                        93112.0859375,
                        7624.11328125,
                        4599.4125976562,
                    ],
                ],
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": [
                    [
                        None,
                        50.2002449036,
                        55.0126533508,
                        57.0280647278,
                        None,
                        None,
                        None,
                        None,
                        68.2973327637,
                        69.0266494751,
                        73.0213851929,
                        None,
                        74.6972732544,
                        80.862159729,
                        82.4692306519,
                        85.0231246948,
                        87.0394363403,
                        92.4544296265,
                        92.4610061646,
                        104.3785171509,
                        115.0390701294,
                        126.1923675537,
                        133.0496368408,
                        136.0618743896,
                        None,
                        None,
                        None,
                        None,
                        144.5760345459,
                        181.1904449463,
                        230.6756896973,
                        268.1039733887,
                    ],
                    [
                        None,
                        87283.4296875,
                        105163.625,
                        246350.078125,
                        None,
                        None,
                        None,
                        None,
                        81607.3046875,
                        107886.640625,
                        150512.90625,
                        None,
                        99324.7109375,
                        80050.4375,
                        108701.53125,
                        278198.71875,
                        95401.265625,
                        92632.890625,
                        111341.5625,
                        119245.7734375,
                        170358.671875,
                        103961.4296875,
                        226297.9375,
                        48576460.0,
                        None,
                        None,
                        None,
                        None,
                        98098.609375,
                        100016.9296875,
                        119618.1015625,
                        16002674.0,
                    ],
                ],
            },
            "msv_ref_aligned": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": [
                    [
                        57.0345,
                        63.3177,
                        63.3205,
                        69.0344,
                        None,
                        None,
                        71.0499,
                        73.0292,
                        84.9778,
                        None,
                        99.0447,
                        None,
                        None,
                        117.055,
                        118.059,
                        None,
                        136.062,
                        137.066,
                        236.709,
                        253.112,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        252.109,
                        None,
                        None,
                    ],
                    [
                        176328.0,
                        328818.0,
                        274432.0,
                        197637.0,
                        None,
                        None,
                        896360.0,
                        1192020.0,
                        378547.0,
                        None,
                        3921880.0,
                        None,
                        None,
                        15737700.0,
                        266131.0,
                        None,
                        144220000.0,
                        3455270.0,
                        185227.0,
                        1284450.0,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        20960800.0,
                        None,
                        None,
                    ],
                ],
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": [
                    [
                        56.7603,
                        None,
                        None,
                        57.0346,
                        61.0292,
                        61.8182,
                        64.9491,
                        67.9255,
                        None,
                        None,
                        73.0292,
                        82.0663,
                        None,
                        None,
                        None,
                        85.0293,
                        None,
                        None,
                        None,
                        None,
                        115.04,
                        None,
                        133.05,
                        136.062,
                        137.067,
                        183.555,
                        230.198,
                        269.108,
                        None,
                        None,
                        None,
                        268.105,
                    ],
                    [
                        35523.7,
                        None,
                        None,
                        184839.0,
                        43216.2,
                        40066.3,
                        40362.0,
                        41550.6,
                        None,
                        None,
                        93791.1,
                        293258.0,
                        None,
                        None,
                        None,
                        202756.0,
                        None,
                        None,
                        None,
                        None,
                        184050.0,
                        None,
                        364543.0,
                        29646700.0,
                        830130.0,
                        51455.4,
                        51206.7,
                        970064.0,
                        None,
                        None,
                        None,
                        12412800.0,
                    ],
                ],
            },
            "name": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": "2'-deoxyadenosine",
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": "adenosine",
            },
            "adduct": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": "[M+H]+",
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": "[M+H]+",
            },
            "inchi_key": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": "OLXZPDWKRNYJJZ-RRKCRQDMSA-N",
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": "OIRDTQYFTABQOQ-KQYNXXCUSA-N",
            },
            "precursor_mz": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": 252.1091393,
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": 268.1040539,
            },
            "measured_precursor_mz": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": 252.10887146,
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": 268.103729248,
            },
            "measured_precursor_intensity": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": 2872807.5,
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": 75979424.0,
            },
            "copy_index": {
                "('metatlas', 'c7dddd297e104ca79caea72a90150532', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5', 2.2203779220581055)": [
                    "metatlas",
                    "c7dddd297e104ca79caea72a90150532",
                    "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5",
                    2.2203779221,
                ],
                "('metatlas', 'cf5e8df145f64bf0856fbf852d1bdb64', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5', 3.0264527797698975)": [
                    "metatlas",
                    "cf5e8df145f64bf0856fbf852d1bdb64",
                    "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5",
                    3.0264527798,
                ],
            },
        }
    )
    hits_plus.index = pd.MultiIndex.from_tuples(
        hits_plus["copy_index"], names=["database", "id", "file_name", "msms_scan"]
    )
    hits_plus.drop(columns=["copy_index"], inplace=True)
    return hits_plus
