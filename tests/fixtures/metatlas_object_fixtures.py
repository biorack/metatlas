from metatlas.datastructures import metatlas_objects as metob
import pytest


@pytest.fixture
def compound():
    compound = metob.Compound()
    compound.unique_id = '60cd6743e56545c6a6cb066ec3553450'
    compound.mono_isotopic_molecular_weight = 251.101839276
    compound.creation_time = 1466212395
    compound.synonyms = "2'-deoxyadenosine"  # value was pruned down
    compound.inchi_key = 'OLXZPDWKRNYJJZ-RRKCRQDMSA-N'
    compound.chebi_url = 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17256'
    compound.permanent_charge = 0
    compound.img_abc_id = ''
    compound.neutralized_2d_inchi = 'InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)'  # noqa: E501
    compound.lipidmaps_url = ''
    compound.source = 'gnps///chebi///metacyc///hmdb'
    compound.kegg_url = 'http://www.genome.jp/dbget-bin/www_bget?C00559'
    compound.hmdb_url = 'http://www.hmdb.ca/metabolites/HMDB00101'
    compound.wikipedia_url = ''
    compound.head_id = '60cd6743e56545c6a6cb066ec3553450'
    compound.formula = 'C10H13N5O3'
    compound.number_components = 1
    compound.iupac_name = ''
    compound.username = 'wjholtz'
    compound.pubchem_compound_id = '13730'
    compound.description = "A purine 2'-deoxyribonucleoside having adenine as the nucleobase."
    compound.metacyc_id = 'DEOXYADENOSINE'
    compound.kegg_id = 'C00559'
    compound.hmdb_id = 'HMDB00101'
    compound.chebi_id = 'CHEBI:17256'
    compound.inchi = 'InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1'  # noqa: E501
    compound.neutralized_inchi_key = 'OLXZPDWKRNYJJZ-RRKCRQDMSA-N'
    compound.prev_uid = 'origin'
    compound.neutralized_inchi = 'InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1'  # noqa: E501
    compound.name = "2'-deoxyadenosine"
    compound.neutralized_2d_inchi_key = 'OLXZPDWKRNYJJZ-UHFFFAOYSA-N'
    compound.num_free_radicals = 0
    compound.lipidmaps_id = ''
    compound.last_modified = 1612996604
    compound.pubchem_url = 'http://pubchem.ncbi.nlm.nih.gov/compound/13730'
    return compound


@pytest.fixture
def rt_reference():
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


@pytest.fixture
def mz_reference():
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


@pytest.fixture
def compound_identification(compound, rt_reference, mz_reference):
    ident = metob.CompoundIdentification()
    ident.unique_id = "18737c7141cc4efaa4545bead13ac751"
    ident.username = "wjholtz"
    ident.description = "No description"
    ident.creation_time = 1613002849
    ident.last_modified = 1613002979
    ident.identification_grade = None
    ident.compound = [compound]
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
    ident.mz_references = [mz_reference]
    ident.rt_references = [rt_reference]
    return ident


@pytest.fixture
def atlas(compound_identification):
    small_atlas = metob.Atlas()
    small_atlas.compound_identifications = [compound_identification]
    return small_atlas


@pytest.fixture
def atlas_two_compounds(compound_identification):
    small_atlas = metob.Atlas()
    small_atlas.compound_identifications = [compound_identification, compound_identification]
    return small_atlas


@pytest.fixture
def lcmsrun():
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


@pytest.fixture
def group(lcmsrun):
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
