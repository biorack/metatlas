# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long, duplicate-code

from . import utils


def test_targeted_by_line01_with_remove(tmp_path):
    image = "registry.spin.nersc.gov/metatlas_test/metatlas_ci:1.0.0"
    expected = {}
    expected[
        str(
            tmp_path
            / "505892_OakGall_final/Test-HILIC/0/0/Targeted/Test-HILIC/EMA-POS/POS_data_sheets/POS_peak_height.tab"
        )
    ] = """group	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S1	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S2	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S3	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S4
file	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5
short groupname	POS_Cone-S1	POS_Cone-S2	POS_Cone-S3	POS_Cone-S4
sample treatment	Cone-S1	Cone-S2	Cone-S3	Cone-S4
short filename	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1
short samplename	POS_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1	POS_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1	POS_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1	POS_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1
0000_2deoxyadenosine_positive_M+H252p1091_2p20	3.047619062e+05	4.167880312e+05	8.376620625e+05	2.359861250e+06
0001_adenine_positive_M+H136p0618_2p52	1.594753875e+06	1.209648500e+07	5.177495600e+07	9.195548800e+07
0002_adenosine_positive_M+H268p1041_3p02	2.661186800e+07	1.197741840e+08	2.677188800e+08	4.739050240e+08
0003_sucrose_positive_M+Na365p1054_13p41	1.215929500e+07	5.998378500e+06	5.243578125e+05	3.552174750e+06"""
    expected[
        str(
            tmp_path
            / "505892_OakGall_final/Test-HILIC/0/0/Targeted/Test-HILIC/EMA-POS/POS_data_sheets/POS_rt_peak.tab"
        )
    ] = """group	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S1	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S2	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S3	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_Cone-S4
file	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5	20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5
short groupname	POS_Cone-S1	POS_Cone-S2	POS_Cone-S3	POS_Cone-S4
sample treatment	Cone-S1	Cone-S2	Cone-S3	Cone-S4
short filename	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1	20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1
short samplename	POS_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1	POS_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1	POS_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1	POS_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1
0000_2deoxyadenosine_positive_M+H252p1091_2p20	2.277504444e+00	2.280636311e+00	2.283326864e+00	2.292241573e+00
0001_adenine_positive_M+H136p0618_2p52	2.616474867e+00	2.639369249e+00	2.618291378e+00	2.657374620e+00
0002_adenosine_positive_M+H268p1041_3p02	3.098848820e+00	3.125092983e+00	3.117606878e+00	3.139331818e+00
0003_sucrose_positive_M+Na365p1054_13p41	1.339970875e+01	1.339875221e+01	1.342395210e+01	1.340112782e+01"""
    command = """\
                    jq -M '(.cells[] | select(.source[] | contains("compound_idx=0")).source) \
                               += ["\\n", \
                                   "agui.compound_idx = 0\\n", \
                                   "agui.set_msms_flag(\\"1, co-isolated precursor but all reference ions are in sample spectrum\\")\\n", \
                                   "agui.data.set_rt(0, \\"rt_min\\", 2.1245)\\n", \
                                   "agui.data.set_rt(0, \\"rt_max\\", 2.4439)\\n", \
                                   "agui.compound_idx = 1\\n", \
                                   "agui.set_peak_flag(\\"remove\\")\\n", \
                                   "agui.compound_idx = 2\\n", \
                                   "agui.set_msms_flag(\\"1, perfect match to internal reference library\\")\\n", \
                                   "agui.data.set_rt(2, \\"rt_min\\", 2.4361)\\n", \
                                   "agui.data.set_rt(2, \\"rt_max\\", 2.8608)\\n", \
                                   "agui.compound_idx = 3\\n", \
                                   "agui.set_msms_flag(\\"1, perfect match to internal reference library\\")\\n", \
                                   "agui.data.set_rt(3, \\"rt_min\\", 2.8428)\\n", \
                                   "agui.data.set_rt(3, \\"rt_max\\", 3.3081)\\n", \
                                   "agui.compound_idx = 4\\n", \
                                   "agui.set_peak_flag(\\"remove\\")\\n", \
                                   "agui.compound_idx = 5\\n", \
                                   "agui.data.set_rt(5, \\"rt_min\\", 13.319)\\n", \
                                   "agui.data.set_rt(5, \\"rt_max\\", 13.520)\\n" \
                                  ]' /src/notebooks/reference/Targeted.ipynb > /out/Remove.ipynb &&  \
                    papermill -k papermill \
                        -p source_atlas_name HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0 \
                        -p source_atlas_unique_id 4b05837a53494dd8b680e6b5059e1934 \
                        -p experiment 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 \
                        -p project_directory /out \
                        -p max_cpus 2 \
                        -p config_file_name /src/test_config.yaml \
                        -p workflow_name Test-HILIC\
                        -p analysis_name EMA-POS \
                        -p rt_alignment_number 0 \
                        -p analysis_number 0 \
                        -y "exclude_lcmsruns: ['20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_53_Cone-S1_5_Rg70to1050-CE102040-QlobataAkingi-S1_Run187', '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_54_Cone-S1_6_Rg70to1050-CE102040-QlobataAkingi-S1_Run221']" \
                        /out/Remove.ipynb \
                        /out/Remove-done.ipynb
                   """
    utils.exec_docker(image, command, tmp_path, {})
    assert utils.num_files_in(tmp_path) == 63
    utils.assert_files_match(expected)
