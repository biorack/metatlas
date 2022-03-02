# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long, duplicate-code

from . import utils


def test_c18_by_line01_with_remove(tmp_path):
    image = "registry.spin.nersc.gov/metatlas_test/metatlas_ci03:v0.0.5"
    experiment = "20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680"
    expected = {}
    expected[
        str(tmp_path / experiment / "root0/FinalEMA-C18/NEG/NEG_data_sheets/NEG_peak_height.tab")
    ] = """group	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_ExCtrl	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_Neg-D30	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_Neg-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_Neg-D89	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S16-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S16-D89	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S32-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S40-D30	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S40-D89	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S53-D30	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S53-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S53-D89
file	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_75_ExCtrl_C_Rg80to1200-CE102040-soil-S1_Run209.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_57_Neg-D30_C_Rg80to1200-CE102040-soil-S1_Run224.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_60_Neg-D45_C_Rg80to1200-CE102040-soil-S1_Run230.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_63_Neg-D89_C_Rg80to1200-CE102040-soil-S1_Run215.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_12_S16-D45_C_Rg80to1200-CE102040-soil-S1_Run203.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_21_S16-D89_C_Rg80to1200-CE102040-soil-S1_Run221.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_18_S32-D45_C_Rg80to1200-CE102040-soil-S1_Run236.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_33_S40-D30_C_Rg80to1200-CE102040-soil-S1_Run233.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_51_S40-D89_C_Rg80to1200-CE102040-soil-S1_Run227.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_36_S53-D30_C_Rg80to1200-CE102040-soil-S1_Run218.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_45_S53-D45_C_Rg80to1200-CE102040-soil-S1_Run212.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_54_S53-D89_C_Rg80to1200-CE102040-soil-S1_Run206.h5
short groupname	NEG_ExCtrl	NEG_Neg-D30	NEG_Neg-D45	NEG_Neg-D89	NEG_S16-D45	NEG_S16-D89	NEG_S32-D45	NEG_S40-D30	NEG_S40-D89	NEG_S53-D30	NEG_S53-D45	NEG_S53-D89
sample treatment	ExCtrl	Neg-D30	Neg-D45	Neg-D89	S16-D45	S16-D89	S32-D45	S40-D30	S40-D89	S53-D30	S53-D45	S53-D89
short filename	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1
short samplename	NEG_ExCtrl_C_Rg80to1200-CE102040-soil-S1	NEG_Neg-D30_C_Rg80to1200-CE102040-soil-S1	NEG_Neg-D45_C_Rg80to1200-CE102040-soil-S1	NEG_Neg-D89_C_Rg80to1200-CE102040-soil-S1	NEG_S16-D45_C_Rg80to1200-CE102040-soil-S1	NEG_S16-D89_C_Rg80to1200-CE102040-soil-S1	NEG_S32-D45_C_Rg80to1200-CE102040-soil-S1	NEG_S40-D30_C_Rg80to1200-CE102040-soil-S1	NEG_S40-D89_C_Rg80to1200-CE102040-soil-S1	NEG_S53-D30_C_Rg80to1200-CE102040-soil-S1	NEG_S53-D45_C_Rg80to1200-CE102040-soil-S1	NEG_S53-D89_C_Rg80to1200-CE102040-soil-S1
0000_azelaic_acid_negative_M-H187p0976_3p49	1.486730560e+08	2.335685440e+08	1.853421120e+08	3.217044000e+07	9.402386000e+06	5.876234000e+06	3.349810200e+07	1.178001700e+07	6.884600500e+06	1.464379700e+07	7.114720000e+06	4.965849000e+06
0001_vulpinic_acid_negative_M-H321p0768_6p26						2.963439258e+04	9.831337500e+05	1.175015332e+04	1.064087402e+04	1.329986625e+06	4.145107812e+05	1.812218750e+05"""

    expected[
        str(tmp_path / experiment / "root0/FinalEMA-C18/NEG/NEG_data_sheets/NEG_rt_peak.tab")
    ] = """group	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_ExCtrl	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_Neg-D30	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_Neg-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_Neg-D89	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S16-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S16-D89	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S32-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S40-D30	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S40-D89	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S53-D30	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S53-D45	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_root0_S53-D89
file	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_75_ExCtrl_C_Rg80to1200-CE102040-soil-S1_Run209.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_57_Neg-D30_C_Rg80to1200-CE102040-soil-S1_Run224.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_60_Neg-D45_C_Rg80to1200-CE102040-soil-S1_Run230.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_63_Neg-D89_C_Rg80to1200-CE102040-soil-S1_Run215.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_12_S16-D45_C_Rg80to1200-CE102040-soil-S1_Run203.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_21_S16-D89_C_Rg80to1200-CE102040-soil-S1_Run221.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_18_S32-D45_C_Rg80to1200-CE102040-soil-S1_Run236.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_33_S40-D30_C_Rg80to1200-CE102040-soil-S1_Run233.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_51_S40-D89_C_Rg80to1200-CE102040-soil-S1_Run227.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_36_S53-D30_C_Rg80to1200-CE102040-soil-S1_Run218.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_45_S53-D45_C_Rg80to1200-CE102040-soil-S1_Run212.h5	20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680_NEG_MSMS_54_S53-D89_C_Rg80to1200-CE102040-soil-S1_Run206.h5
short groupname	NEG_ExCtrl	NEG_Neg-D30	NEG_Neg-D45	NEG_Neg-D89	NEG_S16-D45	NEG_S16-D89	NEG_S32-D45	NEG_S40-D30	NEG_S40-D89	NEG_S53-D30	NEG_S53-D45	NEG_S53-D89
sample treatment	ExCtrl	Neg-D30	Neg-D45	Neg-D89	S16-D45	S16-D89	S32-D45	S40-D30	S40-D89	S53-D30	S53-D45	S53-D89
short filename	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1	20210915_MK_SoilWaterRep_final_C18_NEG_Rg80to1200-CE102040-soil-S1
short samplename	NEG_ExCtrl_C_Rg80to1200-CE102040-soil-S1	NEG_Neg-D30_C_Rg80to1200-CE102040-soil-S1	NEG_Neg-D45_C_Rg80to1200-CE102040-soil-S1	NEG_Neg-D89_C_Rg80to1200-CE102040-soil-S1	NEG_S16-D45_C_Rg80to1200-CE102040-soil-S1	NEG_S16-D89_C_Rg80to1200-CE102040-soil-S1	NEG_S32-D45_C_Rg80to1200-CE102040-soil-S1	NEG_S40-D30_C_Rg80to1200-CE102040-soil-S1	NEG_S40-D89_C_Rg80to1200-CE102040-soil-S1	NEG_S53-D30_C_Rg80to1200-CE102040-soil-S1	NEG_S53-D45_C_Rg80to1200-CE102040-soil-S1	NEG_S53-D89_C_Rg80to1200-CE102040-soil-S1
0000_azelaic_acid_negative_M-H187p0976_3p49	3.446131468e+00	3.445808887e+00	3.444513321e+00	3.452524185e+00	3.443904161e+00	3.448769808e+00	3.440114975e+00	3.452755451e+00	3.451116800e+00	3.445648670e+00	3.447070122e+00	3.442815304e+00
0001_vulpinic_acid_negative_M-H321p0768_6p26						6.212389946e+00	6.196495533e+00	6.216785908e+00	6.211359978e+00	6.209033012e+00	6.213802814e+00	6.207533360e+00"""

    command = f"""\
                    jq -M '(.cells[] | select(.source[] | contains("compound_idx=0")).source) \
                               += ["\\n", \
                                   "agui.compound_idx = 0\\n", \
                                   "agui.set_msms_flag(\\"1, perfect match to internal reference library\\")\\n", \
                                   "agui.data.set_rt(0, \\"rt_min\\", 3.404)\\n", \
                                   "agui.data.set_rt(0, \\"rt_max\\", 3.528)\\n", \
                                   "agui.compound_idx = 1\\n", \
                                   "agui.set_msms_flag(\\"1, perfect match to internal reference library\\")\\n", \
                                   "agui.data.set_rt(1, \\"rt_min\\", 6.1681)\\n", \
                                   "agui.data.set_rt(1, \\"rt_max\\", 6.2470)\\n" \
                                  ]' /src/notebooks/reference/Targeted.ipynb > /out/Remove.ipynb &&  \
                    papermill -k papermill \
                        -p source_atlas C18_20220215_TPL_EMA_Unlab_NEG \
                        -p experiment {experiment} \
                        -p polarity negative \
                        -p output_type FinalEMA-C18 \
                        -p project_directory /out \
                        -p num_points None \
                        -p peak_height None \
                        -p msms_score  0.5 \
                        -p max_cpus 2 \
                        /out/Remove.ipynb \
                        /out/Remove-done.ipynb
                   """
    utils.exec_docker(image, command, tmp_path)
    assert utils.num_files_in(tmp_path) == 50
    utils.assert_files_match(expected)
