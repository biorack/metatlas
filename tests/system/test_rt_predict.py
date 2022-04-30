# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long, duplicate-code

from . import utils


def test_rt_predict_by_line01(tmp_path):
    image = "registry.spin.nersc.gov/metatlas_test/metatlas_ci02:v1.4.21"
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    expected = {}
    expected[
        str(tmp_path / experiment / "root0/data_QC/rt_model.txt")
    ] = """RANSACRegressor(random_state=42)
Linear model with intercept=0.430 and slope=0.95574
groups = 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_root0_QC, 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root0_QC
atlas = HILICz150_ANT20190824_TPL_QCv3_Unlab_POS

LinearRegression()
Polynomial model with intercept=0.733 and coefficents=[0.00000, 0.81321, 0.00919]
groups = 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_root0_QC, 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root0_QC
atlas = HILICz150_ANT20190824_TPL_QCv3_Unlab_POS
"""
    expected[
        str(tmp_path / experiment / "root0/data_QC/RT_Predicted_Model_Comparison.csv")
    ] = """,RT Measured,RT Reference,RT Linear Pred,RT Polynomial Pred,RT Diff Linear,RT Diff Polynomial
0000_uracil_unlabeled_positive_M+H113p0346_1p39,1.884217e+00,1.393700e+00,1.761967e+00,1.884217e+00,1.222508e-01,-1.865175e-14
0001_cytosine_unlabeled_positive_M+H112p0505_4p83,4.878586e+00,4.833664e+00,5.049671e+00,4.878586e+00,-1.710854e-01,1.776357e-15
0002_sucrose_unlabeled_positive_M+H343p1235_13p45,1.332883e+01,1.344515e+01,1.328000e+01,1.332883e+01,4.883459e-02,1.776357e-14"""

    expected[
        str(tmp_path / experiment / "root0/data_QC/QC_Measured_RTs.csv")
    ] = """,20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Pre_Rg70to1050-CE102040--QC_Run6.h5,20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Pre_Rg70to1050-CE102040--QC_Run7.h5,20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Post_Rg70to1050-CE102040--QC_Run307.h5,20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Post_Rg70to1050-CE102040--QC_Run308.h5,mean,median,min,max,standard deviation,standard error,#NaNs
0000_uracil_unlabeled_positive_M+H113p0346_1p39,1.879153e+00,1.889282e+00,1.878133e+00,1.892506e+00,1.884768e+00,1.884217e+00,1.878133e+00,1.892506e+00,7.206614e-03,3.603307e-03,0
0001_cytosine_unlabeled_positive_M+H112p0505_4p83,4.836015e+00,4.846456e+00,4.910717e+00,4.921309e+00,4.878624e+00,4.878586e+00,4.836015e+00,4.921309e+00,4.359777e-02,2.179888e-02,0
0002_sucrose_unlabeled_positive_M+H343p1235_13p45,,,1.334303e+01,1.331463e+01,1.332883e+01,1.332883e+01,1.331463e+01,1.334303e+01,2.008078e-02,1.419926e-02,0"""

    command = """papermill -k papermill \
                        -p experiment 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 \
                        -p model_only False \
                        -p project_directory /out \
                        -p max_cpus 2 \
                        /src/notebooks/reference/RT_Prediction.ipynb \
                        /out/Remove-done.ipynb
    """
    utils.exec_docker(image, command, tmp_path)
    assert utils.num_files_in(tmp_path) == 63  # this is 11 if model_only is set to True
    utils.assert_files_match(expected)
