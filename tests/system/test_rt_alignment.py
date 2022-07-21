# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long, duplicate-code

import numpy as np

from . import utils


def test_rt_alignment_by_line01(tmp_path):
    image = "registry.spin.nersc.gov/metatlas_test/metatlas_ci02:v1.4.24"
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    expected = {}
    expected[
        str(tmp_path / experiment / "root_0_0/Targeted/Test-QC/RT_Alignment/rt_alignment_model.txt")
    ] = """RANSACRegressor(random_state=42)
Linear model with intercept=0.430 and slope=0.95574
groups = 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_root_0_0_QC, 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_QC
atlas = HILICz150_ANT20190824_TPL_QCv3_Unlab_POS

LinearRegression()
Polynomial model with intercept=0.733 and coefficents=[0.00000, 0.81321, 0.00919]
groups = 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_root_0_0_QC, 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_root_0_0_QC
atlas = HILICz150_ANT20190824_TPL_QCv3_Unlab_POS
"""
    expected_df = {}
    expected_df[
        str(
            tmp_path
            / experiment
            / "root_0_0/Targeted/Test-QC/RT_Alignment/RT_Alignment_Model_Comparison.csv"
        )
    ] = {
        "Unnamed: 0": {
            0: "0000_uracil_unlabeled_positive_M+H113p0346_1p39",
            1: "0001_cytosine_unlabeled_positive_M+H112p0505_4p83",
            2: "0002_sucrose_unlabeled_positive_M+H343p1235_13p45",
        },
        "RT Measured": {0: 1.884217, 1: 4.878586, 2: 13.32883},
        "RT Reference": {0: 1.3937, 1: 4.833664, 2: 13.44515},
        "Relative RT Linear": {0: 1.761967, 1: 5.049671, 2: 13.28},
        "Relative RT Polynomial": {0: 1.884217, 1: 4.878586, 2: 13.32883},
        "RT Diff Linear": {0: 0.1222508, 1: -0.1710854, 2: 0.04883459},
        "RT Diff Polynomial": {0: -1.865175e-14, 1: 1.776357e-15, 2: 1.776357e-14},
    }

    expected_df[str(tmp_path / experiment / "root_0_0/Targeted/Test-QC/QC-POS/POS_QC_Measured_RTs.csv")] = {
        "Unnamed: 0": {
            0: "0000_uracil_unlabeled_positive_M+H113p0346_1p39",
            1: "0001_cytosine_unlabeled_positive_M+H112p0505_4p83",
            2: "0002_sucrose_unlabeled_positive_M+H343p1235_13p45",
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Pre_Rg70to1050-CE102040--QC_Run6.h5": {
            0: 1.879153,
            1: 4.836015,
            2: np.NaN,
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Pre_Rg70to1050-CE102040--QC_Run7.h5": {
            0: 1.889282,
            1: 4.846456,
            2: np.NaN,
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Post_Rg70to1050-CE102040--QC_Run307.h5": {
            0: 1.878133,
            1: 4.910717,
            2: 13.34303,
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Post_Rg70to1050-CE102040--QC_Run308.h5": {
            0: 1.892506,
            1: 4.921309,
            2: 13.31463,
        },
        "mean": {0: 1.884768, 1: 4.878624, 2: 13.32883},
        "median": {0: 1.884217, 1: 4.878586, 2: 13.32883},
        "min": {0: 1.878133, 1: 4.836015, 2: 13.31463},
        "max": {0: 1.892506, 1: 4.921309, 2: 13.34303},
        "standard deviation": {0: 0.007206614, 1: 0.04359777, 2: 0.02008078},
        "standard error": {0: 0.003603307, 1: 0.02179888, 2: 0.01419926},
        "#NaNs": {0: 0, 1: 0, 2: 0},
    }

    command = """papermill -k papermill \
                        -p experiment 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 \
                        -p project_directory /out \
                        -p max_cpus 2 \
                        -p config_file_name /src/metatlas_config.yaml \
                        -p workflow_name Test-QC \
                        -p rt_alignment_number 0 \
                        /src/notebooks/reference/RT_Alignment.ipynb \
                        /out/Remove-done.ipynb
    """
    utils.exec_docker(image, command, tmp_path, utils.PAPERMILL_ENV)
    assert utils.num_files_in(tmp_path) == 21
    utils.assert_files_match(expected)
    utils.assert_dfs_match(expected_df)
