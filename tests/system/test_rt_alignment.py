# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long, duplicate-code

import numpy as np

from . import utils


def test_rt_alignment_by_line01(tmp_path):
    image = "registry.spin.nersc.gov/metatlas_test/metatlas_ci:1.1.0"
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    expected = {}
    expected[
        str(
            tmp_path
            / f"{experiment}/root_Test-QC_0_0/Targeted/Test-QC_{experiment}/RT-Alignment/rt_alignment_model.txt"
        )
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
            / f"{experiment}/root_Test-QC_0_0/Targeted/Test-QC_{experiment}/RT-Alignment/RT-Alignment_Model_Comparison.csv"
        )
    ] = {
        "Unnamed: 0": {
            0: "0000_uracil_unlabeled_positive_M+H113p0346_1p39",
            1: "0001_2deoxyadenosine_unlabeled_positive_M+H252p1091_2p23",
            2: "0002_adenine_unlabeled_positive_M+H136p0618_2p56",
            3: "0003_xanthine_unlabeled_positive_M+H153p0407_2p73",
            4: "0004_adenosine_unlabeled_positive_M+H268p1040_3p09",
            5: "0005_cytosine_unlabeled_positive_M+H112p0505_4p83",
            6: "0006_sucrose_unlabeled_positive_M+H343p1235_13p45",
        },
        "RT Measured": {
            0: 1.884217,
            1: 2.293711,
            2: 2.664679,
            3: 2.770917,
            4: 3.109064,
            5: 4.878586,
            6: 13.32883,
        },
        "RT Reference": {
            0: 1.3937,
            1: 2.234007,
            2: 2.557602,
            3: 2.725344,
            4: 3.091019,
            5: 4.833664,
            6: 13.44515,
        },
        "Relative RT Linear": {
            0: 1.761967,
            1: 2.56508,
            2: 2.874353,
            3: 3.03467,
            4: 3.384159,
            5: 5.049671,
            6: 13.28,
        },
        "Relative RT Polynomial": {
            0: 1.884217,
            1: 2.595589,
            2: 2.872996,
            3: 3.017553,
            4: 3.334478,
            5: 4.878586,
            6: 13.32883,
        },
        "RT Diff Linear": {
            0: 0.1222508,
            1: -0.2713695,
            2: -0.2096735,
            3: -0.263753,
            4: -0.2750952,
            5: -0.1710854,
            6: 0.04883459,
        },
        "RT Diff Polynomial": {
            0: -1.865175e-14,
            1: -0.3018785,
            2: -0.2083166,
            3: -0.2466357,
            4: -0.2254138,
            5: 1.776357e-15,
            6: 1.776357e-14,
        },
    }

    expected_df[
        str(
            tmp_path
            / f"{experiment}/root_Test-QC_0_0/Targeted/Test-QC_{experiment}/QC-POS/POS_QC_Measured_RTs.csv"
        )
    ] = {
        "Unnamed: 0": {
            0: "0000_uracil_unlabeled_positive_M+H113p0346_1p88",
            1: "0001_2deoxyadenosine_unlabeled_positive_M+H252p1091_2p60",
            2: "0002_adenine_unlabeled_positive_M+H136p0618_2p87",
            3: "0003_xanthine_unlabeled_positive_M+H153p0407_3p02",
            4: "0004_adenosine_unlabeled_positive_M+H268p1040_3p33",
            5: "0005_cytosine_unlabeled_positive_M+H112p0505_4p88",
            6: "0006_sucrose_unlabeled_positive_M+H343p1235_13p33",
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Pre_Rg70to1050-CE102040--QC_Run6.h5": {
            0: 1.940126,
            1: 2.24759,
            2: 2.59303,
            3: 2.747974,
            4: 3.075452,
            5: 4.836015,
            6: np.NaN,
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Pre_Rg70to1050-CE102040--QC_Run7.h5": {
            0: 1.938525,
            1: 2.275779,
            2: 2.625245,
            3: 2.775257,
            4: 3.078128,
            5: 4.846456,
            6: np.NaN,
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_FPS_MS1_0_QC_Post_Rg70to1050-CE102040--QC_Run307.h5": {
            0: 1.926795,
            1: 2.322836,
            2: 2.704113,
            3: 2.766577,
            4: 3.146962,
            5: 4.910717,
            6: 13.34303,
        },
        "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_0_QC_Post_Rg70to1050-CE102040--QC_Run308.h5": {
            0: 1.932271,
            1: 2.311643,
            2: 2.719899,
            3: 2.788285,
            4: 3.14,
            5: 4.921309,
            6: 13.31463,
        },
        "mean": {0: 1.934429, 1: 2.289462, 2: 2.660572, 3: 2.769523, 4: 3.110135, 5: 4.878624, 6: 13.32883},
        "median": {0: 1.935398, 1: 2.293711, 2: 2.664679, 3: 2.770917, 4: 3.109064, 5: 4.878586, 6: 13.32883},
        "min": {0: 1.926795, 1: 2.24759, 2: 2.59303, 3: 2.747974, 4: 3.075452, 5: 4.836015, 6: 13.31463},
        "max": {0: 1.940126, 1: 2.322836, 2: 2.719899, 3: 2.788285, 4: 3.146962, 5: 4.921309, 6: 13.34303},
        "standard deviation": {
            0: 0.006114746,
            1: 0.03438155,
            2: 0.06117022,
            3: 0.01691088,
            4: 0.0386243,
            5: 0.04359777,
            6: 0.02008078,
        },
        "standard error": {
            0: 0.003057373,
            1: 0.01719078,
            2: 0.03058511,
            3: 0.008455441,
            4: 0.01931215,
            5: 0.02179888,
            6: 0.01419926,
        },
        "#NaNs": {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0},
    }

    command = f"""papermill -k papermill \
                        -p experiment {experiment} \
                        -p project_directory /out \
                        -p max_cpus 2 \
                        -p config_file_name /src/test_config.yaml \
                        -p workflow_name Test-QC \
                        -p rt_alignment_number 0 \
                        -y "inchi_keys_not_in_model: ['GFFGJBXGBJISGV-UHFFFAOYSA-N', 'LRFVTYWOQMYALW-UHFFFAOYSA-N', 'OIRDTQYFTABQOQ-KQYNXXCUSA-N', 'OLXZPDWKRNYJJZ-RRKCRQDMSA-N']" \
                        /src/notebooks/reference/RT-Alignment.ipynb \
                        /out/Remove-done.ipynb
    """
    utils.exec_docker(image, command, tmp_path, utils.PAPERMILL_ENV)
    assert utils.num_files_in(tmp_path) == 26
    utils.assert_files_match(expected)
    utils.assert_dfs_match(expected_df)
