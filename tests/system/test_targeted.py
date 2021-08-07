# pylint: disable=missing-function-docstring, missing-module-docstring, line-too-long

import os
import subprocess


def test_targeted_by_line01_with_remove(tmp_path):
    image = "registry.spin.nersc.gov/metatlas_test/metatlas_ci01:v1.4.0"
    experiment = "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583"
    out_files = {}
    expected = {}
    out_files["peak_height"] = (
        tmp_path / experiment / "root0/FinalEMA-HILIC/POS_data_sheets/POS_peak_height.tab"
    )
    expected["peak_height"] = [
        f"group\t{experiment}_POS_MSMS_root0_Cone-S1\t{experiment}_POS_MSMS_root0_Cone-S2\t{experiment}_POS_MSMS_root0_Cone-S3\t{experiment}_POS_MSMS_root0_Cone-S4",  # noqa: E501
        f"file\t{experiment}_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\t{experiment}_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5\t{experiment}_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5\t{experiment}_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5",  # noqa: E501
        "short groupname\tPOS_Cone-S1\tPOS_Cone-S2\tPOS_Cone-S3\tPOS_Cone-S4",
        "sample treatment\tCone-S1\tCone-S2\tCone-S3\tCone-S4",
        "short filename\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1",  # noqa: E501
        "short samplename\tPOS_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1",  # noqa: E501
        "0000_2deoxyadenosine_positive_M+H252p1091_2p20\t304761.90625\t416788.03125\t837662.0625\t2359861.25",
        "0001_adenine_positive_M+H136p0618_2p52\t1594753.875\t12096485.0\t51774956.0\t91955488.0",
        "0002_adenosine_positive_M+H268p1041_3p02\t26611868.0\t119774184.0\t267718880.0\t473905024.0",
        "",
    ]
    out_files["rt_peak"] = tmp_path / experiment / "root0/FinalEMA-HILIC/POS_data_sheets/POS_rt_peak.tab"
    expected["rt_peak"] = [
        f"group\t{experiment}_POS_MSMS_root0_Cone-S1\t{experiment}_POS_MSMS_root0_Cone-S2\t{experiment}_POS_MSMS_root0_Cone-S3\t{experiment}_POS_MSMS_root0_Cone-S4",
        f"file\t{experiment}_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\t{experiment}_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5\t{experiment}_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5\t{experiment}_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5",
        "short groupname\tPOS_Cone-S1\tPOS_Cone-S2\tPOS_Cone-S3\tPOS_Cone-S4",
        "sample treatment\tCone-S1\tCone-S2\tCone-S3\tCone-S4",
        "short filename\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1",  # noqa: E501
        "short samplename\tPOS_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1",  # noqa: E501
        "0000_2deoxyadenosine_positive_M+H252p1091_2p20\t2.2775044441223145\t2.2806363105773926\t2.2833268642425537\t2.2922415733337402",
        "0001_adenine_positive_M+H136p0618_2p52\t2.6164748668670654\t2.639369249343872\t2.6182913780212402\t2.657374620437622",
        "0002_adenosine_positive_M+H268p1041_3p02\t3.098848819732666\t3.1250929832458496\t3.1176068782806396\t3.139331817626953",
    ]
    subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{os.getcwd()}:/src",
            "-v",
            f"{tmp_path}:/out",
            image,
            "/bin/bash",
            "-c",
            """\
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
                                   "agui.data.set_rt(3, \\"rt_max\\", 3.3081)\\n" \
                                  ]' /src/notebooks/reference/Targeted.ipynb > /out/Remove.ipynb &&  \
                    papermill \
                        -p source_atlas HILICz150_ANT20190824_PRD_EMA_Unlab_POS_20201106_505892_root0 \
                        -p experiment 20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583 \
                        -p metatlas_repo_path /src \
                        -p project_directory /out \
                        -p max_cpus 2 \
                        /out/Remove.ipynb \
                        /out/Remove-done.ipynb
                   """,
        ],
        check=True,
    )
    files = subprocess.check_output(f"find {str(tmp_path)} -type f", shell=True, text=True).strip()
    print(files)
    num_files_created = int(
        subprocess.check_output(f"find {str(tmp_path)} -type f | wc -l", shell=True, text=True).strip()
    )
    for _, path in out_files.items():
        os.system(f"cat {path}")
    assert num_files_created == 39
    for metric_name, path in out_files.items():
        with open(path, "r") as handle:
            for num, line in enumerate(handle.readlines()):
                clean_line = line.rstrip("\n")
                assert expected[metric_name][num] == clean_line
