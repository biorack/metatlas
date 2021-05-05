import os
import subprocess


def test_targeted_by_line(tmp_path):
    image = 'registry.spin.nersc.gov/metatlas_test/metatlas_ci01:v1.0.0'
    experiment = '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583'
    out_file = tmp_path/experiment/'root0/FinalEMA-HILIC/POS_data_sheets/POS_peak_height.tab'
    expected = [
        f"group\t{experiment}_POS_MSMS_root0_Cone-S1\t{experiment}_POS_MSMS_root0_Cone-S2\t{experiment}_POS_MSMS_root0_Cone-S3\t{experiment}_POS_MSMS_root0_Cone-S4",  # noqa: E501
        f"file\t{experiment}_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5\t{experiment}_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5\t{experiment}_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5\t{experiment}_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5",  # noqa: E501
        'short groupname\tPOS_Cone-S1\tPOS_Cone-S2\tPOS_Cone-S3\tPOS_Cone-S4',
        'sample treatment\tCone-S1\tCone-S2\tCone-S3\tCone-S4',
        'short filename\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1\t20201106_PS-KM_OakGall_final_HILICZ_POS_Rg70to1050-CE102040-QlobataAkingi-S1',  # noqa: E501
        'short samplename\tPOS_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1\tPOS_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1',  # noqa: E501
        '0000_2deoxyadenosine_positive_M+H252p1091_2p20\t304761.90625\t416788.03125\t837662.0625\t2359861.25',
        '0001_adenine_positive_M+H136p0618_2p52\t1645281.875\t12096485.0\t51774956.0\t91955488.0',
        '0002_adenine_positive_M+H136p0618_2p52\t1880780.125\t12096485.0\t51774956.0\t91955488.0',
        '0003_xanthine_positive_M+H153p0407_2p70\t72926.875\t60128.625\t231272.640625\t317968.03125',
        '0004_4-pyridoxic_acid_positive_M+H184p0605_2p84\t44113.671875\t66073.203125\t94702.390625\t214180.296875',  # noqa: E501
        '0005_adenosine_positive_M+H268p1041_3p02\t26611868.0\t119774184.0\t267718880.0\t473905024.0',
        ''
    ]
    subprocess.run(["docker", "run",
                    "--rm",
                    "-v", f"{os.getcwd()}:/src",
                    "-v", f"{tmp_path}:/out",
                    image,
                    "papermill",
                    "-p", "experiment", experiment,
                    "-p", "metatlas_repo_path", "/src",
                    "-p", "project_directory", "/out",
                    "-p", "max_cpus", "2",
                    "/src/notebooks/reference/Targeted.ipynb",
                    "/out/Targeted.ipynb"],
                   check=True)
    with open(out_file, 'r') as handle:
        for num, line in enumerate(handle.readlines()):
            clean_line = line.rstrip('\n')
            assert expected[num] == clean_line
