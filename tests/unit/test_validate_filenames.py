# pylint: disable=line-too-long,missing-function-docstring,missing-module-docstring
# pylint: disable=use-implicit-booleaness-not-comparison

from pathlib import Path

import metatlas.tools.validate_filenames as valid

# most of the functions under test here return a list of error messages
# so returning [] is the non error condition and returning a non-zero length
# array is the error condition


def test_extract_batch01():
    path = Path("/one_two_three_four_five_six_seven_eight_nine_ten/foo.raw")
    assert valid.extract_batch(path) == "one_two_three_four_five_six_seven_eight_nine"


def test_parent_dir_is_prefix01():
    path = Path("/foobar/foobar_zoop.raw")
    assert valid.parent_dir_is_prefix(path) == []


def test_parent_dir_is_prefix02():
    path = Path("/nope/foobar_zoop.raw")
    assert valid.parent_dir_is_prefix(path) != []


def test_has_minimum_num_fields01():
    path = Path("/foo/1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.has_minimum_num_fields(path) == []


def test_has_minimum_num_fields02():
    path = Path("/foo/1_2_3_4_5_6_7_8_9_10.raw")
    assert valid.has_minimum_num_fields(path) != []


def test_has_exactly_num_fields01():
    path = Path("/foo/1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.has_exactly_num_fields(path) == []


def test_has_exactly_num_fields02():
    path = Path("/foo/1_2_3_4_5_6_7_8_9_10_11_12_13.raw")
    assert valid.has_exactly_num_fields(path) != []


def test_valid_string01():
    path = Path("/foo/a_bbb_c123_333_(4)_5_6_7_8_9_10_11_12_13_♥.raw")
    assert valid.valid_string(path, field_num=0, min_len=1, max_len=1, alpha_only=True) == []
    assert valid.valid_string(path, field_num=0, min_len=2, max_len=2, alpha_only=True) != []
    assert valid.valid_string(path, field_num=1, min_len=4, max_len=9, alpha_only=True) != []
    assert valid.valid_string(path, field_num=2, min_len=4, max_len=9, alpha_only=True) != []
    assert valid.valid_string(path, field_num=2, min_len=4, max_len=9, alphanumeric_only=True) == []
    assert valid.valid_string(path, field_num=4, min_len=1, max_len=9, alphanumeric_only=True) != []
    assert valid.valid_string(path, field_num=4, min_len=1, max_len=9, alphanumeric_only=False) == []
    assert valid.valid_string(path, field_num=4, min_len=1, max_len=9, sub_field_num=9) != []
    assert valid.valid_string(path, field_num=14) != []


def test_get_alpha_prefix():
    assert valid.get_alpha_prefix("foobar0") == "foobar"
    assert valid.get_alpha_prefix("100") == ""
    assert valid.get_alpha_prefix("zoop") == "zoop"


def test_valid_int():
    path = Path("/foo/a_bbb_c12_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_int(path, field_num=0, min_int=0, max_int=99, allowed_prefixes=[]) != []
    assert valid.valid_int(path, field_num=2, min_int=0, max_int=99, allowed_prefixes=["c"]) == []
    assert valid.valid_int(path, field_num=2, min_int=0, max_int=99, allowed_prefixes=["foo"]) != []
    assert valid.valid_int(path, field_num=2, min_int=0, max_int=99, allowed_prefixes=[]) != []
    assert (
        valid.valid_int(
            path,
            field_num=2,
        )
        != []
    )
    assert valid.valid_int(path, field_num=3, min_int=0, max_int=999, allowed_prefixes=[]) == []
    assert valid.valid_int(path, field_num=3, min_int=0, max_int=2, allowed_prefixes=[]) != []
    assert valid.valid_int(path, field_num=3, min_int=9999) != []
    assert valid.valid_int(path, field_num=99) != []


def test_valid_field0():
    path = Path("/foo/20220101_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) == []
    path = Path("/foo/.raw")
    assert valid.valid_field0(path) != []
    path = Path("/foo/20220101X_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) != []
    path = Path("/foo/0220101_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) != []
    path = Path("/foo/19840101_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) != []
    path = Path("/foo/43210101_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) != []
    path = Path("/foo/20209901_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) != []
    path = Path("/foo/20200199_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field0(path) != []


def test_valid_field1():
    path = Path("/foo/20220101_WH_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) == []
    path = Path("/foo/20220101_Will_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) != []
    path = Path("/foo/20220101_W_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field1(path) != []


def test_valid_field2():
    path = Path("/foo/20220101_WH_WH_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) == []
    path = Path("/foo/20220101_WH_William_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) != []
    path = Path("/foo/20220101_WH_W_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) != []
    path = Path("/foo/20220101_WH_a-b-c_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field2(path) != []


def test_valid_field3():
    path = Path("/foo/20220101_WH_WH_123456_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field3(path) == []
    path = Path("/foo/20220101_WH_WH_project99_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field3(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field3(path) != []


def test_valid_field4():
    path = Path("/foo/20220101_WH_WH_123456_exp_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field4(path) == []
    path = Path("/foo/20220101_WH_WH_project99_(exp)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field4(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field4(path) != []


def test_valid_field5():
    path = Path("/foo/20220101_WH_WH_123456_5_sampset_6_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field5(path) == []
    path = Path("/foo/20220101_WH_WH_123456_5_thisisaverylongsampleset_6_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field5(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field5(path) != []


def test_valid_field6():
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_IDX_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field6(path) == []
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_idx_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field6(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field6(path) != []


def test_valid_field7():
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field7(path) == []
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_idx_C99-bar_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field7(path) != []
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_idx_012345678901234567_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field7(path) != []
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_IDX_a-b-c_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field7(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field7(path) != []


def test_valid_field8():
    path = Path("/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field8(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col*_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field8(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field8(path) != []


def test_parent_dir_is_prefix():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_is_prefix(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_extra/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_is_prefix(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_DIFFERENT_IDX_HILICZ-foo_a-x99-col_extra/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_is_prefix(path) != []


def test_parent_dir_num_fields():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_num_fields(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_extra/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_num_fields(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_extra_too_many/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_num_fields(path) != []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_9_10_11_12_13_14_15_16.raw"
    )
    assert valid.parent_dir_num_fields(path) != []


def test_valid_field9():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_10_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field9(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_NEG_10_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field9(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_POS_10_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field9(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_pos_10_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field9(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field9(path) != []


def test_valid_field10():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS2-a745_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS10-t123_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MSn_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS3-t1-a2_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) != []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_ms1_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) != []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS2-01234567890123456789_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) != []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_x_11_12_13_14_15_16.raw"
    )
    assert valid.valid_field10(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field10(path) != []


def test_valid_field11():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_12_13_14_15_16.raw"
    )
    assert valid.valid_field11(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01234567890123456789234567890123456789_12_13_14_15_16.raw"
    )
    assert valid.valid_field11(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field11(path) != []


def test_valid_field12():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_13_14_15_16.raw"
    )
    assert valid.valid_field12(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp#01_13_14_15_16.raw"
    )
    assert valid.valid_field12(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field12(path) != []


def test_valid_field13():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_14_15_16.raw"
    )
    assert valid.valid_field13(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_thisiswaytoolong_14_15_16.raw"
    )
    assert valid.valid_field13(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field13(path) != []


def test_valid_field14():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_15_16.raw"
    )
    assert valid.valid_field14(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_emptysubfields--_15_16.raw"
    )
    assert valid.valid_field14(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field14(path) != []


def test_valid_field15():
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_Run0.raw"
    )
    assert valid.valid_field15(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_000.raw"
    )
    assert valid.valid_field15(path) == []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_run0.raw"
    )
    assert valid.valid_field15(path) != []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_R0001.raw"
    )
    assert valid.valid_field15(path) != []
    path = Path(
        "/foo/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_12345678900123456.raw"
    )
    assert valid.valid_field15(path) != []
    path = Path("/foo/.raw")
    assert valid.valid_field15(path) != []


def test_field_not_found():
    path = Path("/foo/.raw")
    assert valid.valid_field0(path) != []
    assert valid.valid_field1(path) != []
    assert valid.valid_field2(path) != []
    assert valid.valid_field3(path) != []
    assert valid.valid_field4(path) != []
    assert valid.valid_field5(path) != []
    assert valid.valid_field6(path) != []
    assert valid.valid_field7(path) != []
    assert valid.valid_field8(path) != []
    assert valid.valid_field9(path) != []
    assert valid.valid_field10(path) != []
    assert valid.valid_field11(path) != []
    assert valid.valid_field12(path) != []
    assert valid.valid_field13(path) != []
    assert valid.valid_field14(path) != []
    assert valid.valid_field15(path) != []


def test_validate_file_name():
    path = Path(
        "/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_Run0.raw"
    )
    assert valid.validate_file_name(path, minimal=False)


def test_get_validation_messages():
    path = Path(
        "/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_Run0.raw"
    )
    assert valid.get_validation_messages(path, minimal=False) == ([], [])


def test_is_valid_file_name():
    path = Path(
        "/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_Run0.raw"
    )
    assert valid.is_valid_file_name(path, [valid.valid_field0], [lambda x: ["warn message"]])
    assert valid.is_valid_file_name(path, [valid.valid_field0], [])
    path = Path(
        "/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_foo_bar/20220101_WH_WH_123456_exp_samptset_IDX_HILICZ-foo_a-x99-col_FPS_MS1_Samp-01-a_grp-01_BR1-TR2_optional-dash-delim_Run0.raw"
    )
    assert not valid.is_valid_file_name(path, [valid.parent_dir_num_fields], [])
