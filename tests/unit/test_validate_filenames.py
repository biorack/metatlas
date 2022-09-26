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


def test_num_subfields01():
    path = Path("/foo/a_1-2-3_b_c.raw")
    assert valid.num_subfields(path, 1, 2, 3) == []


def test_num_subfields02():
    path = Path("/foo/a_1-2-3_b_c.raw")
    assert valid.num_subfields(path, 1, 2, 2) != []


def test_valid_string01():
    path = Path("/foo/a_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_string(path, field_num=0, min_len=1, max_len=1, alpha_only=True) == []
    assert valid.valid_string(path, field_num=0, min_len=2, max_len=2, alpha_only=True) != []
    assert valid.valid_string(path, field_num=1, min_len=4, max_len=9, alpha_only=True) != []
    assert valid.valid_string(path, field_num=2, min_len=4, max_len=9, alpha_only=True) != []
    assert valid.valid_string(path, field_num=2, min_len=4, max_len=9, alphanumeric_only=True) == []
    assert valid.valid_string(path, field_num=4, min_len=1, max_len=9, alphanumeric_only=True) != []
    assert valid.valid_string(path, field_num=4, min_len=1, max_len=9, alphanumeric_only=False) == []


def test_get_alpha_prefix():
    assert valid.get_alpha_prefix("foobar0") == "foobar"
    assert valid.get_alpha_prefix("100") == ""
    assert valid.get_alpha_prefix("zoop") == "zoop"


def test_valid_int():
    path = Path("/foo/a_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_int(path, field_num=0, min_int=0, max_int=99, allowed_prefixes=[]) != []
    assert valid.valid_int(path, field_num=2, min_int=0, max_int=99, allowed_prefixes=["c"]) == []
    assert valid.valid_int(path, field_num=2, min_int=0, max_int=99, allowed_prefixes=["foo"]) != []
    assert valid.valid_int(path, field_num=2, min_int=0, max_int=99, allowed_prefixes=[]) != []
    assert valid.valid_int(path, field_num=3, min_int=0, max_int=999, allowed_prefixes=[]) == []
    assert valid.valid_int(path, field_num=3, min_int=0, max_int=2, allowed_prefixes=[]) != []


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
    path = Path("/foo/20200231_bbb_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")


def test_valid_field1():
    path = Path("/foo/20220101_WH_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) == []
    path = Path("/foo/20220101_W1_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) != []
    path = Path("/foo/20220101_Will_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) != []
    path = Path("/foo/20220101_W_c123_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field1(path) != []


def test_valid_field2():
    path = Path("/foo/20220101_WH_WH_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) == []
    path = Path("/foo/20220101_WH_W1_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) != []
    path = Path("/foo/20220101_WH_Will_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field2(path) != []
    path = Path("/foo/20220101_WH_W_444_(5)_7_8_9_10_11_12_13_14_15_16.raw")


def test_valid_field3():
    path = Path("/foo/20220101_WH_WH_project1_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field3(path) == []
    path = Path("/foo/20220101_WH_WH_project99_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field3(path) != []


def test_valid_field4():
    path = Path("/foo/20220101_WH_WH_project1_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field3(path) == []
    path = Path("/foo/20220101_WH_WH_project99_(5)_7_8_9_10_11_12_13_14_15_16.raw")
    assert valid.valid_field3(path) != []
