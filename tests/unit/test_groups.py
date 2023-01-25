""" Group creation tests """

import pytest

from metatlas.datastructures import groups

VOCAB = ["QC", "InjBl", "ISTD"]


def test_get_groups_and_runs(lcmsrun):
    groups_list, runs_list = groups.get_groups_and_runs('exec', VOCAB, [], [], [], [], [lcmsrun])
    assert len(runs_list) == 1
    assert len(groups_list) == 1
    assert groups_list[0].name == '20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_exec_Cone-S1'


def test_filter_lcmsruns01(lcmsrun):
    with pytest.raises(ValueError):
        groups.filter_lcmsruns([lcmsrun], ['NotFound'], [])


def test_filter_lcmsruns02(lcmsrun):
    runs = groups.filter_lcmsruns([lcmsrun], ['Cone'], [])
    assert len(runs) == 1


def test_filter_lcmsruns03(lcmsrun):
    runs = groups.filter_lcmsruns([lcmsrun], [], ['foobar'])
    assert len(runs) == 1


def test_group_name(lcmsrun):
    out_dict = groups.group_name(lcmsrun.name, 'exec', VOCAB)
    assert out_dict["short_name"] == 'POS_Cone-S1'


def test_create_groups(lcmsrun):
    groups_list = groups.create_groups('exec', VOCAB, [lcmsrun])
    assert len(groups_list) == 1
    assert groups_list[0].name == "20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_exec_Cone-S1"


def test_filter_groups01(group):
    out = groups.filter_groups([group], [], [])
    assert out == [group]


def test_filter_groups02(group):
    out = groups.filter_groups([group], ['Cone'], [])
    assert out == [group]


def test_filter_groups03(group):
    out = groups.filter_groups([group], [], ['Cone'])
    assert len(out) == 0


def get_lcmsruns(group):
    runs = groups.get_lcmsruns([group])
    assert len(runs) == 1
