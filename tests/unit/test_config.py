"""Test configuration"""
import pytest

from metatlas.tools.config import Config, OutputLists, Workflow


def test_duplicate_workflow_names():
    with pytest.raises(ValueError):
        Config.parse_obj({"workflows": [{"name": "foo"}, {"name": "foo"}]})


def test_workflow_name_restrictions():
    with pytest.raises(ValueError):
        Workflow.parse_obj({"atlas_name": "x", "atlas_username": "y", "name": "Cannot have spaces"})


def test_output_list_update01():
    out = OutputLists(rt_alignment=["foo", "bar"])
    overrides = {"rt_alignment": ["zoop"]}
    out.update(overrides)
    assert out.rt_alignment == ["zoop"]
    assert out.qc_outputs == []


def test_distribute_always_values_double():
    out1 = OutputLists(always=["all1", "all2"], rt_alignment=["foo", "bar"])
    out2 = OutputLists(always=["all3", "all4"], rt_alignment=["foo2", "bar2"])
    out1.distribute_always_values()
    out2.distribute_always_values()
    assert out1.always == []
    assert out1.rt_alignment == ["foo", "bar", "all1", "all2"]
    assert out2.always == []
    assert out2.rt_alignment == ["foo2", "bar2", "all3", "all4"]


def test_config_distribute_always_values01(config):
    config.workflows[0].rt_alignment.parameters.include_groups.always = ["all-a", "all-b"]
    config.distribute_always_values()
    params = config.workflows[0].rt_alignment.parameters
    assert params.include_groups.always == []
    assert params.include_groups.rt_alignment == ["QC", "all-a", "all-b"]


def test_config_distribute_always_values02(config):
    config.workflows[1].analyses[0].parameters.exclude_groups.always = ["x", "y"]
    params = config.workflows[1].analyses[0].parameters
    assert params.exclude_groups.chromatograms == ['QC', 'NEG', 'FPS']
    config.distribute_always_values()
    assert params.exclude_groups.always == []
    assert params.exclude_groups.chromatograms == ['QC', 'NEG', 'FPS', "x", "y"]


def test_config_distribute_always_values03(config):
    config.workflows[1].analyses[0].parameters.exclude_groups.always = ["x", "y"]
    params = config.workflows[1].analyses[0].parameters
    assert params.exclude_groups.chromatograms == ['QC', 'NEG', 'FPS']
    params.distribute_always_values()
    assert params.exclude_groups.always == []
    assert params.exclude_groups.chromatograms == ['QC', 'NEG', 'FPS', "x", "y"]
