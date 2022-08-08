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
