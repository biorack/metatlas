"""Test configuration"""
import pytest

from metatlas.tools.config import Config, Workflow


def test_duplicate_workflow_names():
    with pytest.raises(ValueError):
        Config.parse_obj({"workflows": [{"name": "foo"}, {"name": "foo"}]})


def test_workflow_name_restrictions():
    with pytest.raises(ValueError):
        Workflow.parse_obj({"atlas_name": "x", "atlas_username": "y", "name": "Cannot have spaces"})
