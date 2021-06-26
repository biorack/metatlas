"""Test of environment setup functions"""
# pylint: disable=missing-function-docstring

import json

from metatlas.tools import environment


def test_create_notebook_with_parameters01():
    orig_data = {
        "cells": [
            None,
            {
                "source": [
                    "# this is a comment\n",
                    "param1 = 0\n",
                    "\n",
                    "param2 = []\n",
                    'param3 = "REPLACE ME"\n',
                ]
            },
        ]
    }
    with open("test.json", "w") as out_fh:
        json.dump(orig_data, out_fh)
    environment.create_notebook_with_parameters(
        "test.json", "out.json", {"param1": 1, "param2": ["foo", "bar"], "param3": "My_Exp_Name"}
    )
    with open("out.json") as in_fh:
        data = json.load(in_fh)
    assert data["cells"][1]["source"][1] == "param1 = 1\n"
    assert data["cells"][1]["source"][3] == "param2 = ['foo', 'bar']\n"
    assert data["cells"][1]["source"][4] == "param3 = 'My_Exp_Name'\n"


def test_create_notebook_with_parameters02():
    orig_data = {
        "cells": [
            None,
            {
                "source": [
                    "# this is a comment\n",
                    "param1 = 0\n",
                    "\n",
                    "param2 = []\n",
                    'param3 = "REPLACE ME"\n',
                ]
            },
        ]
    }
    with open("test.json", "w") as out_fh:
        json.dump(orig_data, out_fh)
    environment.create_notebook_with_parameters("test.json", "out.json", {})
    with open("out.json") as in_fh:
        data = json.load(in_fh)
    assert data["cells"][1]["source"][1] == "param1 = 0\n"
    assert data["cells"][1]["source"][3] == "param2 = []\n"
    assert data["cells"][1]["source"][4] == 'param3 = "REPLACE ME"\n'


def test_create_notebook_with_parameters03():
    orig_data = {
        "cells": [None, {"source": ["# this is a comment\n", "param1 = True\n", "\n", "param2 = None\n"]}]
    }
    with open("test.json", "w") as out_fh:
        json.dump(orig_data, out_fh)
    environment.create_notebook_with_parameters("test.json", "out.json", {"param1": None, "param2": True})
    with open("out.json") as in_fh:
        data = json.load(in_fh)
    assert data["cells"][1]["source"][1] == "param1 = None\n"
    assert data["cells"][1]["source"][3] == "param2 = True\n"
