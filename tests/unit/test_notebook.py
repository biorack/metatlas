"""Test of notebook functions"""
# pylint: disable=missing-function-docstring

import json

from metatlas.tools import notebook


def test_create_notebook_with_parameters01():
    orig_data = {
        "cells": [
            {
                "metadata": {"tags": ["parameters"]},
                "source": [
                    "# this is a comment\n",
                    "param1 = 0\n",
                    "\n",
                    "param2 = []\n",
                    'param3 = "REPLACE ME"\n',
                ],
            },
        ]
    }
    with open("test.json", "w", encoding="utf8") as out_fh:
        json.dump(orig_data, out_fh)
    notebook.create_notebook(
        "test.json", "out.json", {"param1": 1, "param2": ["foo", "bar"], "param3": "My_Exp_Name"}
    )
    with open("out.json", encoding="utf8") as in_fh:
        data = json.load(in_fh)
    assert data["cells"][0]["source"][1] == "param1 = 1\n"
    assert data["cells"][0]["source"][3] == "param2 = ['foo', 'bar']\n"
    assert data["cells"][0]["source"][4] == "param3 = 'My_Exp_Name'\n"


def test_create_notebook_with_parameters02():
    orig_data = {
        "cells": [
            {
                "metadata": {"tags": ["parameters"]},
                "source": [
                    "# this is a comment\n",
                    "param1 = 0\n",
                    "\n",
                    "param2 = []\n",
                    'param3 = "REPLACE ME"\n',
                ],
            },
        ]
    }
    with open("test.json", "w", encoding="utf8") as out_fh:
        json.dump(orig_data, out_fh)
    notebook.create_notebook("test.json", "out.json", {})
    with open("out.json", encoding="utf8") as in_fh:
        data = json.load(in_fh)
    assert data["cells"][0]["source"][1] == "param1 = 0\n"
    assert data["cells"][0]["source"][3] == "param2 = []\n"
    assert data["cells"][0]["source"][4] == 'param3 = "REPLACE ME"\n'


def test_create_notebook_with_parameters03():
    orig_data = {
        "cells": [
            {
                "metadata": {"tags": ["parameters"]},
                "source": [
                    "# this is a comment\n",
                    "param1 = 0\n",
                    "\n",
                    "param2 = []\n",
                    'param3 = "REPLACE ME"\n',
                ],
            },
        ]
    }
    with open("test.json", "w", encoding="utf8") as out_fh:
        json.dump(orig_data, out_fh)
    notebook.create_notebook("test.json", "out.json", {"param1": None, "param2": True})
    with open("out.json", encoding="utf8") as in_fh:
        data = json.load(in_fh)
    assert data["cells"][0]["source"][1] == "param1 = None\n"
    assert data["cells"][0]["source"][3] == "param2 = True\n"
