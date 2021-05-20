# /pylint: disable=line-too-long, missing-function-docstring, missing-module-docstring

import getpass
import os
import sqlite3
import pytest
from metatlas.datastructures import metatlas_objects as metob


@pytest.fixture(name="sqlite")
def fixture_sqlite(tmp_path):
    # make sure we don't accidently pollute the production MySQL DB
    assert os.environ.get("METATLAS_LOCAL") == "TRUE"
    os.chdir(tmp_path)  # don't reuse the sqlite DB
    username = getpass.getuser()
    sqlite3.connect(f"{username}_workspace.db").close()
    dummy = metob.Atlas()
    dummy.name = "this is a dummy atlas to initialize sqlite db"
    metob.store(dummy)
    # do I need to store each type of object?
