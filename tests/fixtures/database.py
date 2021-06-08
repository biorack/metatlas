# /pylint: disable=line-too-long, missing-function-docstring, missing-module-docstring

import getpass
import os
import sqlite3
import pytest
from metatlas.datastructures import metatlas_objects as metob


@pytest.fixture(name="sqlite")
def fixture_sqlite(tmp_path, monkeypatch):
    # make sure we don't accidently pollute the production MySQL DB
    monkeypatch.setenv("METATLAS_LOCAL", "TRUE")
    db_path = tmp_path / "workspace.db"
    monkeypatch.setenv("METATLAS_SQLITE", str(db_path))
    sqlite3.connect(db_path).close()
    dummy = metob.Atlas()
    dummy.name = "this is a dummy atlas to initialize sqlite db"
    metob.store(dummy)
