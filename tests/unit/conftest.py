"""
per-directory pytest configuration
This makes the fixtures available to tests within this directory
"""
from glob import glob


def refactor(string: str) -> str:
    """python file path to module name converter"""
    return string.replace("/", ".").replace("\\", ".").replace(".py", "")


pytest_plugins = [
    refactor(fixture) for fixture in glob("tests/fixtures/*.py") if "__" not in fixture
]
