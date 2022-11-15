#!/usr/bin/env python3

import argparse
from pathlib import Path

from metatlas.tools.config import load_config

DEFAULT_FILE = Path("/global/cfs/cdirs/m2650/targeted_analysis/metatlas_config.yaml")

parser = argparse.ArgumentParser()
parser.add_argument("config_file", type=Path, nargs='?', default=DEFAULT_FILE, help="Metatlas configuration file to test")
args = parser.parse_args()
load_config(args.config_file)
print("Passed")
