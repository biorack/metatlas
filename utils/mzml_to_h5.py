#!/usr/bin/env python3

import argparse
import logging
import sys

from metatlas.io.file_converter import mzml_to_h5_and_add_to_db

parser = argparse.ArgumentParser()
parser.add_argument("mzml_file", help="mzML mass spec file to convert")
args = parser.parse_args()

logging.basicConfig(format="%(levelname)s, %(message)s", level=logging.INFO)
if mzml_to_h5_and_add_to_db(args.mzml_file):
    sys.exit(0)
sys.exit(1)
