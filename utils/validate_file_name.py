#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

from metatlas.tools.validate_filenames import validate_file_name

logging.basicConfig(format='%(levelname)s, %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('msms_file', help='mass spec file name to validate')
args = parser.parse_args()

if not validate_file_name(Path(args.msms_file), minimal=True):
    logging.error('Invalid file name: %s', args.msms_file)
    sys.exit(128)
logging.debug('Valid file name: %s', args.msms_file)
sys.exit(0)
