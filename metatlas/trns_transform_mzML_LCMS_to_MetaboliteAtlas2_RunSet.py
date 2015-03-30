#!/usr/bin/env python

# standard library imports
import os
import sys
import logging
import datetime

# 3rd party imports
import simplejson

# KBase imports
try:
    import biokbase.Transform.script_utils as script_utils
except ImportError:
    from . import script_utils

# Metatlas Imports
from metatlas import mzml_loader


if sys.version.startswith('3'):
    unicode = str


# transformation method that can be called if this module is imported
# Note the logger has different levels it could be run.
# See: https://docs.python.org/2/library/logging.html#logging-levels
#
# The default level is set to INFO which includes everything except DEBUG


def transform(shock_service_url=None, handle_service_url=None,
              output_file_name=None, input_directory=None,
              working_directory=None, shock_id=None, handle_id=None,
              input_mapping=None,
              level=logging.INFO, logger=None):
    """
    Converts mzML file to MetaboliteAtlas2_RunSet json string.

    Parameters
    ----------
        shock_service_url: A url for the KBase SHOCK service.
        handle_service_url: A url for the KBase Handle Service.
        output_file_name: A file name where the output JSON string should be
                          stored.
                          If the output file name is not specified the name
                          will default
                          to the name of the input file appended with
                           '_finfo
        input_directory: The directory where files will be read from.
        working_directory: The directory the resulting json file will be
                           written to.
        shock_id: Shock id for the hdf file if it already exists in shock
        handle_id: Handle id for the hdf file if it already exists as a
                    handle
        input_mapping: JSON string mapping of input files to expected types.
                       If you don't get this you need to scan the input
                       directory and look for your files.
        level: Logging level, defaults to logging.INFO.

    Returns:
        JSON files on disk that can be saved as a KBase workspace objects.

    Author:
        Steven Silvester
    """

    if logger is None:
        logger = script_utils.stderrlogger(__file__)

    logger.info("Starting conversion of mzML to MetaboliteAtlas2.RunSet")
    token = os.environ.get('KB_AUTH_TOKEN')

    if not os.path.isdir(args.working_directory):
        raise Exception("The working directory {0} is not a valid directory!"
                        .format(working_directory))

    logger.info("Scanning for mzML files.")

    valid_extensions = [".mzML"]

    files = os.listdir(input_directory)
    mzml_files = [x for x in files
                         if os.path.splitext(x)[-1] in valid_extensions]

    assert len(mzml_files) != 0

    logger.info("Found {0} files".format(len(mzml_files)))

    for fname in files:
        path = os.path.join(input_directory, fname)

        if not os.path.isfile(path):
            raise Exception("The input file name {0} is not a file!"
                            .format(path))

        hdf_file = mzml_loader(path)

        shock_info = script_utils.upload_file_to_shock(logger,
                shock_service_url, hdf_file, token=token)

        run_info = dict()
        run_info['name'] = fname.replace('.mzML', '')
        run_info['polarity'] = 0
        run_info['group'] = ''
        run_info['inclusion_order'] = ''
        run_info['normalization_factor'] = ''
        run_info['retention_correction'] = ''
        run_info['creator'] = ''
        run_info['creation_date'] = datetime.datetime.utcnow()
        run_info["data"] = shock_info["id"]

        output_file_name = fname.replace('.mzML', '_finfo')

        # This generates the json for the object
        objectString = simplejson.dumps(run_info, sort_keys=True, indent=4)

        output_file_path = os.path.join(working_directory, output_file_name)
        with open(output_file_path, "w") as outFile:
            outFile.write(objectString)

    logger.info("Conversion completed.")


# called only if script is run from command line
if __name__ == "__main__":
    script_details = script_utils.parse_docs(transform.__doc__)

    import argparse

    parser = argparse.ArgumentParser(prog=__file__,
                                     description=script_details["Description"],
                                     epilog=script_details["Authors"])

    parser.add_argument('--shock_service_url',
                        help=script_details["Args"]["shock_service_url"],
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--handle_service_url',
                        help=script_details["Args"]["handle_service_url"],
                        action='store', type=str, nargs='?', default=None,
                        required=False)
    parser.add_argument('--input_directory',
                        help=script_details["Args"]["input_directory"],
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--working_directory',
                        help=script_details["Args"]["working_directory"],
                        action='store', type=str, nargs='?', required=True)
    parser.add_argument('--output_file_name',
                        help=script_details["Args"]["output_file_name"],
                        action='store', type=str, nargs='?', default=None,
                        required=False)
    parser.add_argument('--shock_id',
                        help=script_details["Args"]["shock_id"],
                        action='store', type=str, nargs='?', default=None,
                        required=False)
    parser.add_argument('--handle_id',
                        help=script_details["Args"]["handle_id"],
                        action='store', type=str, nargs='?', default=None,
                        required=False)

    parser.add_argument('--input_mapping',
                        help=script_details["Args"]["input_mapping"],
                        action='store', type=unicode, nargs='?', default=None,
                        required=False)

    args, unknown = parser.parse_known_args()

    logger = script_utils.stderrlogger(__file__)

    logger.debug(args)
    try:
        transform(shock_service_url=args.shock_service_url,
                  handle_service_url=args.handle_service_url,
                  output_file_name=args.output_file_name,
                  input_directory=args.input_directory,
                  working_directory=args.working_directory,
                  shock_id=args.shock_id,
                  handle_id=args.handle_id,
                  input_mapping=args.input_mapping,
                  logger=logger)
    except Exception as e:
        logger.exception(e)
        sys.exit(1)

    sys.exit(0)
