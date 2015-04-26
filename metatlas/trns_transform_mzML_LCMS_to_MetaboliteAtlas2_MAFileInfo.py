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
    from . import kbase_utils as script_utils

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
              input_mapping=None, name=None, polarity='', atlases=None,
              group='', inclusion_order='', normalization_factor='',
              retention_correction='',
              level=logging.INFO, logger=None):
    """
    Converts mzML file to MetaboliteAtlas2_MAFileInfo json string.

    Args:
        shock_service_url: A url for the KBase SHOCK service.
        handle_service_url: A url for the KBase Handle Service.
        output_file_name: A file name where the output JSON string should be
                          stored.
                          If the output file name is not specified the name
                          will default
                          to the name of the input file appended with
                           '_finfo'.
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
        atlases: List of MetaboliteAtlas atlas IDs.
        name: Name of the file, optional.  Defaults to the file name.
        polarity: Run polarity.
        group: Run group.
        inclusion_order: Run inclusion_order.
        retention_correction: Run retention_correction.
        normalization_factor: Run normalization factor.

    Returns:
        JSON files on disk that can be saved as a KBase workspace objects.

    Authors:
        Steven Silvester
    """

    if logger is None:
        logger = script_utils.stderrlogger(__file__)

    logger.info("Starting conversion of mzML to MetaboliteAtlas2.RunSet")
    token = os.environ.get('KB_AUTH_TOKEN')

    if not working_directory or not os.path.isdir(working_directory):
        raise Exception("The working directory {0} is not a valid directory!"
                        .format(working_directory))

    logger.info("Scanning for mzML files.")

    valid_extensions = [".mzML"]

    files = os.listdir(input_directory)
    mzml_files = [x for x in files
                             if os.path.splitext(x)[-1] in valid_extensions]

    assert len(mzml_files) != 0

    logger.info("Found {0} files".format(len(mzml_files)))

    for fname in mzml_files:
        path = os.path.join(input_directory, fname)

        if not os.path.isfile(path):
            raise Exception("The input file name {0} is not a file!"
                            .format(path))

        hdf_file = mzml_loader.mzml_to_hdf(path)

        if shock_service_url:
            shock_info = script_utils.upload_file_to_shock(logger,
                    shock_service_url, hdf_file, token=token)

        run_info = dict()
        run_info['name'] = name or fname.replace('.mzML', '')
        run_info['atlases'] = atlases or []
        run_info['polarity'] = polarity
        run_info['group'] = group
        run_info['inclusion_order'] = inclusion_order
        run_info['normalization_factor'] = normalization_factor
        run_info['retention_correction'] = retention_correction

        if shock_service_url:
            run_info["run_file_id"] = shock_info["id"]
        else:
            run_info['run_file_id'] = 'dummy_shock_id'

        output_file_name = fname.replace('.mzML', '_finfo.json')

        # This generates the json for the object
        objectString = simplejson.dumps(run_info, sort_keys=True, indent=4)

        output_file_path = os.path.join(working_directory, output_file_name)
        with open(output_file_path, "w") as outFile:
            outFile.write(objectString)

    logger.info("Conversion completed.")


def main():
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

    # custom arguments specific to this uploader
    parser.add_argument('--polarity',
                        help=script_details["Args"]["polarity"],
                        action='store', type=str, default='', required=False)
    parser.add_argument('--group',
                        help=script_details["Args"]["group"],
                        action='store', type=str, default='', required=False)
    parser.add_argument('--inclusion_order',
                        help=script_details["Args"]["inclusion_order"],
                        action='store', type=str, default='', required=False)
    parser.add_argument('--retention_correction',
                        help=script_details["Args"]["retention_correction"],
                        action='store', type=str, default='', required=False)
    parser.add_argument('--atlases',
                        help=script_details["Args"]["atlases"],
                        action='store', type=str, nargs='?', default=None,
                        required=False)
    parser.add_argument('--name',
                        help=script_details["Args"]["name"],
                        action='store', type=str, default='', required=False)
    parser.add_argument('--normalization_factor',
                        help=script_details["Args"]["normalization_factor"],
                        action='store', type=str, default='', required=False)

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


# called only if script is run from command line
if __name__ == "__main__":  # pragma: no cover
    main()
