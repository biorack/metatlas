#!/usr/bin/env python

# standard library imports
import os
import sys
import logging
import hashlib

# 3rd party imports
import simplejson

# KBase imports
try:
    import biokbase.Transform.script_utils as script_utils
except ImportError:
    script_utils = None

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
              input_mapping=None, fasta_reference_only=False,
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
                           '_runset'
        input_directory: The directory where files will be read from.
        working_directory: The directory the resulting json file will be
                           written to.
        shock_id: Shock id for the fasta file if it already exists in shock
        handle_id: Handle id for the fasta file if it already exists as a
                    handle
        input_mapping: JSON string mapping of input files to expected types.
                       If you don't get this you need to scan the input
                       directory and look for your files.
        level: Logging level, defaults to logging.INFO.

    Returns:
        JSON file on disk that can be saved as a KBase workspace object.

    Author:
        Steven Silvester
    """

    if logger is None:
        logger = script_utils.stderrlogger(__file__)

    logger.info("Starting conversion of mzML to MetaboliteAtlas2.RunSet")
    token = os.environ.get('KB_AUTH_TOKEN')

    if input_mapping is None:
        logger.info("Scanning for mzML files.")

        valid_extensions = [".mzML"]

        files = os.listdir(input_directory)
        mzml_files = [x for x in files
                             if os.path.splitext(x)[-1] in valid_extensions]

        assert len(mzml_files) != 0

        logger.info("Found {0}".format(str(mzml_files)))

        input_file_name = os.path.join(input_directory, files[0])

        if len(mzml_files) > 1:
            logger.warning("Not sure how to handle multiple mzML files in this"
                           " context. Using {0}".format(input_file_name))
    else:
        join = os.path.join
        input_file_name = join(join(input_directory, "mzML.LCMS"),
                               simplejson.loads(input_mapping)["mzML.LCMS"])

    logger.info("Building Object.")

    if not os.path.isfile(input_file_name):
        raise Exception("The input file name {0} is not a file!"
                        .format(input_file_name))

    if not os.path.isdir(args.working_directory):
        raise Exception("The working directory {0} is not a valid directory!"
                        .format(working_directory))

    fasta_filesize = os.stat(input_file_name).st_size
    if fasta_filesize > 1000000000:
        # Fasta file too large to save sequences into the ContigSet object.
        contigset_warn = """The FASTA input file seems to be too large. A ContigSet
                            object will be created without sequences, but will
                            contain a reference to the file."""
        logger.warning(contigset_warn)

    # TODO: parse the input file here
    run_dict = dict()
    run_dict["md5"] = hashlib.md5()
    run_dict["id"] = output_file_name
    run_dict["name"] = output_file_name
    run_dict["source"] = "KBase"
    run_dict["source_id"] = os.path.basename(input_file_name)

    if output_file_name is None:
        # default to input file name minus file extenstion adding "_run"
        # to the end
        base = os.path.basename(input_file_name)
        output_file_name = "{0}_run.json".format(os.path.splitext(base)[0])

    if shock_id is None:
        shock_info = script_utils.upload_file_to_shock(logger,
                shock_service_url, input_file_name, token=token)
        shock_id = shock_info["id"]

    run_dict["Run_data_ref"] = shock_id

    # For future development if the type is updated to the handle_reference
    # instead of a shock_reference

    # This generates the json for the object
    objectString = simplejson.dumps(run_dict, sort_keys=True, indent=4)

    logger.info("Run data structure creation completed.  Writing out JSON.")

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
