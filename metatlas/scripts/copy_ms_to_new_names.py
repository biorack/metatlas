#!/usr/bin/env python3
"""Utility script for performing batch copy+rename operations"""

import argparse
import getpass
import shutil
import sys

from pathlib import Path
from typing import List

from metatlas.tools.logging import activate_module_logging
from metatlas.tools.validate_filenames import validate_file_name

LOG_DIR = Path("/global/cfs/cdirs/m2650/copy_rename_logs")

DIRECTORY_NUM_FIELDS = 9
FILE_NUM_FIELDS = 16

EPILOG = """Order of operations:
    1. --format (first)
    2. --replace
    3. --after
    4. --before
    5. --copy   (last)

    --replace, --after --before and --copy can be supplied multiples times.
    Within an operation, they are applied in the order they are supplied on the command line from the left to the right

    Therefore, with an example set of options:
    --before 4 foo --after 4 bar --replace 1 zoop --copy 1 2 --after 8 lala --format '{"_".join([x.lower() for x in fields])}'

    They would be applied in this order:
    --format '{"_".join([x.lower() for x in fields])}'
    --replace 1 zoop
    --after 4 bar
    --after 8 lala
    --before 4 foo
    --copy 1 2
"""


class FStr:  # pylint: disable=too-few-public-methods
    """Delayed evaluation of string for f-string-like behavior on template strings"""

    def __init__(self, template, variables):
        self._template = template
        self._variables = variables

    def __repr__(self):
        return eval(f"f'{self._template}'", {}, self._variables)  # pylint: disable=eval-used


def get_args() -> argparse.Namespace:
    """Configure and parse command line arguments"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Batch copying of .raw and .mzML files to new file names",
        epilog=EPILOG,
    )
    parser.add_argument(
        "directory",
        type=lambda p: Path(p).absolute(),
        help="directory containing the files to copy to new names",
    )
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        help="python f-string that describes the filenames to generate. List 'fields' is available",
    )
    parser.add_argument(
        "-r",
        "--replace",
        type=str,
        nargs=2,
        action="append",
        default=[],
        help="takes a 0-indexed field number and the new field text",
    )
    parser.add_argument(
        "-a",
        "--after",
        type=str,
        nargs=2,
        action="append",
        default=[],
        help="takes a 0-indexed field number to insert after and new field text",
    )
    parser.add_argument(
        "-b",
        "--before",
        type=str,
        nargs=2,
        action="append",
        default=[],
        help="takes a 0-indexed field number to insert before and new field text",
    )
    parser.add_argument(
        "-c",
        "--copy",
        type=int,
        nargs=2,
        action="append",
        default=[],
        help="takes a 0-indexed source field number and 0-indexed destination field number",
    )
    parser.add_argument(
        "-k",
        "--keep-dir",
        help="copy files within the directory they already reside in",
        action="store_true",
    )
    parser.add_argument(
        "--force",
        help="Allow destination file names that do not pass validation",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--dry-run",
        help="only print out the filenames that would be created",
        action="store_true",
    )
    parser.add_argument(
        "-v", "--verbose", help="Increase output. Can be repeated for more output", action="count", default=0
    )
    return parser.parse_args()


def get_source_file_names(directory: Path) -> List[Path]:
    """Gets all .raw and .mzML (case insensitive) files in directory"""
    return list(directory.glob("*.[Rr][Aa][Ww]")) + list(directory.glob("*.[Mm][Zz][Mm][Ll]"))


def get_to_field(stem: str, field_num: int) -> str:
    """Returns the fraction of stem that goes up to field_num with "_" as the delimiter"""
    return "_".join(stem.split("_")[:field_num])


def get_dest_parent(source: Path, dest_stem: str, keep_dir: bool) -> Path:
    """Returns the directory to copy the files to"""
    if keep_dir:
        return source.parent
    return source.parent.parent / get_to_field(dest_stem, DIRECTORY_NUM_FIELDS)


def replace_field(stem: str, field_num: int, value: str) -> str:
    """Replaces a field with a constant value"""
    fields = stem.split("_")
    fields[field_num] = value
    return "_".join(fields)


def insert_field_after(stem: str, field_num: int, value: str) -> str:
    """Adds a new field with a constant value"""
    fields = stem.split("_")
    fields.insert(field_num + 1, value)
    return "_".join(fields)


def insert_field_before(stem: str, field_num: int, value: str) -> str:
    """Adds a new field with a constant value"""
    return insert_field_after(stem, field_num - 1, value)


def copy_field(stem: str, source: int, dest: int) -> str:
    """Copies a field value over an existing field"""
    fields = stem.split("_")
    fields[dest] = fields[source]
    return "_".join(fields)


def get_dest_file_name(source: Path, args: argparse.Namespace) -> Path:
    """Generates a destination path from a source path"""
    if args.format is not None:
        dest_stem = str(FStr(args.format, {"fields": source.stem.split("_")}))
    else:
        dest_stem = source.stem
    for field_num, value in args.replace:
        dest_stem = replace_field(dest_stem, int(field_num), value)
    for field_num, value in args.after:
        dest_stem = insert_field_after(dest_stem, int(field_num), value)
    for field_num, value in args.before:
        dest_stem = insert_field_before(dest_stem, int(field_num), value)
    for source_field_num, dest_field_num in args.copy:
        dest_stem = copy_field(dest_stem, source_field_num, dest_field_num)
    return (get_dest_parent(source, dest_stem, args.keep_dir) / dest_stem).with_suffix(source.suffix)


def main():
    """Main body of script"""
    args = get_args()
    logger = activate_module_logging(
        __name__,
        console_level=["WARN", "INFO", "DEBUG"][args.verbose],
        console_format="{color}{levelname:8}{reset} {message}",
        file_level="INFO",
        filename=LOG_DIR / f"{getpass.getuser( )}.log",
    )
    source_file_names = get_source_file_names(args.directory)
    files = [{"source": src, "dest": get_dest_file_name(src, args)} for src in source_file_names]
    num_invalid = sum(validate_file_name(file["dest"]) for file in files)
    if num_invalid > 0:
        logger.critical("%d destination file names failed one or more required checks.", num_invalid)
        if args.force:
            logger.info("Performing copy and rename anyways because the 'force' argument is set...")
        else:
            logger.critical(
                "Terminating without copying or renaming any files due to invalid destination file names."
            )
            sys.exit(128)
    for source, dest in files:
        if args.dry_run or args.verbose > 0:
            print(f"{source} -> {dest}")
        if not args.dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, dest)


if __name__ == "__main__":
    main()
