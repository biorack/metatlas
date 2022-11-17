"""Tests that filenames meet specifications"""

import logging

from datetime import datetime
from pathlib import Path

from typing import Callable, List, Optional, Sequence, Tuple

NUM_REQUIRED_DIRECTORY_FIELDS = 9
NUM_OPTIONAL_DIRECTORY_FIELDS = 1
FILE_NUM_FIELDS = 16

FIELDS = [
    {"code": "DATE", "label": "Date"},
    {"code": "NORTHENLABINITIALS", "label": "Initials of LCMS user"},
    {"code": "COLLABINITIALS", "label": "Initials of sample sumbitter"},
    {"code": "PROJ", "label": "Project"},
    {"code": "EXP", "label": "Experiment"},
    {"code": "SAMPSET", "label": "Sample set"},
    {"code": "SYSTEM", "label": "System"},
    {"code": "COLUMN-method", "label": "Chromatography and optional method parameters"},
    {"code": "SERIAL", "label": "Column serial number"},
    {"code": "POL", "label": "Polarity"},
    {"code": "ACQ", "label": "Acquisition type"},
    {"code": "SAMPLE#", "label": "Sample number"},
    {"code": "SAMPLEGROUP", "label": "Sample group"},
    {"code": "REP", "label": "Replicate number"},
    {"code": "OPTIONAL", "label": "Additonal parameters"},
    {"code": "SEQ", "label": "Sequence injection #"},
]

SYSTEMS = [
    "GCMS",
    "5800",
    "QE119",
    "QE144",
    "QEHF",
    "4800JBEI",
    "EXP120A",
    "IDX",
    "HPLC1",
    "EXP120B",
    "LTQXL",
    "4800TN",
    "QE139-UV",
    "QE139",
    "QQQ",
    "QTOF",
    "Exploris",
    "Exp1",
]

COLUMN_TYPES = ["HILICZ", "HILIC", "C18", "PGC"]

POLARITIES = ["FPS", "POS", "NEG", "EI"]

ACQ_TYPES = ["MS1", "MS2", "MS3", "MS4", "MS5", "MS6", "MS7", "MS8", "MS9", "MS10", "MSn", "MSMS"]

Check = Callable[[Path], List[str]]

logger = logging.getLogger(__name__)


def get_message_prefix(field_num) -> str:
    """Return a information about the field can can prefix an error message"""
    return f"{FIELDS[field_num]['label']} (field {field_num})"


def extract_batch(file_name: Path) -> str:
    """Return the batch porition of the file name"""
    parent = file_name.parent.stem
    return "_".join(parent.split("_")[:NUM_REQUIRED_DIRECTORY_FIELDS])


def parent_dir_is_prefix(file_name: Path) -> List[str]:
    """Test that parent directory matches the start of the filename"""
    passing = file_name.stem.startswith(extract_batch(file_name))
    return [] if passing else ["Filename and parent directory do not contain the same batch fields."]


def parent_dir_num_fields(file_name: Path) -> List[str]:
    """Test that parent directory contains a legal number of fields"""
    min_fields = NUM_REQUIRED_DIRECTORY_FIELDS
    max_fields = min_fields + NUM_OPTIONAL_DIRECTORY_FIELDS
    parent = file_name.parent.stem
    num_fields = len(parent.split("_"))
    if num_fields > max_fields:
        return [f"Parent directory contains {num_fields} fields and the maximum allowed is {max_fields}."]
    if num_fields < min_fields:
        return [f"Parent directory contains {num_fields} fields but the minimum allowed is {min_fields}."]
    return []


def has_minimum_num_fields(file_name: Path) -> List[str]:
    """Test for minimum number of fields"""
    num_fields = len(file_name.stem.split("_"))
    passing = num_fields >= FILE_NUM_FIELDS
    if passing:
        return []
    return [f"Filename contains {num_fields} fields but a minimum of {FILE_NUM_FIELDS} are required."]


def has_exactly_num_fields(file_name: Path) -> List[str]:
    """Test for exact number of fields"""
    num_fields = len(file_name.stem.split("_"))
    passing = num_fields >= FILE_NUM_FIELDS
    if passing:
        return []
    return [f"Filename contains {num_fields} fields but should have {FILE_NUM_FIELDS}."]


def field_exists(file_name: Path, field_num: int) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Output tuple is (prefix, field, error_message)"""
    try:
        prefix = get_message_prefix(field_num)
    except IndexError:
        return (None, None, f"Invalid field number: {field_num}.")
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return (None, None, f"{prefix} not found")
    return (prefix, field, None)


def valid_string(
    file_name: Path,
    field_num: int,
    min_len: int = 1,
    max_len: Optional[int] = None,
    alpha_only: bool = False,
    alphanumeric_only: bool = True,
    numeric_only: bool = False,
    allow_dashes: bool = False,
    sub_field_num: Optional[int] = None,
) -> List[str]:
    """
    Generic validation for string fields and optionally sub-fields
    if allow dashes is True, then the presence of dashes does not invalidate
    the constraints from the *_only inputs
    """
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    if sub_field_num is not None:
        prefix = f"{prefix}, sub field {sub_field_num}"
        try:
            field = field.split("-")[sub_field_num]
        except IndexError:
            return [f"{prefix} not found"]
    out = []
    if max_len is not None and len(field) > max_len:
        out.append(f"contains more than {max_len} characters.")
    if len(field) < min_len:
        out.append(f"contains less than {max_len} characters.")
    if allow_dashes:
        sub_fields = field.split("-")
        if "" in sub_fields:
            out.append("contains empty sub field.")
        field = "".join(sub_fields)
    if not field.isascii():
        out.append("contains non-ascii characters.")
    if alpha_only and not field.isalpha():
        out.append("contains non-letter characters.")
    if alphanumeric_only and not field.isalnum():
        out.append("contains non-alphanumeric characters.")
    if numeric_only and not field.isdigit():
        out.append("contains non-digit characters.")
    return [f"{prefix} {message}" for message in out]


def get_alpha_prefix(text: str) -> str:
    """Return text up to the first digit"""
    for i, character in enumerate(text):
        if character.isdigit():
            return text[:i]
    return text


def valid_int(
    file_name: Path,
    field_num: int,
    min_int: int = 0,
    max_int: Optional[int] = None,
    allowed_prefixes: Optional[List[str]] = None,
) -> List[str]:
    """Only integer with optional prefix"""
    message_prefix, field, err_msg = field_exists(file_name, field_num)
    if message_prefix is None or not field:
        return [err_msg] if err_msg else []
    out = []
    prefix = get_alpha_prefix(field)
    if prefix:
        if allowed_prefixes is None:
            out.append("does not start with a digit")
        elif prefix not in allowed_prefixes:
            out.append(f"prefix can only be one of {' ,'.join(allowed_prefixes)}")
    num_part = field[len(prefix) :]
    if not num_part.isdigit():
        out.append("cannot convert to an integer.")
    else:
        num = int(num_part)
        if max_int is not None and num > max_int:
            out.append(f"cannot be greater than {max_int}.")
        if num < min_int:
            out.append(f"cannot be less than {min_int}.")
    if out:
        return [f"{message_prefix} {message}" for message in out]
    return []


def valid_field0(file_name: Path) -> List[str]:
    """Date of LCMS run in YYYYMMDD format"""
    field_num = 0
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or field is None or field == "":
        return [err_msg] if err_msg else []
    if len(field) > 8:
        return [f"{prefix} contains more than 8 characters."]
    if len(field) < 8:
        return [f"{prefix} contains less than 8 characters."]
    out = []
    year = int(field[:4])
    month = int(field[4:6])
    day = int(field[6:8])
    current_time = datetime.now()
    if year < 2000:
        out.append(f"Year (in field {field_num}) cannot be before 2000.")
    elif year > current_time.year:
        out.append(f"Year (in field {field_num}) cannot be in the future.")
    if month < 1 or month > 12:
        out.append(f"Month (in field {field_num}) must be in range 01-12.")
    if day < 1 or day > 31:
        out.append(f"Day (in field {field_num}) must be in range 01-31.")
    try:
        datetime.strptime(field, "%Y%m%d")
    except ValueError:
        out.append(f"Day (in field {field_num}) is not consistent with the year and month.")
    return out


def valid_field1(file_name: Path) -> List[str]:
    """Initials of department generating the data"""
    return valid_string(file_name, field_num=1, min_len=2, max_len=3)


def valid_field2(file_name: Path) -> List[str]:
    """Initials of sample submitter"""
    field_num = 2
    max_len = 5
    min_len = 2
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    out = []
    if len(field) > max_len:
        out.append(f"{prefix} contains more than {max_len} characters")
    if len(field) < min_len:
        out.append(f"{prefix} contains fewer than {min_len} characters")
    sub_fields = field.split("-")
    if len(sub_fields) > 2:
        out.append(f"{prefix} contains more than 2 sub fields")
    for idx, _ in enumerate(sub_fields):
        out.extend(valid_string(file_name, field_num, min_len=1, max_len=max_len, sub_field_num=idx))
    return out


def valid_field3(file_name: Path) -> List[str]:
    """Project"""
    return valid_string(file_name, field_num=3, min_len=4, max_len=10, numeric_only=True)


def valid_field4(file_name: Path) -> List[str]:
    """Experiment"""
    return valid_string(file_name, field_num=4, min_len=3, max_len=15)


def valid_field5(file_name: Path) -> List[str]:
    """Sample Set"""
    return valid_string(file_name, field_num=5, min_len=1, max_len=15)


def valid_field6(file_name: Path) -> List[str]:
    """System"""
    field_num = 6
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    if field not in SYSTEMS:
        return [f"{prefix} must be one of {','.join(SYSTEMS)}."]
    return []


def valid_field7(file_name: Path) -> List[str]:
    """Column-method"""
    # This could be updated to test for the known chromatography types
    field_num = 7
    max_len = 15
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    out = []
    if len(field) > max_len:
        out.append(f"{prefix} contains more than {max_len} characters")
    sub_fields = field.split("-")
    if sub_fields[0] not in COLUMN_TYPES:
        out.append(f"{prefix} first sub field must be one of {','.join(COLUMN_TYPES)}.")
    if len(sub_fields) > 2:
        out.append(f"{prefix} contains more than 2 sub fields")
    if len(sub_fields) == 2:
        out.extend(valid_string(file_name, field_num, sub_field_num=1))
    return out


def valid_field8(file_name: Path) -> List[str]:
    """Serial"""
    # need to allow arbitrary number of dashes
    return valid_string(file_name, field_num=8, min_len=1, max_len=20, allow_dashes=True)


def valid_field9(file_name: Path) -> List[str]:
    """Polarity"""
    field_num = 9
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    # This could be improved to test that the chromatography type matches with the polarity
    if field not in POLARITIES:
        return [f"{prefix} must be one of {','.join(POLARITIES)}."]
    return []


def valid_field10(file_name: Path) -> List[str]:
    """Acquisition type"""
    field_num = 10
    min_len = 3
    max_len = 15
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    out = []
    if len(field) > max_len:
        out.append(f"{prefix} contains more than {max_len} characters")
    if len(field) < min_len:
        out.append(f"{prefix} contains fewer than {min_len} characters")
    sub_fields = field.split("-")
    # This could be improved to test that the chromatography type matches with the polarity
    if len(sub_fields) == 0 or sub_fields[0] not in ACQ_TYPES:
        out.append(f"{prefix} must start with one of {','.join(ACQ_TYPES)}.")
    if len(sub_fields) > 2:
        out.append(f"{prefix} contains more than 2 sub fields")
    if len(sub_fields) == 2:
        out.extend(valid_string(file_name, field_num, sub_field_num=1))
    return out


def valid_field11(file_name: Path) -> List[str]:
    """Sample number"""
    return valid_string(file_name, field_num=11, min_len=1, max_len=10, allow_dashes=True)


def valid_field12(file_name: Path) -> List[str]:
    """Sample group"""
    return valid_string(file_name, field_num=12, min_len=1, max_len=40, allow_dashes=True)


def valid_field13(file_name: Path) -> List[str]:
    """Replicate number"""
    return valid_string(file_name, field_num=13, min_len=1, max_len=8, allow_dashes=True)


def valid_field14(file_name: Path) -> List[str]:
    """Additional parameters"""
    return valid_string(file_name, field_num=14, min_len=1, max_len=40, allow_dashes=True)


def valid_field15(file_name: Path) -> List[str]:
    """Sequence injection number"""
    field_num = 15
    min_len = 1
    max_len = 7
    prefix, field, err_msg = field_exists(file_name, field_num)
    if prefix is None or not field:
        return [err_msg] if err_msg else []
    out = []
    if len(field) > max_len:
        out.append(f"{prefix} contains more than {max_len} characters")
    if len(field) < min_len:
        out.append(f"{prefix} contains fewer than {min_len} characters")
    out.extend(valid_int(file_name, field_num=15, min_int=0, max_int=9999999, allowed_prefixes=["Run"]))
    return out


# the set of checks needed for analysis to complete
minimal_checks = [
    has_minimum_num_fields,
    valid_field9,
]
# the union of minimal_checks and beyond_minimal_checks is set of checks derived
# from the Standard Operating Procedure doc
beyond_minimal_checks = [
    parent_dir_is_prefix,
    parent_dir_num_fields,
    valid_field0,
    valid_field1,
    valid_field2,
    valid_field3,
    valid_field4,
    valid_field5,
    valid_field6,
    valid_field7,
    valid_field8,
    valid_field10,
    valid_field11,
    valid_field12,
    valid_field13,
    valid_field14,
    valid_field15,
    has_exactly_num_fields,
]
warn_only_checks: List[Callable] = []


def validate_file_name(file_name: Path, minimal: bool = False) -> bool:
    """Run set of validation functions on a file name"""
    required_checks = minimal_checks + ([] if minimal else beyond_minimal_checks)
    warn_checks = (beyond_minimal_checks if minimal else []) + warn_only_checks
    return is_valid_file_name(file_name, required_checks, warn_checks)


def get_validation_messages(file_name: Path, minimal: bool = False) -> Tuple[List[str], List[str]]:
    """Runs checks on file_name and returns tuple of (warning_messages, error_messages)"""
    required_checks = minimal_checks + ([] if minimal else beyond_minimal_checks)
    warn_checks = (beyond_minimal_checks if minimal else []) + warn_only_checks
    return (run_check_list(file_name, warn_checks), run_check_list(file_name, required_checks))


def run_check_list(file_name: Path, checks: Sequence[Check]) -> List[str]:
    """Run checks on a file and return list of failure messages"""
    return [message for check_fun in checks for message in check_fun(file_name)]


def is_valid_file_name(
    file_name: Path, required_checks: Sequence[Check], warn_checks: Sequence[Check]
) -> bool:
    """
    inputs:
       file_name: file_name to evaluate
       required_checks: list of checks that must pass for the file to valid
       warn_checks: list of checks that display a message on failure, but do not cause the file to fail

       checks are functions that return [] if valid or a list of messages if not valid
    """
    logger.info("Validating file name for: %s", file_name)
    warnings = run_check_list(file_name, warn_checks)
    for warn in warnings:
        logger.warning(warn)
    errors = run_check_list(file_name, required_checks)
    for error in errors:
        logger.error(error)
    num_warn = len(warnings)
    num_error = len(errors)
    if num_error == 0:
        if num_warn == 0:
            logger.info("Passed filename validation: %s", file_name)
        else:
            logger.info("Passed filename validation, but had %d warnings: %s", num_warn, file_name)
        return True
    logger.info(
        "Failed filename validation with %d errors and %d warnings: %s", num_error, num_warn, file_name
    )
    return False
