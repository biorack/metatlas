"""Tests that filenames meet specifications"""

import logging

from datetime import datetime
from pathlib import Path

from typing import Callable, List, Optional, Sequence

DIRECTORY_NUM_FIELDS = 9
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
    "QTOF",
    "5800",
    "QQQ",
    "QE119",
    "QE144",
    "QEHF",
    "QE139",
    "QE139-UV",
    "LTQXL",
    "4800TN",
    "4800JBEI",
]

POLARITIES = ["FPS", "POS", "NEG", "EI"]

ACQ_TYPES = ["MS1", "MSMS"]

Check = Callable[[Path], List[str]]

logger = logging.getLogger(__name__)


def validate_file_name(file_name: Path, minimal: bool = False) -> bool:
    """Run set of validation functions on a file name"""
    minimal_checks = [  # the set of checks needed for analysis to complete
        has_minimum_num_fields,
        valid_field9,
    ]
    # the union of minimal_checks and beyond_minimal_checks is set of checks derived
    # from the Standard Operating Procedure doc
    beyond_minimal_checks = [
        parent_dir_is_prefix,
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
    ]
    warn_only_checks = [
        has_exactly_num_fields,
        valid_num_subfields_field15,
    ]
    required_checks = minimal_checks + ([] if minimal else beyond_minimal_checks)
    warn_checks = (beyond_minimal_checks if minimal else []) + warn_only_checks
    return is_valid_file_name(file_name, required_checks, warn_checks)


def run_check_list(file_name: Path, checks: Sequence[Check], log_func: Callable) -> int:
    """Run checks on a file, log failures, and return number of failures"""
    num_fail = 0
    for check_func in checks:
        messages = check_func(file_name)
        if messages:
            num_fail += 1
            for text in messages:
                log_func(text)
    return num_fail


def is_valid_file_name(
    file_name: Path, required_checks: Sequence[Check], warn_only_checks: Sequence[Check]
) -> bool:
    """
    inputs:
       file_name: file_name to evaluate
       required_checks: list of checks that must pass for the file to valid
       warn_only_checks: list of checks that display a message on failure, but do not cause the file to fail

       checks are functions that return [] if valid or a list of messages if not valid
    """
    logger.info("Validating file name for: %s", file_name)
    num_warn = run_check_list(file_name, warn_only_checks, logger.warning)
    num_error = run_check_list(file_name, required_checks, logger.error)
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


def extract_batch(file_name: Path) -> str:
    """Return the batch porition of the file name"""
    parent = file_name.parent.stem
    return "_".join(parent.split("_")[:DIRECTORY_NUM_FIELDS])


def parent_dir_is_prefix(file_name: Path) -> List[str]:
    """Test that parent directory matches the start of the filename"""
    passing = file_name.stem.startswith(extract_batch(file_name))
    return [] if passing else ["Filename and parent directory do not contain the same batch fields."]


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


def num_subfields(file_name: Path, field_num: int, min_num: int, max_num: Optional[int]) -> List[str]:
    """Test for range of number of sub-fields"""
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    out = []
    sub_fields = field.split("-")
    if (field == "" and min_num > 0) or (min_num > len(sub_fields)):
        out.append("number of sub-fields is too small.")
    if max_num is not None and max_num < len(sub_fields):
        out.append("number of sub-fields is too large.")
    if out:
        prefix = f"{FIELDS[field_num]['label']} (field {field_num})"
        return [f"{prefix} {message}" for message in out]
    return []


def valid_string(
    file_name: Path,
    field_num: int,
    min_len: int,
    max_len: int,
    alpha_only: bool = False,
    alphanumeric_only: bool = False,
) -> List[str]:
    """Only test length of field"""
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    out = []
    if len(field) > max_len:
        out.append(f"contains more than {max_len} characters.")
    if len(field) < min_len:
        out.append(f"contains less than {max_len} characters.")
    if alpha_only and not field.isalpha():
        out.append("contains non-letter characters.")
    if alphanumeric_only and not field.isalnum():
        out.append("contains non-alphanumeric characters.")
    if out:
        prefix = f"{FIELDS[field_num]['label']} (field {field_num})"
        return [f"{prefix} {message}" for message in out]
    return []


def get_alpha_prefix(text: str) -> str:
    """Return text up to the first digit"""
    for i, character in enumerate(text):
        if character.isdigit():
            return text[:i]
    return text


def valid_int(
    file_name: Path,
    field_num: int,
    min_int: int,
    max_int: Optional[int] = None,
    allowed_prefixes: Optional[List[str]] = None,
) -> List[str]:
    """Only integer with optional prefix"""
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
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
        prefix = f"{FIELDS[field_num]['label']} (field {field_num})"
        return [f"{prefix} {message}" for message in out]
    return []


def valid_field0(file_name: Path) -> List[str]:
    """Date of LCMS run in YYYYMMDD format"""
    field_num = 0
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    if len(field) > 8:
        return ["Date (field 0) contains more than 8 characters."]
    if len(field) < 8:
        return ["Date (field 0) contains less than 8 characters."]
    out = []
    year = int(field[:4])
    month = int(field[4:6])
    day = int(field[6:8])
    current_time = datetime.now()
    if year < 2000:
        out.append("Year (in field 0) cannot be before 2000.")
    elif year > current_time.year:
        out.append("Year (in field 0) cannot be in the future.")
    if month < 1 or month > 12:
        out.append("Month (in field 0) must be in range 01-12.")
    if day < 1 or day > 31:
        out.append("Day (in field 0) must be in range 01-31.")
    try:
        datetime.strptime(field, "%Y%m%d")
    except ValueError:
        out.append("Day (in field 0) is not consistent with the year and month.")
    return out


def valid_field1(file_name: Path) -> List[str]:
    """Initials of LCMS user"""
    return valid_string(file_name, field_num=1, min_len=2, max_len=3, alpha_only=True)


def valid_field2(file_name: Path) -> List[str]:
    """Initials of sample submitter"""
    return valid_string(file_name, field_num=2, min_len=2, max_len=4, alpha_only=True)


def valid_field3(file_name: Path) -> List[str]:
    """Project"""
    return valid_string(file_name, field_num=3, min_len=1, max_len=8)


def valid_field4(file_name: Path) -> List[str]:
    """Experiment"""
    return valid_string(file_name, field_num=4, min_len=1, max_len=8)


def valid_field5(file_name: Path) -> List[str]:
    """Sample Set"""
    return valid_string(file_name, field_num=5, min_len=1, max_len=8)


def valid_field6(file_name: Path) -> List[str]:
    """System"""
    field_num = 6
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    if field not in SYSTEMS:
        return [f"{FIELDS[field_num]['label']} (field {field_num}) must be one of {','.join(SYSTEMS)}."]
    return []


def valid_field7(file_name: Path) -> List[str]:
    """Column-method"""
    # This could be updated to test for the known chromatography types
    return valid_string(file_name, field_num=7, min_len=1, max_len=100)


def valid_field8(file_name: Path) -> List[str]:
    """Serial"""
    return valid_string(file_name, field_num=8, min_len=1, max_len=100, alphanumeric_only=True)


def valid_field9(file_name: Path) -> List[str]:
    """Polarity"""
    field_num = 9
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    # This could be improved to test that the chromatography type matches with the polarity
    if field not in POLARITIES:
        return [f"{FIELDS[field_num]['label']} (field {field_num}) must be one of {','.join(POLARITIES)}."]
    return []


def valid_field10(file_name: Path) -> List[str]:
    """Acquisition type"""
    field_num = 10
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    sub_fields = field.split("-")
    # This could be improved to test that the chromatography type matches with the polarity
    if sub_fields[0] not in ACQ_TYPES:
        return [
            f"{FIELDS[field_num]['label']} (field {field_num}) must start with one of {','.join(ACQ_TYPES)}."
        ]
    if sub_fields[0] == "MS1" and len(sub_fields) > 1:
        return [f"{FIELDS[field_num]['label']} (field {field_num}) takes no optional sub-fields for MS1."]
    if len(sub_fields) > 2:
        return [f"{FIELDS[field_num]['label']} (field {field_num}) takes a maximum of 2 sub-fields."]
    return []


def valid_field11(file_name: Path) -> List[str]:
    """Sample number"""
    return valid_int(file_name, field_num=11, min_int=0)


def valid_field12(file_name: Path) -> List[str]:
    """Sample group"""
    return valid_string(file_name, field_num=12, min_len=1, max_len=20)


def valid_field13(file_name: Path) -> List[str]:
    """Replicate number"""
    return valid_int(file_name, field_num=13, min_int=0, max_int=999, allowed_prefixes=["Rep", "R"])


def valid_field14(file_name: Path) -> List[str]:
    """Additional parameters"""
    field_num = 14
    try:
        field = file_name.stem.split("_")[field_num]
    except IndexError:
        return [f"Field {field_num} not found"]
    out = []
    if field == "NA":
        return []
    sub_fields = field.split("-")
    if any(x == "" for x in sub_fields):
        out.append("empty sub-fields are not allowed.")
    if any(x == "NA" for x in sub_fields):
        out.append("NA only allowed as first and only sub-field.")
    if out:
        prefix = f"{FIELDS[field_num]['label']} (field {field_num})"
        return [f"{prefix} {message}" for message in out]
    return []


def valid_field15(file_name: Path) -> List[str]:
    """Sequence injection number"""
    return valid_int(file_name, field_num=15, min_int=0, max_int=999, allowed_prefixes=["Seq", "S", "Run"])


def valid_num_subfields_field15(file_name: Path) -> List[str]:
    """Additional parameters sub-field number check"""
    return num_subfields(file_name, field_num=15, min_num=1, max_num=3)
