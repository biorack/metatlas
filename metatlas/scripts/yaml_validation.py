#!/usr/bin/env python
"""YAML validation with extra check for issues with key"""

import itertools
import sys

from typing import Any, List

import numpy as np
import yaml


def flatten(list_of_lists: List[Any]) -> List[Any]:
    """flatten nested lists of arbitrary depth"""
    return list(itertools.chain.from_iterable(list_of_lists))


def get_keys(in_data: Any) -> List:
    """Recursively get all keys used within a data structure"""
    if np.isscalar(in_data) or in_data is None:
        return []
    try:
        return list(in_data.keys()) + flatten([get_keys(v) for v in in_data.values()])
    except AttributeError:
        # some sort of list like iterable
        return flatten([get_keys(x) for x in in_data])


def check_for_colon(in_keys: List[str]) -> None:
    """die with error if any of in_keys contains a ':'"""
    all_valid = True
    for k in in_keys:
        if ':' in k:
            print(
                f"ERROR: Bad key '{k}'. Either quote the key or add a space after the ':'.",
                file=sys.stderr
            )
            all_valid = False
    if not all_valid:
        sys.exit(128)


if __name__ == "__main__":
    try:
        data = yaml.safe_load(sys.stdin)
    except yaml.parser.ParserError:
        print("ERROR: Could not parse YAML.", file=sys.stderr)
        sys.exit(128)
    keys = get_keys(data)
    check_for_colon(keys)
