"""Utility functions for working with metatlas data structures"""

import logging

from typing import NewType, Optional

from metatlas.datastructures import metatlas_objects as metob

logger = logging.getLogger(__name__)

AtlasName = NewType("AtlasName", str)
Username = NewType("Username", str)


def get_atlas(
    name: Optional[AtlasName] = None, username: Optional[Username] = None, unique_id: Optional[str] = None
) -> metob.Atlas:
    """Loads atlas from database, throws error if multiple atlases match query"""
    args = dict(name=name, username=username, unique_id=unique_id)
    args_no_none = {k: v for k, v in args.items() if v is not None}
    atlases = metob.retrieve("Atlas", **args_no_none)
    try:
        if len(atlases) == 0:
            raise ValueError(f"Database does not contain an atlas with {args_no_none}.")
    except ValueError as err:
        logger.exception(err)
        raise err
    try:
        if len(atlases) > 1:
            raise ValueError(f"Database contains more than one atlas with {args_no_none}.")
    except ValueError as err:
        logger.exception(err)
        raise err
    return atlases[0]
