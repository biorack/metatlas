"""Utility functions for working with metatlas data structures"""

import logging

from typing import NewType

from metatlas.datastructures import metatlas_objects as metob

logger = logging.getLogger(__name__)

AtlasName = NewType("AtlasName", str)
Username = NewType("Username", str)


def get_atlas(name: AtlasName, username: Username) -> metob.Atlas:
    """Load atlas from database"""
    atlases = metob.retrieve("Atlas", name=name, username=username)
    try:
        if len(atlases) == 0:
            raise ValueError(f"Database does not contain an atlas {name} owned by {username}.")
    except ValueError as err:
        logger.exception(err)
        raise err
    try:
        if len(atlases) > 1:
            raise ValueError(f"Database contains more than one atlas {name} owned by {username}.")
    except ValueError as err:
        logger.exception(err)
        raise err
    return atlases[0]
