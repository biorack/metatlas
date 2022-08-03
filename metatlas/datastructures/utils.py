"""Utility functions for working with metatlas data structures"""

import logging

from typing import NewType, Optional

from metatlas.datastructures import metatlas_objects as metob
from metatlas.tools.util import or_default

logger = logging.getLogger(__name__)

AtlasName = NewType("AtlasName", str)
Username = NewType("Username", str)


def get_atlas(unique_id: str, name: Optional[AtlasName] = None) -> metob.Atlas:
    """Loads atlas from database, throws error if multiple atlases match query"""
    atlases = metob.retrieve("Atlas", name=or_default(name, "%"), unique_id=unique_id, username="*")
    query_description = f"{or_default(name, '') + ' '}with unique_id '{unique_id}'."
    try:
        if len(atlases) == 0:
            raise ValueError(f"Database does not contain an atlas {query_description}")
    except ValueError as err:
        logger.exception(err)
        raise err
    try:
        if len(atlases) > 1:
            raise ValueError(f"Database contains more than one atlas {query_description}")
    except ValueError as err:
        logger.exception(err)
        raise err
    return atlases[0]
