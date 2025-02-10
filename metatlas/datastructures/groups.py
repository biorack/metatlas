""" Group generation """

import logging

from typing import cast, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from metatlas.datastructures.id_types import (
    FileMatchList,
    GroupList,
    GroupMatchList,
    LcmsRunDict,
    LcmsRunsList,
)
import metatlas.datastructures.metatlas_objects as metob
import metatlas.plots.dill2plots as dp

from metatlas.tools.util import or_default


logger = logging.getLogger(__name__)


# pylint: disable=too-many-arguments
def get_groups_and_runs(
    execution: str,
    groups_controlled_vocab: Sequence[str],
    include_lcmsruns: FileMatchList,
    exclude_lcmsruns: FileMatchList,
    include_groups: GroupMatchList,
    exclude_groups: GroupMatchList,
    all_lcmsruns: LcmsRunsList,
    all_groups: Optional[GroupList] = None,
) -> Tuple[LcmsRunsList, GroupList]:
    """Generate list of currently used LCMS runs and groups"""
    filtered_lcmsruns = filter_lcmsruns(all_lcmsruns, include_lcmsruns, exclude_lcmsruns)
    if all_groups is None:
        all_groups = create_groups(execution, groups_controlled_vocab, filtered_lcmsruns)
    groups = filter_groups(all_groups, include_groups, exclude_groups)
    return (groups, get_lcmsruns(groups))


def filter_lcmsruns(
    all_lcmsruns: LcmsRunsList, include: FileMatchList, exclude: FileMatchList
) -> LcmsRunsList:
    """filter lcmsruns by include and exclude lists"""
    if include:
        post_include = [r for r in all_lcmsruns if any(map(r.name.__contains__, include))]
        logger.debug(
            "Filtered out %d LCMS runs for not matching within include containing: %s",
            len(all_lcmsruns) - len(post_include),
            include,
        )
    else:
        post_include = all_lcmsruns
    if exclude:
        post_exclude = [r for r in post_include if not any(map(r.name.__contains__, exclude))]
        logger.debug(
            "Filtered out %d LCMS runs for matching within exclude containing: %s",
            len(post_include) - len(post_exclude),
            exclude,
        )
    else:
        post_exclude = post_include
    #for run in post_exclude:
    #    logger.debug("Run: %s", run.name)
    logger.debug("Runs: %s", len(post_exclude))
    logger.debug("After filtering, %s LCMS output files remain.", len(post_exclude))
    try:
        if len(post_exclude) == 0:
            raise ValueError("At least 1 LCMS run is required for analysis.")
    except ValueError as err:
        logger.exception(err)
        raise err
    return post_exclude


def group_name(base_filename: str, execution: str, groups_controlled_vocab: Sequence[str]) -> Dict[str, str]:
    """Returns dict with keys group and short_name corresponding to base_filename"""
    tokens = base_filename.split("_")
    prefix = "_".join(tokens[:11])
    empty_list: List[str] = []  # to define type for s below
    indices = [
        i
        for i, s in enumerate(or_default(groups_controlled_vocab, empty_list))
        if s.lower() in base_filename.lower()
    ]
    suffix = groups_controlled_vocab[indices[0]].lstrip("_") if indices else tokens[12]
    full_name = f"{prefix}_{execution}_{suffix}"
    short_name = f"{tokens[9]}_{suffix}"  # Prepending POL to short_name
    return {"group": full_name, "short_name": short_name}


def create_groups(
    execution: str, groups_controlled_vocab: Sequence[str], lcmsruns: LcmsRunsList
) -> GroupList:
    """Generate the full set of groups for lcmsruns"""
    files_dict: Dict[str, LcmsRunDict] = {}
    for lcms_file in lcmsruns:
        base_name: str = lcms_file.name.split(".")[0]
        files_dict[base_name] = cast(
            LcmsRunDict, {"object": lcms_file, **group_name(base_name, execution, groups_controlled_vocab)}
        )
    groups_df = pd.DataFrame(files_dict).T
    groups_df.drop(columns=["object"], inplace=True)
    groups_df.index.name = "filename"
    groups_df.reset_index(inplace=True)
    unique_groups_df = groups_df[["group", "short_name"]].drop_duplicates()
    return [
        metob.Group(
            name=values["group"],
            short_name=values["short_name"],
            items=[
                file_value["object"]
                for file_value in files_dict.values()
                if file_value["group"] == values["group"]
            ],
        )
        for values in unique_groups_df.to_dict("index").values()
    ]


def filter_groups(groups: GroupList, include: GroupMatchList, exclude: GroupMatchList) -> GroupList:
    """filter and sorts groups"""
    recent = dp.filter_metatlas_objects_to_most_recent(groups, "name")
    post_include = dp.filter_metatlas_objects_by_list(recent, "name", include)
    if include:
        logger.debug(
            "Filtered out %d groups for not matching within include containing: %s",
            len(recent) - len(post_include),
            include,
        )
    post_exclude = dp.remove_metatlas_objects_by_list(post_include, "name", exclude)
    if exclude:
        logger.debug(
            "Filtered out %d groups for matching within exclude containing: %s",
            len(post_include) - len(post_exclude),
            exclude,
        )
    not_empty = sorted(dp.filter_empty_metatlas_objects(post_exclude, "items"), key=lambda x: x.name)
    if include or exclude:
        logger.debug(
            "Filtered out %d groups because no LCMS runs remained in them",
            len(post_exclude) - len(not_empty),
        )
    return not_empty


def get_lcmsruns(groups: GroupList) -> LcmsRunsList:
    return list({run for group in groups for run in group.items})
