#!/usr/bin/env python

import sys, os
import subprocess
import glob
import json
import argparse
import datetime
import time
import zipfile
import re
from distutils import dir_util, file_util
import pwd
import grp
from ast import literal_eval

import numpy as np
import pandas as pd

sys.path.insert(0, "/global/homes/d/dgct/Repos/metatlas/")
import metatlas.metatlas_objects as metob
from metatlas.helpers import dill2plots as dp

SCRATCH = os.environ["SCRATCH"]
TOP_FOLDER = ""
uid = pwd.getpwnam(os.environ["USER"]).pw_uid
gid = grp.getgrnam("genome").gr_gid


def get_files_to_submit(jamo_submission, args):
    """
    Finds and returns all files based on jamo_submission.type, jamo_submission.groups,
    and jamo_submission.file_filters. Files that are not raw data are copied into the
    user's scratch folder and given proper permissions and analyses are renamed using
    soft links and jamo_submission.prepend to assure unique names between analyses.
    """

    if jamo_submission.type in ["raw_data", "pactolus", "spectral_hits"]:

        if jamo_submission.is_group:
            files = set()

            for g in jamo_submission.groups:
                groups = dp.select_groups_for_analysis(
                    name=g,
                    do_print=False,
                    most_recent=True,
                    remove_empty=True,
                    include_list=[],
                    exclude_list=jamo_submission.file_filters,
                )
                for g in groups:
                    for f in g.items:
                        files.add(f.mzml_file)
            print len(files)

        else:
            assert all(
                os.path.isfile(path) or os.path.isdir(path)
                for path in jamo_submission.groups
            )

            files = set(path for path in jamo_submission.groups if os.path.isfile(path))
            for folder in [
                path for path in jamo_submission.groups if os.path.isdir(path)
            ]:
                for f in os.listdir(folder):
                    if f.lower().endswith(".mzml"):
                        if not any(
                            filt.lower() in f.lower()
                            for filt in jamo_submission.file_filters
                        ):
                            files.add(os.path.realpath(os.path.join(folder, f)))

        if jamo_submission.type == "raw_data":
            files |= set(
                os.path.splitext(f)[0] + ".h5"
                for f in files
                if os.path.isfile(os.path.splitext(f)[0] + ".h5")
            )
            print len(files)
        if jamo_submission.type == "pactolus":
            files = set(
                os.path.splitext(f)[0] + ".pactolus.gz"
                for f in files
                if os.path.isfile(os.path.splitext(f)[0] + ".pactolus.gz")
            )
        if jamo_submission.type == "spectral_hits":
            files = set(
                os.path.splitext(f)[0] + "_spectral-hits.tab.gz"
                for f in files
                if os.path.isfile(os.path.splitext(f)[0] + "_spectral-hits.tab.gz")
            )

    if jamo_submission.type in ["sop", "readme"]:
        assert all(
            os.path.isfile(path) or os.path.isdir(path)
            for path in jamo_submission.groups
        )

        files = set()

        for f in [path for path in jamo_submission.groups if os.path.isfile(path)]:
            new_f = os.path.join(SCRATCH, os.path.basename(f))

            file_util.copy_file(f, new_f, update=1)
            os.chown(new_f, uid, gid)
            os.chmod(new_f, 0o640)

            files.add(new_f)

        for folder in [path for path in jamo_submission.groups if os.path.isdir(path)]:

            new_folder = os.path.join(
                SCRATCH, os.path.basename(os.path.dirname(folder))
            )

            dir_util.copy_tree(folder, new_folder, update=1)
            os.chown(new_folder, uid, gid)
            os.chmod(new_folder, 0o750)

            for root, sub_dirs, dir_files in os.walk(new_folder):
                for sub_dir in sub_dirs:
                    os.chown(os.path.join(new_folder, root, sub_dir), uid, gid)
                    os.chmod(os.path.join(new_folder, root, sub_dir), 0o750)

                for f in dir_files:
                    os.chown(os.path.join(new_folder, root, f), uid, gid)
                    os.chmod(os.path.join(new_folder, root, f), 0o640)

                    files.add(os.path.join(new_folder, root, f))

    if jamo_submission.type in ["untargeted_analysis", "targeted_analysis"]:

        assert all(
            os.path.isfile(path) or os.path.isdir(path)
            for path in jamo_submission.groups
        )

        files = set()

        for f in [path for path in jamo_submission.groups if os.path.isfile(path)]:
            new_f = os.path.join(SCRATCH, os.path.basename(f))

            file_util.copy_file(f, new_f, update=1)
            os.chown(new_f, uid, gid)
            os.chmod(new_f, 0o640)

            if not os.path.islink(new_f):

                while True:
                    name, ext = os.path.splitext(new_f)
                    link_f = (
                        name
                        + "-"
                        + re.sub(
                            r"[^0-9]",
                            "",
                            str(datetime.datetime.now().replace(microsecond=0)),
                        )
                        + ext
                    )

                    if os.path.basename(link_f).lower() in [
                        os.path.basename(f).lower() for f in files
                    ]:
                        time.sleep(1)
                        continue
                    else:
                        os.symlink(
                            os.path.join(SCRATCH, new_f), os.path.join(SCRATCH, link_f)
                        )
                        new_f = link_f
                        break

                files.add(os.path.join(SCRATCH, new_f))

        for folder in [path for path in jamo_submission.groups if os.path.isdir(path)]:

            new_folder = os.path.join(SCRATCH, os.path.basename(folder))

            dir_util.copy_tree(folder, new_folder, update=1)
            os.chown(new_folder, uid, gid)
            os.chmod(new_folder, 0o750)

            for root, sub_dirs, dir_files in os.walk(new_folder):
                for sub_dir in sub_dirs:
                    os.chown(os.path.join(new_folder, root, sub_dir), uid, gid)
                    os.chmod(os.path.join(new_folder, root, sub_dir), 0o750)

                for new_f in dir_files:
                    if not os.path.islink(os.path.join(new_folder, root, new_f)):
                        try:
                            os.chown(os.path.join(new_folder, root, new_f), uid, gid)
                        except OSError:
                            print os.path.join(new_folder, root, new_f)
                            raise OSError
                        os.chmod(os.path.join(new_folder, root, new_f), 0o640)

                        while True:
                            name, ext = os.path.splitext(new_f)
                            link_f = (
                                name
                                + "-"
                                + re.sub(
                                    r"[^0-9]",
                                    "",
                                    str(datetime.datetime.now().replace(microsecond=0)),
                                )
                                + ext
                            )

                            if os.path.basename(link_f).lower() in [
                                os.path.basename(f).lower() for f in files
                            ]:
                                time.sleep(1)
                                continue
                            else:
                                os.symlink(
                                    os.path.join(new_folder, root, new_f),
                                    os.path.join(new_folder, root, link_f),
                                )
                                break

                        new_f = link_f

                        files.add(
                            os.path.join(
                                os.path.realpath(new_folder),
                                os.path.realpath(root),
                                new_f,
                            )
                        )

    return files


def make_jamo_json(jamo_submissions_df, args):
    """
    Takes each row in jamo_submissions_df, finds all file locations in the submission,
    and creates a json file formatted to submit those files to genome portal.
    """

    # Loop over all submissions not yet created in jamo_submissions_df
    for idx, jamo_submission in jamo_submissions_df[
        jamo_submissions_df.created.isnull()
        & ~pd.isnull(jamo_submissions_df).all(axis=1)
    ].iterrows():
        assert jamo_submission.type in [
            "raw_data",
            "sop",
            "readme",
            "pactolus",
            "spectral_hits",
            "untargeted_analysis",
            "targeted_analysis",
        ]
        assert ~pd.isnull(jamo_submission.id)

        # Print info on submission and its files
        print jamo_submission.id, jamo_submission.type
        files = get_files_to_submit(jamo_submission, args)
        print str(len(files)) + " files"

        # Creates json entry according to the submission's type
        if jamo_submission.type == "raw_data":

            json_submission = {
                "metadata": {"sequencing_project_id": int(jamo_submission.id)}
            }
            json_submission["outputs"] = [
                {
                    "file": os.path.realpath(f),
                    "label": jamo_submission.type + "_mzML"
                    if "mzml" in os.path.splitext(f)[1].lower()
                    else jamo_submission.type + "_h5",
                    "metadata": {
                        "portal_display_location": os.path.normpath(
                            os.path.join(jamo_submission.display_location, "mzml")
                        )
                        if "mzml" in os.path.splitext(f)[1].lower()
                        else os.path.normpath(
                            os.path.join(jamo_submission.display_location, "hdf5")
                        )
                    },
                }
                for f in files
            ]

        elif jamo_submission.type in ["sop", "readme"]:

            json_submission = {
                "metadata": {"sequencing_project_id": int(jamo_submission.id)}
            }
            json_submission["outputs"] = [
                {"file": os.path.realpath(f), "label": jamo_submission.type}
                for f in files
            ]

        elif jamo_submission.type in ["pactolus", "spectral_hits"]:

            json_submission = {
                "metadata": {"analysis_project_id": int(jamo_submission.id)}
            }
            json_submission["outputs"] = [
                {
                    "file": os.path.realpath(f),
                    "label": jamo_submission.type,
                    "metadata": {
                        "portal_display_location": jamo_submission.display_location
                    },
                }
                for f in files
            ]

        elif jamo_submission.type in ["untargeted_analysis", "targeted_analysis"]:

            json_submission = {
                "metadata": {"analysis_project_id": int(jamo_submission.id)}
            }
            json_submission["outputs"] = [
                {
                    "file": f,
                    "label": jamo_submission.type,
                    "file_format": os.path.splitext(f)[1][1:],
                    "metadata": {
                        "portal_display_location": os.path.join(
                            jamo_submission.display_location,
                            *os.path.relpath(f, SCRATCH).split(os.path.sep)[1:-1]
                        )
                    },
                }
                for f in files
            ]

        # Create or update json file with json entry created
        with os.fdopen(
            os.open(
                "./" + str(jamo_submission.id) + "_" + jamo_submission.type + ".json",
                os.O_RDWR | os.O_CREAT,
            ),
            "w+",
        ) as json_file:
            try:
                json_file_data = json.load(json_file)
                json_submission["outputs"][:0] = json_file_data["outputs"]
                json_submission["outputs"] = {
                    js["file"]: js for js in json_submission["outputs"]
                }.values()
            except ValueError:
                pass

            json_file.seek(0)
            json_file.write(
                json.dumps(json_submission, indent=4, separators=(",", ":"))
            )
            json_file.truncate()

        # Update jamo_submissions_df
        jamo_submissions_df.ix[idx, "created"] = datetime.datetime.now()

    # Save updated jamo_submissions_df
    jamo_submissions_df.fillna("").to_csv(
        args.jamo_submissions,
        index=False,
        columns=[
            "created",
            "submitted",
            "id",
            "type",
            "metadata",
            "display_location",
            "prepend",
            "groups",
            "is_group",
            "file_filters",
        ],
    )


# def zip_analyses(jamo_submissions_df, args):
#
#     zip_dir = os.path.realpath(args.zip)
#
#     assert os.path.isdir(zip_dir)
#
#     for idx, jamo_submission in jamo_submissions_df[(jamo_submissions_df.type == 'targeted_analysis') | (jamo_submissions_df.type == 'untargeted_analysis')].iterrows():
#
#         # shutil.make_archive(os.path.join(zip_dir, os.path.basename(jamo_submission.groups[0])),
#         #                     'zip',
#         #                     jamo_submission.groups[0])
#
#         with zipfile.ZipFile(os.path.join(zip_dir, os.path.basename(jamo_submission.groups[0]) + '.zip'), 'w', zipfile.ZIP_DEFLATED, allowZip64 = True) as zip_file:
#             for root, dirs, files in os.walk(jamo_submission.groups[0]):
#                 for f in files:
#                     zip_file.write(os.path.join(root, f), os.path.join(os.path.relpath(root, jamo_submission.groups[0]), f))
#
#         jamo_submissions_df.ix[idx, 'groups'].insert(0, os.path.join(zip_dir, os.path.basename(jamo_submission.groups[0]) + '.zip'))
#         print os.path.basename(jamo_submission.groups[0])
#
#     jamo_submissions_df.fillna('').to_csv(args.jamo_submissions, index=False,
#                                     columns=['created', 'submitted',
#                                              'project_id', 'analysis_id',
#                                              'type', 'template', 'metadata', 'display_location',
#                                              'groups', 'is_group', 'file_filters'])

# def test_jamo_json(jamo_submissions_df, args):
#     for idx, jamo_submission in jamo_submissions_df[jamo_submissions_df.submitted.isnull() & ~pd.isnull(jamo_submissions_df).all(axis=1)].iterrows():
#         assert not all(jamo_submission[['project_id', 'analysis_id']].notnull())
#
#         file = ''
#
#         if pd.notnull(jamo_submission.project_id):
#             file = '_'.join(jamo_submission[['project_id', 'type']].astype(str).tolist()) + '.json'
#         if pd.notnull(jamo_submission.analysis_id):
#             file = '_'.join(jamo_submission[['analysis_id', 'type']].astype(str).tolist()) + '.json'
#
#         if os.path.isfile(os.path.join(os.path.split(__file__)[0], file)):
#
#
# def submit_jamo_json(jamo_submissions_df, args):
#     pass
#
# def update_jamo_json(jamo_submissions_df, args):
#     pass


def arg_parser():
    parser = argparse.ArgumentParser()
    submission_type = parser.add_mutually_exclusive_group(required=True)

    submission_type.add_argument("-m", "--make", action="store_true")
    # submission_type.add_argument('-t', '--test', action='store_true')

    # submission_type.add_argument('-s', '--submit', action='store_true')
    # submission_type.add_argument('-u', '--update', action='store_true')

    parser.add_argument("jamo_submissions", type=str)

    return parser


def main():
    parser = arg_parser()
    args, cmd_kwargs = parser.parse_known_args()

    jamo_submissions_df = pd.read_csv(
        args.jamo_submissions,
        dtype={
            "created": str,
            "submitted": str,
            "id": str,
            "type": str,
            "metadata": bool,
            "display_location": str,
            "prepend": str,
            "groups": str,
            "is_group": bool,
            "file_filters": str,
        },
        parse_dates=["created", "submitted"],
    )
    jamo_submissions_df.file_filters = jamo_submissions_df.file_filters.apply(
        lambda f: literal_eval(f) if not pd.isnull(f) else f
    )
    jamo_submissions_df.groups = jamo_submissions_df.groups.apply(
        lambda g: literal_eval(g) if not pd.isnull(g) else g
    )

    if args.make:
        make_jamo_json(jamo_submissions_df, args)

    # if args.test:
    #     test_jamo_json(jamo_submissions_df, args)

    # if args.submit:
    #     submit_jamo_json(jamo_submissions_df, args)
    #
    # if args.update:
    #     update_jamo_json(jamo_submissions_df, args)


if __name__ == "__main__":
    main()
