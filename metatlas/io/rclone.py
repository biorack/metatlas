""" copy files to Google Drive using rclone """

import configparser
import json
import logging
import subprocess

from pathlib import Path
from subprocess import Popen, PIPE
from typing import List, Optional, Tuple, Union

from tqdm.notebook import tqdm

from metatlas.tools.notebook import in_papermill

logger = logging.getLogger(__name__)


class RClone:
    """Access to Google Drive"""

    def __init__(self, rclone_path: Union[Path, str]) -> None:
        self.rclone_path = Path(rclone_path)

    def config_file(self) -> Optional[str]:
        """Returns path to config file or None"""
        try:
            result = subprocess.check_output([self.rclone_path, "config", "file"], text=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None
        return result.split("\n")[1]

    def get_name_for_id(self, identifier: str) -> Optional[str]:
        """
        Inputs:
            identifer: unique folder identifier from Google Drive URL
        if identifier is in the config file, then return the name assigned to the identifier
        otherwise return None
        """
        ini_file = self.config_file()
        if ini_file is None:
            return None
        config = configparser.ConfigParser()
        config.read(ini_file)
        for name in config.sections():
            props = config[name]
            if "type" in props and props["type"] == "drive":
                if "root_folder_id" in props and props["root_folder_id"] == identifier:
                    return name
        return None

    def copy_to_drive(
        self, source: Path, drive: str, dest_path: Optional[Path] = None, progress: bool = False
    ) -> None:
        """
        Inputs:
            source: file or directory to copy to drive
            drive: name in the RClone configuration for a location in Google Drive
            dest_path: location under drive to copy to, will create folders if needed
            progress: display a tqdm progress bar
        """
        dest = Path(f"{drive}:" if dest_path is None else f"{drive}:{dest_path}")
        cmd = [self.rclone_path, "copy", source, dest]
        try:
            if progress:
                cmd.append("--progress")
                with tqdm(
                    total=100, desc="Percentage of total number of files", unit="%", disable=in_papermill()
                ) as pbar:
                    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                        for line in proc.stdout or []:
                            line = line.strip()
                            if line.startswith("Transferred:") and line.endswith("%"):
                                percent = float(line.split(",")[1].split("%")[0])
                                pbar.n = percent
                                pbar.refresh()
            else:
                subprocess.check_output(cmd, text=True)
        except subprocess.CalledProcessError as err:
            logger.exception(err)
            raise err
        except FileNotFoundError:
            logger.info("rclone not found. Skipping transfer to Google Drive")

    def get_id_for_path(self, path_string: str) -> str:
        """
        Inputs:
            path_string: a string containing drive_name a colon and one or more folders like:
                         'my_drive:folder1/folder2'
        returns an ID string which can be used in a Google Drive URL
        """
        drive, folders = parse_path(path_string)
        assert isinstance(folders, list)
        assert isinstance(folders[:-1], list)
        all_but_last = f"{drive}:{'/'.join(folders[:-1])}"
        command_list = [self.rclone_path, "lsjson", "--dirs-only", all_but_last]
        try:
            result = subprocess.check_output(command_list, text=True)
        except subprocess.CalledProcessError as err:
            logger.exception(err)
            raise err
        returned_folders = json.loads(result)
        for folder in returned_folders:
            if folder["Name"] == folders[-1]:
                return folder["ID"]
        raise FileNotFoundError(f"Could not find a file or folder at {path_string}")

    def path_to_url(self, path_string: str) -> str:
        """
        Inputs:
            path_string: a string containing drive_name a colon and one or more folders like:
                         'my_drive:folder1/folder2'
        returns an URL for opening the object at path_string
        """
        drive_id = self.get_id_for_path(path_string)
        return f"https://drive.google.com/drive/folders/{drive_id}"


def parse_path(path_string: str) -> Tuple[str, List[str]]:
    """
    Inputs:
        path_string: a string containing drive_name a colon and one or more folders like:
                     'my_drive:folder1/folder2'
    returns a tuple of the drive_name, folder_list
    """
    drive = path_string.split(":")[0]
    remainder = ":".join(path_string.split(":")[1:])
    return drive, remainder.split("/")
