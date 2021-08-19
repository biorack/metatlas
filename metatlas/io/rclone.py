""" copy files to Google Drive using rclone """

import configparser
import logging
import subprocess


logger = logging.getLogger(__name__)


class RClone:
    """Access to Google Drive"""

    def __init__(self, rclone_path):
        self.rclone_path = rclone_path

    def config_file(self):
        """Returns path to config file or None"""
        try:
            result = subprocess.check_output([self.rclone_path, "config", "file"], text=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None
        return result.split("\n")[1]

    def get_name_for_id(self, identifier):
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

    def copy_to_drive(self, source, drive, dest_path=None):
        """
        Inputs:
            source: file or directory to copy to drive
            drive: name in the RClone configuration for a location in Google Drive
            dest_path: location under drive to copy to, will create folders if needed
        """
        dest = f"{drive}:" if dest_path is None else f"{drive}:{dest_path}"
        try:
            subprocess.check_output([self.rclone_path, "copy", source, dest], text=True)
        except subprocess.CalledProcessError as err:
            logger.exception(err)
            raise err
        except FileNotFoundError:
            logger.info("rclone not found. Skipping transfer to Google Drive")
