import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from metatlas.untargeted.gnps2_inputs import everything_bagel_input_spectra


class EverythingBagelInputsTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.local_data = Path(self.directory.name)
        self.raw_data = "USERUPLOAD/bpbowen/raw_data/egsb/project"
        self.prefix = "20250211_EB_JK_107419-001_project_20250103_EXP_HILICZ_SERIAL"
        self.sftp = Mock()

    def add_file(self, suffix, subdirectory=None):
        directory = self.local_data / subdirectory if subdirectory else self.local_data
        directory.mkdir(exist_ok=True)
        filename = f"{self.prefix}_{suffix}.mzML"
        (directory / filename).touch()
        return filename

    def select(self, polarity):
        return everything_bagel_input_spectra(
            self.sftp, self.local_data, self.raw_data, "bpbowen", polarity
        )

    def test_selects_target_polarity_and_fps_for_each_mode(self):
        pos = self.add_file("POS_MS2_0_QC")
        neg = self.add_file("NEG_MS2_0_sample")
        fps = self.add_file("FPS_MS1_0_InjBl", "nested")
        ei = self.add_file("EI_MS1_0_sample")
        self.sftp.listdir.return_value = [pos, neg, fps, ei, "stale_NEG_file.mzML"]
        for polarity, expected in (("positive", pos), ("negative", neg)):
            with self.subTest(polarity=polarity):
                self.assertEqual(
                    self.select(polarity).split(";"),
                    [f"{self.raw_data}/{name}" for name in sorted((fps, expected))],
                )
        self.sftp.listdir.assert_called_with("/raw_data/egsb/project")

    def test_fps_only_and_empty_selection(self):
        fps = self.add_file("FPS_MS1_0_sample")
        self.sftp.listdir.return_value = [fps]
        self.assertEqual(self.select("negative"), f"{self.raw_data}/{fps}")
        (self.local_data / fps).unlink()
        self.add_file("POS_MS1_0_sample")
        with self.assertRaisesRegex(ValueError, "No negative or FPS"):
            self.select("negative")

    def test_missing_uploaded_file_and_duplicate_basename(self):
        fps = self.add_file("FPS_MS1_0_sample")
        self.sftp.listdir.return_value = []
        with self.assertRaisesRegex(ValueError, "missing at GNPS2"):
            self.select("negative")
        self.sftp.listdir.return_value = [fps]
        self.add_file("FPS_MS1_0_sample", "nested")
        with self.assertRaisesRegex(ValueError, "Duplicate mzML filename"):
            self.select("negative")

    def test_rejects_unclassifiable_files_and_invalid_input(self):
        (self.local_data / "unrecognizable.mzML").touch()
        with self.assertRaisesRegex(ValueError, "Cannot determine polarity"):
            self.select("negative")
        with self.assertRaisesRegex(ValueError, "Unsupported polarity"):
            self.select("unknown")
        with self.assertRaisesRegex(ValueError, "Invalid GNPS2 raw data path"):
            everything_bagel_input_spectra(
                self.sftp, self.local_data, "USERUPLOAD/other/project", "bpbowen", "positive"
            )

    def test_remote_listing_failure_is_not_hidden(self):
        self.add_file("NEG_MS1_0_sample")
        self.sftp.listdir.side_effect = OSError("directory unavailable")
        with self.assertRaisesRegex(OSError, "directory unavailable"):
            self.select("negative")


if __name__ == "__main__":
    unittest.main()
