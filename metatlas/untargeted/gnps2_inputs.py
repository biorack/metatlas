"""Helpers for selecting GNPS2 workflow inputs."""

from pathlib import Path


def everything_bagel_input_spectra(sftp, local_directory, raw_data, username, polarity):
    """Select local FPS and polarity-specific mzML files that are present at GNPS2."""
    polarities = {"positive": "POS", "negative": "NEG"}
    if polarity not in polarities:
        raise ValueError(f"Unsupported polarity: {polarity}")

    prefix = f"USERUPLOAD/{username}/"
    if not raw_data.startswith(prefix) or ";" in raw_data:
        raise ValueError(f"Invalid GNPS2 raw data path: {raw_data}")
    local_path = Path(local_directory)
    if not local_path.is_dir():
        raise ValueError(f"Local raw data directory not found: {local_directory}")
    selected = set()
    for path in local_path.rglob("*.mzML"):
        if not path.is_file():
            continue
        filename = path.name
        fields = filename.split("_")
        if len(fields) <= 9 or fields[9] not in {"POS", "NEG", "FPS", "EI"}:
            raise ValueError(f"Cannot determine polarity of GNPS2 file: {filename}")
        if ";" in filename:
            raise ValueError(f"GNPS2 filename contains a semicolon: {filename}")
        if fields[9] in {polarities[polarity], "FPS"}:
            if filename in selected:
                raise ValueError(f"Duplicate mzML filename in local raw data: {filename}")
            selected.add(filename)

    if not selected:
        raise ValueError(f"No {polarity} or FPS mzML files found in {local_directory}")
    remote_directory = "/" + raw_data[len(prefix):]
    missing = selected.difference(sftp.listdir(remote_directory))
    if missing:
        raise ValueError(f"Selected mzML files missing at GNPS2: {', '.join(sorted(missing))}")
    return ";".join(f"{raw_data}/{filename}" for filename in sorted(selected))
