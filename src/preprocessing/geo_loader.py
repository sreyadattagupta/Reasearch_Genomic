"""Load GEO series-matrix files.

Extracted and consolidated from the original exploratory notebooks
(sreya/sikta/eshmita/deepsikha) so the loading logic lives in one place.
"""
from __future__ import annotations

import gzip
from io import StringIO
from pathlib import Path

import pandas as pd


def load_geo_series_matrix(path: str | Path) -> pd.DataFrame:
    """Parse a GEO ``*_series_matrix.txt.gz`` file into a DataFrame.

    Returns a frame indexed by probe ID (``ID_REF``) with one column per
    sample (GSM accession). This matches the orientation used throughout the
    original notebooks; transpose for a samples-x-probes matrix.
    """
    path = Path(path)
    # Force UTF-8 (GEO files are UTF-8); avoids Windows cp1252 decode errors.
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    start = end = None
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif line.startswith("!series_matrix_table_end"):
            end = i
            break

    if start is None or end is None:
        raise ValueError(f"Series matrix table markers not found in {path}")

    data = "".join(lines[start:end])
    df = pd.read_csv(StringIO(data), sep="\t", index_col=0)
    return df


def get_platform_id(path: str | Path) -> str | None:
    """Return the GPL platform id declared in a series-matrix file."""
    path = Path(path)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_platform_id"):
                return line.split("\t", 1)[1].strip().strip('"')
    return None
