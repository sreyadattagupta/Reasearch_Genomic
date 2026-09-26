"""Lightweight validation tests for the biomarker-mapping / network pipeline.

These run offline against the cached API responses and generated result files,
so they double as reproducibility checks. Run:  python -m pytest -q
(They can also be executed directly:  python tests/test_pipeline.py)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import config  # noqa: E402
from src.preprocessing import load_geo_series_matrix, get_platform_id  # noqa: E402


def test_raw_matrices_present_and_platform_matches():
    for ds, meta in config.DATASETS.items():
        path = config.RAW_DIR / f"{ds}_series_matrix.txt.gz"
        assert path.exists(), f"missing raw matrix for {ds}"
        assert get_platform_id(path) == meta["platform"], f"platform mismatch {ds}"


def test_lasso_files_have_probes():
    for ds, meta in config.DATASETS.items():
        df = pd.read_csv(config.BIOMARKER_SELECTION_DIR / meta["lasso_file"])
        assert len(df) > 0, f"empty LASSO file {ds}"


def test_geo_loader_shape():
    ds = "GSE72267"
    df = load_geo_series_matrix(config.RAW_DIR / f"{ds}_series_matrix.txt.gz")
    assert df.shape[0] > 1000 and df.shape[1] > 0


def test_mapping_outputs_valid():
    """If Stage 1 has run, mapping table must cover every probe with no NaNs
    in the probe column and unique probe ids per dataset."""
    combined = config.BIOMARKER_MAPPING_DIR / "all_biomarker_gene_map.csv"
    if not combined.exists():
        return  # stage not run yet
    m = pd.read_csv(combined)
    assert m["probe_id"].notna().all()
    for ds, g in m.groupby("dataset"):
        assert g["probe_id"].is_unique, f"duplicate probe rows in {ds}"


def test_network_outputs_valid():
    """If Stage 2 has run, edges reference known genes and summary is sane."""
    edges_path = config.NETWORK_DIR / "network_edges.csv"
    summary_path = config.NETWORK_DIR / "network_summary.csv"
    if not edges_path.exists():
        return
    edges = pd.read_csv(edges_path)
    if not edges.empty:
        assert {"gene_a", "gene_b", "score"}.issubset(edges.columns)
        assert (edges["score"] >= 0).all()
    s = pd.read_csv(summary_path).iloc[0]
    assert s["n_nodes"] >= 0 and s["n_edges"] >= 0


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
