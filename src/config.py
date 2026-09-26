"""Central configuration: paths and the dataset/platform registry.

All paths are derived relative to the repository root so the project is
portable across machines (no hard-coded absolute paths).
"""
from __future__ import annotations

from pathlib import Path

# --- Repository layout -----------------------------------------------------
# src/config.py -> src -> <repo root>
ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"          # cached API responses

RESULTS_DIR = ROOT / "results"
BIOMARKER_SELECTION_DIR = RESULTS_DIR / "biomarker_selection"   # existing LASSO output
BIOMARKER_MAPPING_DIR = RESULTS_DIR / "biomarker_mapping"       # new probe->gene maps
NETWORK_DIR = RESULTS_DIR / "network"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"
MODELS_DIR = RESULTS_DIR / "models"

# --- Dataset registry ------------------------------------------------------
# Each GEO series with its platform (GPL) and the identifier type used by the
# probe IDs stored in the LASSO Top-32 result files. `mygene_scopes` tells the
# mygene.info query API which id namespaces to search.
DATASETS = {
    "GSE99039": {
        "platform": "GPL570",
        "platform_name": "Affymetrix Human Genome U133 Plus 2.0",
        "id_type": "affy_probe",
        "mygene_scopes": "reporter",
        "lasso_file": "GSE99039_Top32_LASSO.csv",
    },
    "GSE6613": {
        "platform": "GPL96",
        "platform_name": "Affymetrix Human Genome U133A",
        "id_type": "affy_probe",
        "mygene_scopes": "reporter",
        "lasso_file": "GSE6613_Top32_LASSO.csv",
    },
    "GSE72267": {
        "platform": "GPL571",
        "platform_name": "Affymetrix Human Genome U133A 2.0",
        "id_type": "affy_probe",
        "mygene_scopes": "reporter",
        "lasso_file": "GSE72267_Top32_LASSO.csv",
    },
    "GSE18838": {
        "platform": "GPL5175",
        "platform_name": "Affymetrix Human Exon 1.0 ST (transcript cluster IDs)",
        "id_type": "affy_transcript_cluster",
        "mygene_scopes": "reporter,accession.rna,ensembl.transcript",
        "lasso_file": "GSE18838_Top32_LASSO.csv",
    },
    "GSE57475": {
        "platform": "GPL6947",
        "platform_name": "Illumina HumanHT-12 V3.0",
        "id_type": "illumina_probe",
        "mygene_scopes": "reporter",
        # Illumina ILMN_ probe IDs are not indexed by mygene.info; fall back to
        # the GEO platform (GPL) annotation table, which carries a Symbol column.
        "geo_fallback": True,
        "lasso_file": "GSE57475_Top32_LASSO.csv",
    },
}

SPECIES = "human"  # taxid 9606 for all datasets

# External service endpoints
MYGENE_URL = "https://mygene.info/v3/query"
GEO_PLATFORM_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc={gpl}&targ=self&form=text&view=data"
)
STRING_API = "https://string-db.org/api"
STRING_SPECIES = 9606  # Homo sapiens


def ensure_dirs() -> None:
    """Create all output directories if missing (safe to call repeatedly)."""
    for d in (
        EXTERNAL_DIR,
        BIOMARKER_MAPPING_DIR,
        NETWORK_DIR,
        FIGURES_DIR,
        TABLES_DIR,
    ):
        d.mkdir(parents=True, exist_ok=True)
