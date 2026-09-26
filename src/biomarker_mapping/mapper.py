"""Biomarker (probe) -> gene mapping via the mygene.info REST API.

Input  : LASSO Top-32 result files (one per GEO dataset). The ``Gene`` column
         holds *platform probe identifiers*, not gene symbols.
Output : a tidy mapping table plus an explicit report of probes that are
         unmapped, ambiguous (probe -> many genes) or duplicated
         (many probes -> one gene).

Source & methodology
--------------------
Identifiers are resolved with the mygene.info query API
(https://mygene.info/v3/query), which maps microarray reporter/probe IDs to
NCBI/Ensembl gene records. The raw JSON response is cached under
``data/external`` so re-runs are reproducible and work offline.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

from src import config


def _cache_path(dataset: str) -> Path:
    return config.EXTERNAL_DIR / f"mygene_{dataset}.json"


def _query_mygene(probes: list[str], scopes: str) -> list[dict]:
    """POST a batch query to mygene.info and return the raw hit list."""
    resp = requests.post(
        config.MYGENE_URL,
        data={
            "q": ",".join(probes),
            "scopes": scopes,
            "fields": "symbol,name,entrezgene,ensembl.gene",
            "species": config.SPECIES,
        },
        timeout=60,
    )
    resp.raise_for_status()
    return resp.json()


def fetch_mapping(dataset: str, probes: list[str], use_cache: bool = True) -> list[dict]:
    """Return mygene.info hits for ``probes``, caching the raw response."""
    cache = _cache_path(dataset)
    if use_cache and cache.exists():
        return json.loads(cache.read_text())

    scopes = config.DATASETS[dataset]["mygene_scopes"]
    hits: list[dict] = []
    # Batch to stay well within API limits; Top-32 fits one batch but this
    # keeps the function safe for larger biomarker sets.
    for i in range(0, len(probes), 100):
        batch = probes[i : i + 100]
        hits.extend(_query_mygene(batch, scopes))
        time.sleep(0.34)  # be polite to the public endpoint

    config.EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(hits, indent=2))
    return hits


def geo_platform_map(dataset: str, probes: list[str],
                     use_cache: bool = True) -> dict[str, dict]:
    """Fallback: map probes via the GEO platform (GPL) annotation table.

    Used for platforms (e.g. Illumina GPL6947) whose probe IDs mygene.info does
    not index. Only the subset of rows matching ``probes`` is cached (a small
    CSV) rather than the full multi-MB platform file.
    """
    gpl = config.DATASETS[dataset]["platform"]
    cache = config.EXTERNAL_DIR / f"geo_{gpl}_subset.csv"
    wanted = set(probes)

    if use_cache and cache.exists():
        sub = pd.read_csv(cache, dtype=str)
        if wanted.issubset(set(sub["ID"])):
            return _geo_rows_to_dict(sub, wanted)

    from io import StringIO
    url = config.GEO_PLATFORM_URL.format(gpl=gpl)
    text = requests.get(url, timeout=120).text
    start = None
    lines = text.splitlines(keepends=True)
    for i, ln in enumerate(lines):
        if ln.startswith("!platform_table_begin"):
            start = i + 1
            break
    if start is None:
        return {}
    end = len(lines)
    for i in range(start, len(lines)):
        if lines[i].startswith("!platform_table_end"):
            end = i
            break
    table = pd.read_csv(StringIO("".join(lines[start:end])), sep="\t", dtype=str)
    table["ID"] = table["ID"].astype(str).str.strip()
    sub = table[table["ID"].isin(wanted)].copy()

    config.EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
    keep = [c for c in ["ID", "Symbol", "Entrez_Gene_ID", "RefSeq_ID",
                        "Definition", "ILMN_Gene"] if c in sub.columns]
    sub[keep].to_csv(cache, index=False)
    return _geo_rows_to_dict(sub, wanted)


def _geo_rows_to_dict(sub: pd.DataFrame, wanted: set[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for _, r in sub.iterrows():
        pid = str(r["ID"]).strip()
        if pid not in wanted:
            continue
        sym = str(r.get("Symbol", "") or "").strip()
        if not sym or sym.lower() in {"nan", ""}:
            continue
        out[pid] = {
            "gene_symbol": sym,
            "gene_name": (str(r.get("Definition", "")).strip() or None),
            "entrezgene": (str(r.get("Entrez_Gene_ID", "")).strip() or None),
            "ensembl_gene": None,
        }
    return out


def _read_probes(dataset: str) -> pd.DataFrame:
    """Load a LASSO Top-32 file, normalising column names.

    Some files carry only a ``Gene`` column; others add Coefficient/Importance.
    """
    path = config.BIOMARKER_SELECTION_DIR / config.DATASETS[dataset]["lasso_file"]
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "probe_id"})
    df["probe_id"] = df["probe_id"].astype(str).str.strip()
    return df


def build_mapping(dataset: str, use_cache: bool = True) -> dict[str, pd.DataFrame]:
    """Map one dataset's biomarker probes to genes.

    Returns a dict with three DataFrames:
      - ``mapping``  : one row per probe with the resolved gene (best hit)
      - ``unmapped`` : probes with no gene hit
      - ``ambiguous``: probes that resolve to more than one gene
    """
    probes_df = _read_probes(dataset)
    probes = probes_df["probe_id"].tolist()
    hits = fetch_mapping(dataset, probes, use_cache=use_cache)

    # Group hits by the queried probe id.
    by_probe: dict[str, list[dict]] = {}
    for h in hits:
        by_probe.setdefault(h.get("query"), []).append(h)

    # Optional GEO platform-annotation fallback for probes mygene can't resolve.
    fallback: dict[str, dict] = {}
    if config.DATASETS[dataset].get("geo_fallback"):
        need = [p for p in probes
                if not [h for h in by_probe.get(p, [])
                        if not h.get("notfound") and h.get("symbol")]]
        if need:
            fallback = geo_platform_map(dataset, need, use_cache=use_cache)

    rows = []
    unmapped = []
    ambiguous = []
    for probe in probes:
        ph = by_probe.get(probe, [])
        valid = [h for h in ph if not h.get("notfound") and h.get("symbol")]
        if not valid:
            if probe in fallback:
                fb = fallback[probe]
                rows.append({"probe_id": probe, **fb, "n_gene_hits": 1,
                             "source": "geo_platform"})
                continue
            unmapped.append(probe)
            rows.append({"probe_id": probe, "gene_symbol": None,
                         "gene_name": None, "entrezgene": None,
                         "ensembl_gene": None, "n_gene_hits": 0,
                         "source": None})
            continue

        symbols = sorted({h["symbol"] for h in valid})
        if len(symbols) > 1:
            ambiguous.append({"probe_id": probe,
                              "gene_symbols": "|".join(symbols)})

        # Best hit = highest mygene _score (first is best per API ordering).
        best = max(valid, key=lambda h: h.get("_score", 0))
        ens = best.get("ensembl")
        if isinstance(ens, list):
            ens_gene = "|".join(e.get("gene", "") for e in ens)
        elif isinstance(ens, dict):
            ens_gene = ens.get("gene")
        else:
            ens_gene = None
        rows.append({
            "probe_id": probe,
            "gene_symbol": best.get("symbol"),
            "gene_name": best.get("name"),
            "entrezgene": best.get("entrezgene"),
            "ensembl_gene": ens_gene,
            "n_gene_hits": len(symbols),
            "source": "mygene",
        })

    mapping = probes_df.merge(pd.DataFrame(rows), on="probe_id", how="left")
    mapping.insert(0, "dataset", dataset)

    # Duplicates: multiple probes collapsing to the same gene symbol.
    mapped_only = mapping.dropna(subset=["gene_symbol"])
    dup_mask = mapped_only["gene_symbol"].duplicated(keep=False)
    duplicates = (mapped_only[dup_mask]
                  .sort_values("gene_symbol")[["probe_id", "gene_symbol"]])

    return {
        "mapping": mapping,
        "unmapped": pd.DataFrame({"dataset": dataset, "probe_id": unmapped}),
        "ambiguous": pd.DataFrame(ambiguous),
        "duplicates": duplicates,
    }
