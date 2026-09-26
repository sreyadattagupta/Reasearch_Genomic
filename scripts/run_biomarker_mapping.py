"""Stage 1: map LASSO biomarker probes to gene symbols for every dataset.

Run from the repo root:  python -m scripts.run_biomarker_mapping
Outputs (results/biomarker_mapping/):
  <GSE>_biomarker_gene_map.csv   per-dataset probe->gene map
  all_biomarker_gene_map.csv     combined map
  unmapped_probes.csv            probes with no gene hit
  ambiguous_probes.csv           probe -> multiple genes
  duplicate_genes.csv            multiple probes -> same gene
  mapping_summary.csv            per-dataset counts
"""
from __future__ import annotations

import pandas as pd

from src import config
from src.biomarker_mapping import build_mapping
from src.utils.io import save_table


def main(use_cache: bool = True) -> None:
    config.ensure_dirs()

    all_maps, all_unmapped, all_ambiguous, all_dups, summary = [], [], [], [], []

    for dataset in config.DATASETS:
        print(f"[mapping] {dataset} ({config.DATASETS[dataset]['platform']}) ...")
        res = build_mapping(dataset, use_cache=use_cache)
        m = res["mapping"]

        save_table(m, config.BIOMARKER_MAPPING_DIR / f"{dataset}_biomarker_gene_map.csv")
        all_maps.append(m)
        all_unmapped.append(res["unmapped"])
        if not res["ambiguous"].empty:
            res["ambiguous"].insert(0, "dataset", dataset)
            all_ambiguous.append(res["ambiguous"])
        if not res["duplicates"].empty:
            d = res["duplicates"].copy()
            d.insert(0, "dataset", dataset)
            all_dups.append(d)

        n = len(m)
        n_mapped = m["gene_symbol"].notna().sum()
        summary.append({
            "dataset": dataset,
            "platform": config.DATASETS[dataset]["platform"],
            "n_probes": n,
            "n_mapped": int(n_mapped),
            "n_unmapped": int(n - n_mapped),
            "n_ambiguous": int(len(res["ambiguous"])),
            "n_duplicate_genes": int(res["duplicates"]["gene_symbol"].nunique())
            if not res["duplicates"].empty else 0,
            "pct_mapped": round(100 * n_mapped / n, 1) if n else 0.0,
        })

    combined = pd.concat(all_maps, ignore_index=True)
    save_table(combined, config.BIOMARKER_MAPPING_DIR / "all_biomarker_gene_map.csv")
    save_table(pd.concat(all_unmapped, ignore_index=True),
               config.BIOMARKER_MAPPING_DIR / "unmapped_probes.csv")
    save_table(pd.concat(all_ambiguous, ignore_index=True) if all_ambiguous
               else pd.DataFrame(columns=["dataset", "probe_id", "gene_symbols"]),
               config.BIOMARKER_MAPPING_DIR / "ambiguous_probes.csv")
    save_table(pd.concat(all_dups, ignore_index=True) if all_dups
               else pd.DataFrame(columns=["dataset", "probe_id", "gene_symbol"]),
               config.BIOMARKER_MAPPING_DIR / "duplicate_genes.csv")
    summary_df = pd.DataFrame(summary)
    save_table(summary_df, config.BIOMARKER_MAPPING_DIR / "mapping_summary.csv")

    print("\n=== mapping summary ===")
    print(summary_df.to_string(index=False))
    print(f"\nUnique mapped genes across all datasets: "
          f"{combined['gene_symbol'].dropna().nunique()}")


if __name__ == "__main__":
    main()
