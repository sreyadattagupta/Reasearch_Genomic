"""Stage 2: build and analyse the biomarker gene-interaction network.

Run from the repo root:  python -m scripts.run_network_analysis
Depends on Stage 1 output (results/biomarker_mapping/all_biomarker_gene_map.csv).

Outputs:
  results/network/network_edges.csv        STRING edges (gene_a, gene_b, score)
  results/network/node_metrics.csv         degree/centrality per gene
  results/network/hub_genes.csv            top hub genes by degree
  results/network/network_summary.csv      global network statistics
  results/figures/biomarker_network.png    network visualization
"""
from __future__ import annotations

import pandas as pd

from src import config
from src.network_analysis import (analyse, build_graph, draw, fetch_string_edges,
                                   summarise)
from src.utils.io import save_table

REQUIRED_SCORE = 400  # STRING medium confidence
TOP_HUBS = 20


def main(use_cache: bool = True) -> None:
    config.ensure_dirs()

    map_path = config.BIOMARKER_MAPPING_DIR / "all_biomarker_gene_map.csv"
    if not map_path.exists():
        raise SystemExit("Run scripts.run_biomarker_mapping first (mapping file missing).")

    mapping = pd.read_csv(map_path)
    genes = sorted(mapping["gene_symbol"].dropna().unique().tolist())
    print(f"[network] {len(genes)} unique mapped genes -> querying STRING ...")

    edges = fetch_string_edges(genes, required_score=REQUIRED_SCORE,
                               tag="combined", use_cache=use_cache)
    save_table(edges if not edges.empty
               else pd.DataFrame(columns=["gene_a", "gene_b", "score"]),
               config.NETWORK_DIR / "network_edges.csv")

    G = build_graph(edges, all_genes=genes)
    metrics = analyse(G)
    save_table(metrics, config.NETWORK_DIR / "node_metrics.csv")
    save_table(metrics.head(TOP_HUBS), config.NETWORK_DIR / "hub_genes.csv")

    stats = summarise(G)
    save_table(pd.DataFrame([stats]), config.NETWORK_DIR / "network_summary.csv")

    draw(G, metrics, config.FIGURES_DIR / "biomarker_network.png",
         title="Parkinson's biomarker gene-interaction network (STRING, score>=400)")

    print("\n=== network summary ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")
    print("\n=== top hub genes ===")
    if not metrics.empty:
        print(metrics.head(TOP_HUBS)[["gene", "degree", "betweenness",
                                       "component"]].to_string(index=False))


if __name__ == "__main__":
    main()
