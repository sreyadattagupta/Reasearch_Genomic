"""Gene-interaction network construction and analysis.

Interactions are retrieved from the STRING database REST API
(https://string-db.org/api) for the mapped biomarker gene symbols. STRING
edges are functional-association scores (0-1000) integrating experimental,
database, co-expression and text-mining evidence. The raw TSV is cached under
``data/external`` for reproducibility.

The network is analysed with NetworkX (degree/centrality, hub genes,
connected components) and drawn with matplotlib.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless / reproducible rendering
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import requests

from src import config


def fetch_string_edges(genes: list[str], required_score: int = 400,
                       tag: str = "combined", use_cache: bool = True) -> pd.DataFrame:
    """Fetch STRING functional-association edges among ``genes``.

    ``required_score`` is STRING's confidence threshold (400 = medium).
    """
    genes = sorted({g for g in genes if g})
    cache = config.EXTERNAL_DIR / f"string_network_{tag}.tsv"
    if use_cache and cache.exists():
        text = cache.read_text()
    else:
        resp = requests.post(
            f"{config.STRING_API}/tsv/network",
            data={
                "identifiers": "\r".join(genes),
                "species": config.STRING_SPECIES,
                "required_score": required_score,
                "caller_identity": "reasearch_genomic_pipeline",
            },
            timeout=120,
        )
        resp.raise_for_status()
        text = resp.text
        config.EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
        cache.write_text(text)

    from io import StringIO
    df = pd.read_csv(StringIO(text), sep="\t")
    if df.empty:
        return df
    # Keep the useful columns; STRING scores come as 0-1 floats here.
    edges = df[["preferredName_A", "preferredName_B", "score"]].copy()
    edges.columns = ["gene_a", "gene_b", "score"]
    return edges


def build_graph(edges: pd.DataFrame, all_genes: list[str] | None = None) -> nx.Graph:
    """Build an undirected weighted graph from STRING edges.

    If ``all_genes`` is given, isolated (unconnected) genes are added as nodes
    so the network reflects the full biomarker set, not just connected ones.
    """
    G = nx.Graph()
    if all_genes:
        G.add_nodes_from(g for g in all_genes if g)
    for _, r in edges.iterrows():
        G.add_edge(r["gene_a"], r["gene_b"], weight=float(r["score"]))
    return G


def analyse(G: nx.Graph) -> pd.DataFrame:
    """Per-node centrality metrics, sorted by degree (hub genes first)."""
    if G.number_of_nodes() == 0:
        return pd.DataFrame()
    deg = dict(G.degree())
    metrics = pd.DataFrame({
        "gene": list(deg.keys()),
        "degree": list(deg.values()),
        "degree_centrality": pd.Series(nx.degree_centrality(G)),
        "betweenness": pd.Series(nx.betweenness_centrality(G, weight="weight")),
        "closeness": pd.Series(nx.closeness_centrality(G)),
    })
    try:
        metrics["eigenvector"] = pd.Series(
            nx.eigenvector_centrality_numpy(G, weight="weight"))
    except Exception:
        metrics["eigenvector"] = float("nan")
    # component membership
    comp_id = {}
    for i, comp in enumerate(nx.connected_components(G)):
        for n in comp:
            comp_id[n] = i
    metrics["component"] = metrics["gene"].map(comp_id)
    return (metrics.sort_values(["degree", "betweenness"], ascending=False)
            .reset_index(drop=True))


def summarise(G: nx.Graph) -> dict:
    components = list(nx.connected_components(G))
    return {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "n_components": len(components),
        "largest_component_size": max((len(c) for c in components), default=0),
        "density": nx.density(G) if G.number_of_nodes() > 1 else 0.0,
        "n_isolated": len(list(nx.isolates(G))),
    }


def draw(G: nx.Graph, metrics: pd.DataFrame, out_path: str | Path,
         title: str = "Biomarker gene-interaction network",
         label_top: int = 15) -> None:
    """Render the network: node size ~ degree, colour ~ connected component."""
    out_path = Path(out_path)
    if G.number_of_nodes() == 0:
        return

    deg = dict(G.degree())
    comp_id = {n: c for n, c in metrics.set_index("gene")["component"].items()}
    node_sizes = [80 + 120 * deg.get(n, 0) for n in G.nodes()]
    node_colors = [comp_id.get(n, 0) for n in G.nodes()]

    hubs = set(metrics.head(label_top)["gene"])
    labels = {n: n for n in G.nodes() if n in hubs}

    pos = nx.spring_layout(G, seed=42, k=0.5)
    fig, ax = plt.subplots(figsize=(14, 11))
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.6, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                           cmap=plt.cm.tab20, alpha=0.9, ax=ax)
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9,
                            font_weight="bold", ax=ax)
    ax.set_title(title, fontsize=15)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
