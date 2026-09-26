"""Generate README figures from ACTUAL pipeline results (no synthetic data).

Run from repo root:  python -m scripts.make_figures
Produces, in results/figures/:
  preprocessing.png     GSE99039 class balance (Control vs PD)
  biomarker_analysis.png  Top LASSO biomarker importances (GSE57475)
  gene_mapping.png      mapped vs unmapped probes per dataset
  hub_genes.png         top hub genes by network degree
(The main network figure biomarker_network.png is produced by
 scripts/run_network_analysis.py.)
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src import config

BLUE, ORANGE, GREY = "#2b6cb0", "#dd6b20", "#cbd5e0"


def fig_preprocessing() -> None:
    labels = pd.read_csv(config.PROCESSED_DIR / "GSE99039_labels.csv", header=None)
    counts = labels[0].value_counts().sort_index()
    names = {0: "Control", 1: "Parkinson's"}
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.bar([names.get(i, str(i)) for i in counts.index], counts.values,
           color=[GREY, BLUE])
    for i, v in enumerate(counts.values):
        ax.text(i, v + 0.5, str(v), ha="center", fontweight="bold")
    ax.set_title("GSE99039 sample groups (after preprocessing)")
    ax.set_ylabel("Number of samples")
    ax.set_ylim(0, counts.max() + 6)
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "preprocessing.png", dpi=200)
    plt.close(fig)


def fig_biomarker_analysis(dataset: str = "GSE57475", top: int = 15) -> None:
    df = pd.read_csv(config.BIOMARKER_SELECTION_DIR /
                     config.DATASETS[dataset]["lasso_file"])
    df = df.sort_values("Importance", ascending=True).tail(top)
    colors = [ORANGE if c < 0 else BLUE for c in df["Coefficient"]]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(df["Gene"].astype(str), df["Importance"], color=colors)
    ax.set_title(f"Top {top} LASSO biomarkers by importance ({dataset})")
    ax.set_xlabel("Importance (|LASSO coefficient|)")
    ax.text(0.98, 0.02, "blue = ↑ in PD   orange = ↓ in PD",
            transform=ax.transAxes, ha="right", fontsize=8, color="#555")
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "biomarker_analysis.png", dpi=200)
    plt.close(fig)


def fig_gene_mapping() -> None:
    s = pd.read_csv(config.BIOMARKER_MAPPING_DIR / "mapping_summary.csv")
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(s["dataset"], s["n_mapped"], color=BLUE, label="Mapped to gene")
    ax.bar(s["dataset"], s["n_unmapped"], bottom=s["n_mapped"], color=GREY,
           label="Unmapped")
    for i, r in s.iterrows():
        ax.text(i, r["n_probes"] + 0.4, f"{r['pct_mapped']:.0f}%",
                ha="center", fontsize=9, fontweight="bold")
    ax.set_title("Biomarker probe → gene mapping success per dataset")
    ax.set_ylabel("Number of probes (of 32)")
    ax.set_ylim(0, 38)
    ax.legend()
    plt.xticks(rotation=20)
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "gene_mapping.png", dpi=200)
    plt.close(fig)


def fig_hub_genes(top: int = 15) -> None:
    m = pd.read_csv(config.NETWORK_DIR / "node_metrics.csv")
    m = m[m["degree"] > 0].sort_values("degree", ascending=True).tail(top)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(m["gene"], m["degree"], color=BLUE)
    ax.set_title(f"Top {top} hub genes (most connections in the network)")
    ax.set_xlabel("Degree (number of interacting genes)")
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "hub_genes.png", dpi=200)
    plt.close(fig)


def main() -> None:
    config.ensure_dirs()
    fig_preprocessing()
    fig_biomarker_analysis()
    fig_gene_mapping()
    fig_hub_genes()
    print("Figures written to", config.FIGURES_DIR)


if __name__ == "__main__":
    main()
