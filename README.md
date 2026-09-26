# Reasearch_Genomic — Parkinson's Biomarker Discovery, Gene Mapping & Network Analysis

Machine-learning identification of Parkinson's disease (PD) gene-expression
biomarkers from public GEO datasets, extended with **biomarker → gene mapping**
and **gene-interaction network analysis**.

## Objective

1. Select candidate PD biomarkers from microarray gene-expression data
   (LASSO logistic regression) — *existing work*.
2. Map the selected probe-level biomarkers to human gene symbols/annotations.
3. Build and analyse a gene-interaction network of the mapped biomarkers
   (hubs, centrality, connected components) with clear visualizations.

## Datasets

Five NCBI GEO PD series (see `docs/methodology.md` for platform details):
`GSE99039`, `GSE6613`, `GSE72267`, `GSE18838`, `GSE57475`.
Raw series matrices are in `data/raw/`.

## Folder structure

```
Reasearch_Genomic/
├── data/
│   ├── raw/          # GEO *_series_matrix.txt.gz
│   ├── processed/    # GSE99039 expression matrix + labels
│   └── external/     # cached mygene.info / STRING / GPL responses
├── src/
│   ├── config.py               # paths + dataset/platform registry
│   ├── preprocessing/          # GEO series-matrix loader
│   ├── biomarker_mapping/      # probe -> gene mapping (mygene + GEO fallback)
│   ├── network_analysis/       # STRING network + NetworkX metrics + plots
│   └── utils/
├── scripts/
│   ├── run_biomarker_mapping.py
│   └── run_network_analysis.py
├── notebooks/
│   └── preprocessing/          # original exploratory + LASSO notebooks
├── results/
│   ├── biomarker_selection/    # LASSO Top-32 per dataset (existing output)
│   ├── biomarker_mapping/      # probe->gene maps + unmapped/ambiguous/dup reports
│   ├── network/                # edges, node metrics, hub genes, summary
│   ├── figures/                # network visualization
│   └── models/                 # trained .h5 / .pkl classifiers
├── docs/methodology.md
├── tests/test_pipeline.py
├── requirements.txt
└── README.md
```

## Installation

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate     Linux/macOS:  source .venv/bin/activate
pip install -r requirements.txt
```

(`tensorflow` is optional — only needed to load/retrain the Keras `.h5` models.)

## How to run the pipeline

Run from the repository root (module form matters for imports):

```bash
# Stage 1 — map biomarker probes to genes
python -m scripts.run_biomarker_mapping

# Stage 2 — build & analyse the gene-interaction network
python -m scripts.run_network_analysis

# Validation
python tests/test_pipeline.py        # or: python -m pytest -q
```

Both stages cache API responses in `data/external/`, so subsequent runs are
reproducible and work offline.

## Expected outputs

- `results/biomarker_mapping/` — per-dataset `*_biomarker_gene_map.csv`,
  `all_biomarker_gene_map.csv`, `mapping_summary.csv`, and
  `unmapped_probes.csv` / `ambiguous_probes.csv` / `duplicate_genes.csv`.
- `results/network/` — `network_edges.csv`, `node_metrics.csv`,
  `hub_genes.csv`, `network_summary.csv`.
- `results/figures/biomarker_network.png` — the network visualization.

### Latest run summary

| Dataset  | Platform | Probes | Mapped | % |
|----------|----------|--------|--------|-----|
| GSE99039 | GPL570   | 32 | 27 | 84.4 |
| GSE6613  | GPL96    | 32 | 30 | 93.8 |
| GSE72267 | GPL571   | 32 | 27 | 84.4 |
| GSE18838 | GPL5175  | 32 | 28 | 87.5 |
| GSE57475 | GPL6947  | 32 | 29 | 90.6 |

139 unique mapped genes → STRING network: **139 nodes, 44 edges**, largest
connected component 35 genes. Top hubs include **BCL2, BDNF, CCT2, CD74,
BCLAF1**; known PD-associated genes **PACRG** and **SYT11** appear in the main
component.

## Data / API / database sources

- **Gene expression:** NCBI GEO series matrices.
- **Probe → gene mapping:** [mygene.info](https://mygene.info) `/v3/query`;
  fallback to NCBI GEO **GPL platform annotation** tables (Illumina GPL6947).
- **Gene interactions:** [STRING](https://string-db.org) v12 REST API
  (*Homo sapiens*, confidence ≥ 400).

## Methodology & limitations

See `docs/methodology.md`. Key notes: mappings depend on live databases (pinned
via `data/external/` caches); STRING edges are functional associations, not
proven physical PD interactions; the network is sparse because biomarkers span
five platforms/cohorts.

## Reproducibility

- No hard-coded absolute paths — everything derives from the repo root via
  `src/config.py`.
- Deterministic layout seed for the network figure.
- Cached external responses committed under `data/external/`.
