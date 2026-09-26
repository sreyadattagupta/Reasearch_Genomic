# Methodology

## 1. Datasets

Five public Parkinson's disease (PD) blood/brain gene-expression series from NCBI GEO:

| Dataset  | Platform (GPL) | Technology                              | Probe ID format             |
|----------|----------------|-----------------------------------------|-----------------------------|
| GSE99039 | GPL570         | Affymetrix HG-U133 Plus 2.0             | `1554241_at`                |
| GSE6613  | GPL96          | Affymetrix HG-U133A                     | `209950_s_at`               |
| GSE72267 | GPL571         | Affymetrix HG-U133A 2.0                 | `220906_at`                 |
| GSE18838 | GPL5175        | Affymetrix Human Exon 1.0 ST            | numeric transcript clusters |
| GSE57475 | GPL6947        | Illumina HumanHT-12 V3.0                | `ILMN_2359029`              |

Raw series matrices live in `data/raw/`. `data/processed/` holds the GSE99039
expression matrix and its sample labels (0 = control, 1 = PD).

## 2. Existing workflow (preserved)

Implemented in the original notebooks (`notebooks/preprocessing/`):

1. Download GEO series matrices (`00_download_geo_matrices.ipynb`).
2. Load series matrix → transpose to samples × probes.
3. Standardize features (`StandardScaler`).
4. **Biomarker selection:** L1-penalized logistic regression (LASSO) →
   rank probes by |coefficient|, keep top 32 per dataset
   (`results/biomarker_selection/GSE*_Top32_LASSO.csv`).
5. Classification models: Keras MLP (`.h5`) and scikit-learn MLP (`.pkl`),
   stored in `results/models/`.

This methodology and its results are unchanged by this extension.

## 3. Biomarker → gene mapping (new)

The LASSO biomarkers are **platform probe IDs**, not gene symbols. Mapping is
performed by `src/biomarker_mapping`:

- **Primary source — mygene.info** (`/v3/query`): batch-resolves reporter/probe
  IDs to gene `symbol`, `name`, `entrezgene`, and `ensembl.gene`
  (species = human). Scopes are set per platform in `src/config.py`.
- **Fallback — GEO platform (GPL) annotation:** Illumina `ILMN_` IDs (GPL6947)
  are not indexed by mygene.info, so unmapped probes are resolved from the GEO
  platform table's `Symbol` column. Only the needed subset is cached.
- **Caching:** every raw API response is stored under `data/external/`, so
  re-runs are reproducible and offline.

### Reporting (not silent dropping)
For every dataset the pipeline records:
- **unmapped** probes (`unmapped_probes.csv`) — no gene hit,
- **ambiguous** probes (`ambiguous_probes.csv`) — one probe → several genes
  (best hit by mygene `_score` is used for the primary map),
- **duplicate** genes (`duplicate_genes.csv`) — several probes → one gene.

A per-dataset `mapping_summary.csv` gives counts and percent mapped.

## 4. Gene-interaction network (new)

`src/network_analysis` builds a functional-association network from the mapped
gene symbols:

- **Source — STRING v12 REST API** (`/tsv/network`), *Homo sapiens*, confidence
  threshold `required_score = 400` (medium). STRING scores integrate
  experimental, database, co-expression, and text-mining evidence.
- Graph built with **NetworkX** (undirected, edge weight = STRING score). All
  mapped genes are added as nodes so isolated biomarkers remain visible.
- **Analysis:** degree, degree/betweenness/closeness/eigenvector centrality,
  hub-gene ranking, and connected components (`node_metrics.csv`,
  `hub_genes.csv`, `network_summary.csv`).
- **Visualization:** spring layout, node size ∝ degree, colour = component,
  hub labels (`results/figures/biomarker_network.png`).

## 5. Limitations & reproducibility notes

- Probe→gene mapping depends on live annotation databases; cached responses in
  `data/external/` pin the exact results used here.
- GSE18838 (Exon 1.0 ST) transcript-cluster IDs map at the transcript-cluster
  level; a few resolve ambiguously.
- The network is sparse (many isolated nodes): the biomarkers come from five
  different platforms/cohorts and are not expected to be densely
  interconnected at STRING medium confidence. Lower `required_score` in
  `scripts/run_network_analysis.py` for a denser graph.
- STRING edges are *associations*, not experimentally proven physical PD
  interactions.
