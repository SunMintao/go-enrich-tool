# GO Enrichment and Clustering Tool

**GO Enrichment and Clustering Tool** is an integrated command-line tool for Gene Ontology (GO) enrichment analysis and clustering of enrichment results. It performs hypergeometric tests for enrichment and uses gene overlap-based similarity (Jaccard / Ochiai) to reduce dimensionality (UMAP / MDS / t‑SNE) and apply unsupervised clustering (Mean‑Shift / HDBSCAN) to the enriched GO terms. Multiple visualization charts and network files are generated to help researchers explore functional modules from a global perspective.

## Key Features

- **GO annotation file format conversion**  
  
  - Automatically identifies the gene ID column in Excel files and extracts all GO IDs.  
  - Uses the GO basic ontology file (go-basic.obo) to retrieve standard term names for GO IDs, outputting a two‑column format (`ID`, `Annotation`) compatible with this tool.  
  - Supports arbitrary column layouts and various delimiters (comma, pipe, semicolon, etc.) – no manual reformatting required.

- **GO enrichment analysis**  
  
  - Computes hypergeometric test \(p\)-values and Benjamini‑Hochberg correction from a local GO annotation Excel file.  
  - Allows specification of ontology sub‑ontologies (BP / CC / MF / ALL), automatically downloads and parses the GO basic ontology file (go‑basic.obo) to obtain namespace information.  
  - Output includes detailed information such as the list of associated genes (`All_index`).

- **Dynamic filtering and redundancy merging**  
  
  - Automatically filters significant terms based on adjusted \(p\)-value, gene count, and median of negative log‑\(p\).  
  - Merges redundant GO terms that share identical gene sets, keeping the term with the smallest \(p\)-value as the representative.

- **Similarity matrix**  
  
  - Computes Jaccard or Ochiai coefficients to construct a similarity/distance matrix based on actual gene overlap between terms.  
  - Supports multi‑processing parallel computation (automatically enabled when the number of terms > 100).

- **Dimensionality reduction and clustering**  
  
  - Dimensionality reduction methods: UMAP (default), MDS, t‑SNE.  
  - Clustering algorithms: Mean‑Shift (adaptive bandwidth), HDBSCAN.  
  - Silhouette coefficient evaluation for clustering quality, automatic identification of boundary points.

- **Rich visualization outputs**  
  
  - Cluster scatter plot, similarity heatmap, bubble chart, static/interactive network (Pyvis).  
  - Cluster‑specific word clouds (WordCloud).  
  - Export of Cytoscape‑compatible node and edge files.

## Installation

### Install from source

```bash
# Clone the repository
git clone https://github.com/SunMintao/go-enrich-tool.git
cd go-enrich-tool

# It is recommended to use a virtual environment
python -m venv venv
source venv/bin/activate   # Linux/macOS
venv\Scripts\activate      # Windows

# Install the package
pip install -e .
```

### Dependencies

Core dependencies are automatically installed (see `pyproject.toml` for details):

- `pandas`, `numpy`, `scipy`, `statsmodels`, `requests`
- `matplotlib`, `seaborn`, `wordcloud`, `pyvis`
- `scikit-learn`, `umap-learn`, `hdbscan`
- `openpyxl`

## Quick Start

### 1. Prepare input files

- **GO annotation file** (Excel, e.g., `go_annotation.xlsx`):  
  Must contain at least two columns: `ID` (gene identifier) and `Annotation` (GO annotations in the format `GO:xxxxxx,description`; multiple entries separated by semicolons).  
  If your annotation file uses a different format, you can use the `convert` subcommand to standardise it (see below).

- **Target gene list file** (Excel, e.g., `id.xlsx`):  
  The first column should contain gene IDs.

### 2. Convert any annotation file to the standard format (optional)

```bash
go-enrich convert any_format.xlsx standard.xlsx --obo go-basic.obo
```

The generated `standard.xlsx` can then be directly used for enrichment analysis.

### 3. Run enrichment analysis

```bash
go-enrich enrich \
    --go-anno data/go_annotation.xlsx \
    --gene-list data/id.xlsx \
    --output results/enrichment.xlsx \
    --ontology ALL
```

### 4. Cluster the enrichment results

```bash
go-enrich cluster \
    --input results/enrichment.xlsx \
    --output-prefix results/clustering_output
```

After execution, multiple files will be generated in the `results/` directory (see the "Output Files" section below).

## Command‑line Arguments

### `enrich` subcommand

| Argument      | Required | Description                                                                  |
| ------------- | -------- | ---------------------------------------------------------------------------- |
| `--go-anno`   | Yes      | Path to the GO annotation Excel file (standard format).                      |
| `--gene-list` | Yes      | Path to the target gene list Excel file.                                     |
| `--output`    | No       | Output Excel file path (default: `go_enrichment_result.xlsx`).               |
| `--ontology`  | No       | Ontology sub‑ontology: `BP`, `CC`, `MF`, `ALL` (default: `BP`).              |
| `--use-obo`   | No       | Whether to download and use the obo file for namespace (enabled by default). |
| `--obo-path`  | No       | Path to the obo file (default: `go-basic.obo`).                              |

### `cluster` subcommand

| Argument                     | Description                                                                                                                                    |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`                    | Path to the enrichment result Excel file (must contain columns: `GO`, `Description`, `P-value`, `Adjusted P-value`, `Count_all`, `All_index`). |
| `--output-prefix`            | Prefix for output files (default: `go_clustering`).                                                                                            |
| **Filtering**                |                                                                                                                                                |
| `--p-cutoff`                 | Adjusted \(p\)-value cutoff (default: 0.05).                                                                                                   |
| `--min-genes`                | Minimum gene count (default: 3).                                                                                                               |
| `--logp-median-ratio`        | Whether to filter by median of negative log‑\(p\) (enabled by default).                                                                        |
| `--keep-top-n`               | Number of top terms to keep when falling back automatically (default: 15).                                                                     |
| `--no-merge-identical`       | Do not merge terms with identical gene sets (merging is enabled by default).                                                                   |
| **Similarity**               |                                                                                                                                                |
| `--similarity`               | Similarity metric: `jaccard` or `ochiai` (default: `jaccard`).                                                                                 |
| `--parallel`                 | Enable parallel computation (enabled by default).                                                                                              |
| **Dimensionality reduction** |                                                                                                                                                |
| `--embed-method`             | Method: `umap`, `mds`, or `tsne` (default: `umap`).                                                                                            |
| `--random-state`             | Random seed (default: 42).                                                                                                                     |
| **Clustering**               |                                                                                                                                                |
| `--cluster-method`           | Algorithm: `meanshift` or `hdbscan` (default: `meanshift`).                                                                                    |
| `--k-neighbors`              | Number of neighbours for Mean‑Shift bandwidth estimation (default: 8).                                                                         |
| `--min-bandwidth`            | Minimum bandwidth for Mean‑Shift (default: 0.3).                                                                                               |
| `--hdbscan-min-cluster-size` | Minimum cluster size for HDBSCAN (default: 3).                                                                                                 |
| `--auto-fallback`            | Automatically adjust parameters when clustering fails (enabled by default).                                                                    |
| **Visualisation & export**   |                                                                                                                                                |
| `--similarity-threshold`     | Threshold for drawing edges in the network (default: 0.3).                                                                                     |
| `--no-plots`                 | Do not generate plots (only output the Excel table).                                                                                           |
| `--wordcloud-dir`            | Directory for word clouds (default: `wordclouds`).                                                                                             |
| `--cytoscape`                | Export Cytoscape‑compatible files (enabled by default).                                                                                        |

### `convert` subcommand

| Argument | Required | Description                                                                                                                    |
| -------- | -------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `input`  | Yes      | Input file path (supports arbitrary Excel formats).                                                                            |
| `output` | Yes      | Output file path (standard‑format `.xlsx`).                                                                                    |
| `--obo`  | No       | Path to the GO basic ontology file (default: `go-basic.obo`; `.txt` suffix is also allowed if the content matches obo format). |

## Output Files

After running the `cluster` subcommand, the following files are generated in the output prefix directory:

| File type            | Example filename               | Description                                                                                                                                                             |
| -------------------- | ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Excel results        | `{prefix}_results.xlsx`        | Contains 4 sheets: `Clustering_Results` (detailed term results), `Cluster_Summary` (cluster summary), `Parameters` (run parameters), `Merged_Terms` (merged term info). |
| Scatter plot         | `{prefix}_scatter.png`         | 2D embedding coordinates with cluster colours.                                                                                                                          |
| Similarity heatmap   | `{prefix}_heatmap.png`         | Heatmap of the similarity matrix sorted by cluster.                                                                                                                     |
| Bubble chart         | `{prefix}_bubble.png`          | Terms vs. \(-\log_{10}(p)\), bubble size represents gene count.                                                                                                         |
| Static network       | `{prefix}_network.png`         | Network of terms with similarity above the threshold.                                                                                                                   |
| Interactive network  | `{prefix}_network.html`        | HTML interactive network (Pyvis) supporting drag, zoom, and hover details.                                                                                              |
| Word cloud directory | `{prefix}_wordclouds/`         | One PNG word cloud per cluster.                                                                                                                                         |
| Cytoscape nodes      | `{prefix}_cytoscape_nodes.txt` | Node information (tab‑separated).                                                                                                                                       |
| Cytoscape edges      | `{prefix}_cytoscape_edges.txt` | Edge list (GO_A, GO_B, Similarity).                                                                                                                                     |

## Complete Workflow Example

```bash
# 1. Convert any annotation file (if needed)
go-enrich convert my_annotation.xlsx standard.xlsx --obo go-basic.obo

# 2. Enrichment analysis
go-enrich enrich --go-anno standard.xlsx --gene-list genes.xlsx --output enrich.xlsx --ontology ALL

# 3. Clustering with Ochiai similarity, t‑SNE embedding, and HDBSCAN
go-enrich cluster --input enrich.xlsx --output-prefix my_clustering \
    --similarity ochiai --embed-method tsne --cluster-method hdbscan

# 4. Only output the table, skip all plots
go-enrich cluster --input enrich.xlsx --output-prefix my_clustering_no_plots --no-plots
```

## Frequently Asked Questions

**Q1: Missing modules `umap` or `hdbscan`?**  
A: Ensure all dependencies are installed: `pip install -e .` will install them automatically. Alternatively, manually install with `pip install umap-learn hdbscan`.

**Q2: Why are all terms clustered into one cluster?**  
A: This may occur when the number of terms is too small or similarity is generally high. Try using `--cluster-method hdbscan` or adjust `--min-bandwidth` / `--k-neighbors`.

**Q3: Some GO terms lack `GO_level_1` classification in the enrichment results?**  
A: Confirm that `--use-obo` is enabled and your network is working (the tool will download `go-basic.obo` automatically). If it still fails, manually download the file and specify it with `--obo-path`.

**Q4: Descriptions are missing after running `convert`?**  
A: Ensure you provide a valid `--obo` file that contains the mapping from GO IDs to names. If the input file itself contains descriptions, the converter will preserve them; otherwise it will fill missing ones from the obo file.

## Citation

If you use this tool in your research, please cite:

> [Sun Mintao]. *GO Enrichment and Clustering Tool*. GitHub repository. https://github.com/SunMintao/go-enrich-tool

## License

This project is licensed under the [MIT License](LICENSE). Feel free to use and modify it.
