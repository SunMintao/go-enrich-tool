# GO Enrichment and Clustering Tool

**GO Enrichment and Clustering Tool** is an integrated command-line tool for Gene Ontology (GO) enrichment analysis and clustering of enrichment results. It performs hypergeometric tests for enrichment and uses gene-overlap-based similarity (Jaccard / Ochiai) to reduce dimensionality (UMAP / MDS / t-SNE) and apply unsupervised clustering (Mean-Shift / HDBSCAN) to the enriched GO terms. Multiple visualization charts and network files are generated to help researchers explore functional modules from a global perspective.

## Key Features

- **NCBI GO annotation download**
  
  - Downloads the GO annotation (GAF) file for a species by its Latin name from the NCBI RefSeq plant FTP site (`https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/`).
  - Automatically locates the latest `annotation_releases` assembly and its `*_gene_ontology.gaf.gz` file, then converts it to the standard two-column format (`ID`, `Annotation`).
  - Lists all available plant species with `--list-species`; the converted file uses gene symbols (falling back to NCBI gene IDs when a symbol is missing).

- **GO annotation file format conversion**
  
  - Automatically identifies the gene ID column in Excel files and extracts all GO IDs.
  - Uses the built-in GO basic ontology file (`go-basic.obo`) to retrieve standard term names for GO IDs, outputting a two-column format (`ID`, `Annotation`) compatible with this tool.
  - Supports arbitrary column layouts and various delimiters (comma, pipe, semicolon, etc.) — no manual reformatting required.
  - **The GO ontology file is bundled with the package**, so you don't need to download or specify it manually.

- **GO enrichment analysis**
  
  - Computes hypergeometric test *p*-values and Benjamini-Hochberg correction from a local GO annotation Excel file.
  - Allows specification of ontology sub-ontologies (BP / CC / MF / ALL) and parses the GO basic ontology file (`go-basic.obo`) to obtain namespace information.
  - Output includes detailed information such as the list of associated genes (`All_index`).

- **Dynamic filtering and redundancy merging**
  
  - Automatically filters significant terms based on adjusted *p*-value, gene count, and median of negative log-*p*.
  - Merges redundant GO terms that share identical gene sets, keeping the term with the smallest *p*-value as the representative.

- **Similarity matrix and statistics**
  
  - Computes Jaccard or Ochiai coefficients to construct a similarity/distance matrix based on actual gene overlap between terms.
  - Supports multi-processing parallel computation (automatically enabled when the number of terms > 100).
  - Saves the similarity matrix to an Excel file and writes summary statistics (mean, median, 75th percentile, fraction of zero-similarity term pairs, and maximum) to a text file.

- **Dimensionality reduction and clustering**
  
  - Dimensionality reduction methods: UMAP (default), MDS, t-SNE.
  - Clustering algorithms: Mean-Shift (adaptive bandwidth), HDBSCAN.
  - Silhouette coefficient evaluation for clustering quality, automatic identification of boundary points.

- **Rich visualization outputs**
  
  - Cluster scatter plot, similarity heatmap, term-level bubble chart, static network (with per-cluster dashed ellipses and an external legend), and interactive network (Pyvis).
  - Cluster-level bubble chart: x-axis is the number of GO terms per cluster, bubble size is the union of genes in the cluster, and the y-axis is labelled with the 3–5 most frequent words of the cluster.
  - Cluster-specific word clouds (WordCloud).
  - Export of Cytoscape-compatible node and edge files.

- **Top GO term ranking**
  
  - Combines *p*-value and gene count into a single score — `score = -log10(Adjusted P-value) × log1p(Count_all)` — ranks GO terms, and exports the top 10 to an Excel file together with a classic enrichment bubble chart.

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

Core dependencies are installed automatically (see `pyproject.toml` for details):

- `pandas`, `numpy`, `scipy`, `statsmodels`, `requests`
- `matplotlib`, `seaborn`, `wordcloud`, `pyvis`
- `scikit-learn`, `umap-learn`, `hdbscan`
- `openpyxl`

## Quick Start

### 1. Prepare input files

- **GO annotation file** (Excel, e.g., `go_annotation.xlsx`):
  Must contain at least two columns: `ID` (gene identifier) and `Annotation` (GO annotations in the format `GO:xxxxxx,description`; multiple entries separated by semicolons).
  If your annotation file uses a different format, use the `convert` subcommand to standardise it (see below).
  
  If you don't have an annotation file yet, download one directly from NCBI RefSeq and convert it automatically (the result uses gene symbols):
  
  ```bash
  go-enrich download-go Solanum_lycopersicum --output go_annotation.xlsx
  ```

- **Target gene list file** (Excel, e.g., `id.xlsx`):
  The first column should contain gene IDs.

### 2. Convert any annotation file to the standard format (optional)

```bash
go-enrich convert any_format.xlsx standard.xlsx
```

The tool automatically uses the bundled `go-basic.obo` file to fill in GO term names.
If you want to use a custom ontology file, specify it with `--obo`:

```bash
go-enrich convert any_format.xlsx standard.xlsx --obo /path/to/your/go-basic.obo
```

The generated `standard.xlsx` can then be used directly for enrichment analysis.

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

After execution, multiple files are generated in the `results/` directory (see the "Output Files" section below).

## Command-line Arguments

### `enrich` subcommand

| Argument      | Required | Description                                                     |
| ------------- | -------- | --------------------------------------------------------------- |
| `--go-anno`   | Yes      | Path to the GO annotation Excel file (standard format).         |
| `--gene-list` | Yes      | Path to the target gene list Excel file.                        |
| `--output`    | No       | Output Excel file path (default: `go_enrichment_result.xlsx`).  |
| `--ontology`  | No       | Ontology sub-ontology: `BP`, `CC`, `MF`, `ALL` (default: `BP`). |
| `--use-obo`   | No       | Whether to use the obo file for namespace (enabled by default). |
| `--obo-path`  | No       | Path to the obo file (default uses the bundled file).           |

### `cluster` subcommand

| Argument                     | Description                                                                                                                                    |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`                    | Path to the enrichment result Excel file (must contain columns: `GO`, `Description`, `P-value`, `Adjusted P-value`, `Count_all`, `All_index`). |
| `--sheet`                    | Excel sheet name of the input file (default: `Sheet1`).                                                                                        |
| `--output-prefix`            | Prefix for output files (default: `go_clustering`).                                                                                            |
| **Filtering**                |                                                                                                                                                |
| `--p-cutoff`                 | Adjusted *p*-value cutoff (default: 0.05).                                                                                                     |
| `--min-genes`                | Minimum gene count (default: 3).                                                                                                               |
| `--logp-median-ratio`        | Whether to filter by median of negative log-*p* (enabled by default).                                                                          |
| `--keep-top-n`               | Number of top terms to keep when falling back automatically (default: 15).                                                                     |
| `--p-value-type`             | Which *p*-value to use for filtering: `adjusted` (Benjamini-Hochberg) or `raw` (default: `adjusted`).                                          |
| `--no-merge-identical`       | Do not merge terms with identical gene sets (merging is enabled by default).                                                                   |
| **Similarity**               |                                                                                                                                                |
| `--similarity`               | Similarity metric: `jaccard` or `ochiai` (default: `jaccard`).                                                                                 |
| `--parallel`                 | Enable parallel computation (enabled by default).                                                                                              |
| **Dimensionality reduction** |                                                                                                                                                |
| `--embed-method`             | Method: `umap`, `mds`, or `tsne` (default: `umap`).                                                                                            |
| `--random-state`             | Random seed (default: 42).                                                                                                                     |
| **Clustering**               |                                                                                                                                                |
| `--cluster-method`           | Algorithm: `meanshift` or `hdbscan` (default: `meanshift`).                                                                                    |
| `--k-neighbors`              | Number of neighbours for Mean-Shift bandwidth estimation (default: 8).                                                                         |
| `--min-bandwidth`            | Minimum bandwidth for Mean-Shift (default: 0.3).                                                                                               |
| `--hdbscan-min-cluster-size` | Minimum cluster size for HDBSCAN (default: 3).                                                                                                 |
| `--auto-fallback`            | Automatically adjust parameters when clustering fails (enabled by default).                                                                    |
| **Visualisation & export**   |                                                                                                                                                |
| `--similarity-threshold`     | Threshold for drawing edges in the network (default: 0.3).                                                                                     |
| `--no-plots`                 | Do not generate plots (only output the Excel tables and text summaries).                                                                       |
| `--wordcloud-dir`            | Directory for word clouds (default: `wordclouds`).                                                                                             |
| `--cytoscape`                | Export Cytoscape-compatible files (enabled by default).                                                                                        |

### `convert` subcommand

| Argument | Required | Description                                                                                                                 |
| -------- | -------- | --------------------------------------------------------------------------------------------------------------------------- |
| `input`  | Yes      | Input file path (supports arbitrary Excel formats).                                                                         |
| `output` | Yes      | Output file path (standard-format `.xlsx`).                                                                                 |
| `--obo`  | No       | Path to the GO basic ontology file. **By default, the tool uses the bundled `go-basic.obo` file**, so you rarely need this. |

### `download-go` subcommand

Downloads the GO annotation (GAF) file for a species from the NCBI RefSeq plant FTP site and converts it to the standard two-column format (`ID`, `Annotation`). The resulting `ID` column uses gene symbols (falling back to NCBI gene IDs when a symbol is missing).

| Argument         | Required | Description                                                                            |
| ---------------- | -------- | -------------------------------------------------------------------------------------- |
| `species`        | Yes      | Species Latin name with underscores, e.g. `Solanum_lycopersicum`.                      |
| `--output`       | No       | Output standard GO annotation Excel file (default: `go_annotation.xlsx`).              |
| `--obo`          | No       | Path to the GO basic ontology file (default uses the bundled file).                    |
| `--assembly`     | No       | Specific assembly directory name, e.g. `GCF_036512215.1-RS_2024_09` (default: latest). |
| `--keep-gaf`     | No       | Keep the downloaded GAF file after conversion.                                         |
| `--list-species` | No       | List all available plant species on NCBI RefSeq instead of downloading.                |

## Output Files

After running the `cluster` subcommand, the following files are generated using the output prefix:

| File type             | Example filename                  | Description                                                                                                                                                             |
| --------------------- | --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Similarity matrix     | `{prefix}_similarity_matrix.xlsx` | Similarity matrix (Jaccard/Ochiai) with GO IDs as row and column labels.                                                                                                |
| Similarity statistics | `{prefix}_similarity_stats.txt`   | Mean, median, 75th percentile, fraction of zero-similarity term pairs, and maximum similarity.                                                                          |
| Excel results         | `{prefix}_results.xlsx`           | Contains 4 sheets: `Clustering_Results` (detailed term results), `Cluster_Summary` (cluster summary), `Parameters` (run parameters), `Merged_Terms` (merged term info). |
| Top 10 GO terms       | `{prefix}_top10_go_terms.xlsx`    | Top 10 GO terms ranked by the combined score `-log10(Adjusted P-value) × log1p(Count_all)`.                                                                             |
| Scatter plot          | `{prefix}_scatter.png`            | 2D embedding coordinates with cluster colours.                                                                                                                          |
| Similarity heatmap    | `{prefix}_heatmap.png`            | Heatmap of the similarity matrix sorted by cluster.                                                                                                                     |
| Term bubble chart     | `{prefix}_bubble.png`             | Terms vs. -log10(*p*), bubble size represents gene count.                                                                                                               |
| Static network        | `{prefix}_network.png`            | Network of terms with similarity above the threshold; each cluster is circled by a dashed ellipse in its own colour and the legend is placed outside the axes.          |
| Cluster bubble chart  | `{prefix}_cluster_bubble.png`     | Cluster-level bubble chart: x = number of GO terms, bubble size = union of genes, y-axis labelled with 3–5 frequent words per cluster.                                  |
| Top 10 bubble chart   | `{prefix}_top10_bubble.png`       | Classic enrichment bubble chart of the top 10 terms (x = -log10(Adjusted P-value), bubble size = gene count, colour = significance).                                    |
| Interactive network   | `{prefix}_network.html`           | HTML interactive network (Pyvis) supporting drag, zoom, and hover details.                                                                                              |
| Word cloud directory  | `{prefix}_wordclouds/`            | One PNG word cloud per cluster.                                                                                                                                         |
| Cytoscape nodes       | `{prefix}_cytoscape_nodes.txt`    | Node information (tab-separated).                                                                                                                                       |
| Cytoscape edges       | `{prefix}_cytoscape_edges.txt`    | Edge list (GO_A, GO_B, Similarity).                                                                                                                                     |

## Complete Workflow Example

```bash
# 0. (optional) Download a GO annotation file from NCBI RefSeq
go-enrich download-go Solanum_lycopersicum --output standard.xlsx

# 1. Convert any annotation file (the tool uses bundled obo automatically)
go-enrich convert my_annotation.xlsx standard.xlsx

# 2. Enrichment analysis
go-enrich enrich --go-anno standard.xlsx --gene-list genes.xlsx --output enrich.xlsx --ontology ALL

# 3. Clustering with Ochiai similarity, t-SNE embedding, and HDBSCAN
go-enrich cluster --input enrich.xlsx --output-prefix my_clustering \
    --similarity ochiai --embed-method tsne --cluster-method hdbscan

# 4. Only output the tables and summaries, skip all plots
go-enrich cluster --input enrich.xlsx --output-prefix my_clustering_no_plots --no-plots
```

## Frequently Asked Questions

**Q1: Missing modules `umap` or `hdbscan`?**
A: Ensure all dependencies are installed: `pip install -e .` installs them automatically. Alternatively, install manually with `pip install umap-learn hdbscan`.

**Q2: Why are all terms clustered into one cluster?**
A: This may occur when the number of terms is too small or similarity is generally high. Try `--cluster-method hdbscan` or adjust `--min-bandwidth` / `--k-neighbors`.

**Q3: Some GO terms lack `GO_level_1` classification in the enrichment results?**
A: The tool uses the bundled `go-basic.obo` file by default to obtain namespace information. If you specify a custom `--obo-path`, make sure the path is correct and the file contains the namespace records. You can also download `go-basic.obo` manually and pass it with `--obo-path`.

**Q4: Descriptions are missing after running `convert`?**
A: The tool automatically uses the bundled `go-basic.obo` file to fill in descriptions. If descriptions are still missing, ensure that:

- The GO IDs in your input file are valid (e.g., `GO:0001234`).
- If you manually specified `--obo`, the file path is correct and contains the necessary mappings.

## Citation

If you use this tool in your research, please cite:

> [Sun Mintao]. *GO Enrichment and Clustering Tool*. GitHub repository. https://github.com/SunMintao/go-enrich-tool

## License

This project is licensed under the [MIT License](LICENSE). Feel free to use and modify it.
