#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GUI entry point for GO enrichment tool using Gooey.
Run with: python -m go_enrich.gui
or after install: go-enrich-gui
"""

import sys
from gooey import Gooey, GooeyParser

# Import core functions (same as cli.py)
from go_enrich.enrichment import run_enrichment
from go_enrich.core import run_clustering
from go_enrich.convert import convert_any_go_annotation


@Gooey(
    program_name="GO Enrichment Tool",
    program_description="Enrichment, conversion, and clustering of GO terms",
    clear_before_run=True,
    default_size=(1000, 800),
    optional_cols=2,
    required_cols=1,
    auto_start=False,
    show_success_modal=True,
    load_config_from='',
)
def main():
    parser = GooeyParser(description="GO enrichment analysis and clustering tool")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Choose a function")

    # ---------- enrich sub-command ----------
    enrich_parser = subparsers.add_parser("enrich", help="Run GO enrichment analysis")
    enrich_parser.add_argument(
        "--go-anno",
        required=True,
        widget="FileChooser",
        help="GO annotation file (Excel, with ID and Annotation columns)"
    )
    enrich_parser.add_argument(
        "--gene-list",
        required=True,
        widget="FileChooser",
        help="Target gene list file (Excel, first column gene IDs)"
    )
    enrich_parser.add_argument(
        "--output",
        default="go_enrichment_result.xlsx",
        widget="FileSaver",
        help="Output Excel file path"
    )
    enrich_parser.add_argument(
        "--ontology",
        default="BP",
        choices=["BP", "CC", "MF", "ALL"],
        widget="Dropdown",
        help="GO ontology type"
    )
    enrich_parser.add_argument(
        "--use-obo",
        action="store_true",
        default=True,
        help="Use obo file for namespace (enabled by default)"
    )
    enrich_parser.add_argument(
        "--obo-path",
        default="go-basic.obo",
        widget="FileChooser",
        help="Local obo file path (default uses built-in file)"
    )

    # ---------- cluster sub-command ----------
    cluster_parser = subparsers.add_parser("cluster", help="Cluster GO enrichment results")
    cluster_parser.add_argument(
        "--input",
        required=True,
        widget="FileChooser",
        help="GO enrichment result Excel file (must contain required columns)"
    )
    cluster_parser.add_argument(
        "--sheet",
        default="Sheet1",
        help="Excel sheet name"
    )
    cluster_parser.add_argument(
        "--output-prefix",
        default="go_clustering",
        help="Prefix for output files"
    )

    # Filtering parameters
    cluster_parser.add_argument(
        "--p-cutoff",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold"
    )
    cluster_parser.add_argument(
        "--min-genes",
        type=int,
        default=3,
        help="Minimum gene count"
    )
    cluster_parser.add_argument(
        "--logp-median-ratio",
        action="store_true",
        default=True,
        help="Use median of negative log p-value for filtering"
    )
    cluster_parser.add_argument(
        "--keep-top-n",
        type=int,
        default=15,
        help="Number of top terms to keep when falling back automatically"
    )
    cluster_parser.add_argument(
        "--no-merge-identical",
        action="store_true",
        help="Do not merge terms with identical gene sets"
    )

    # Similarity parameters
    cluster_parser.add_argument(
        "--similarity",
        default="jaccard",
        choices=["jaccard", "ochiai"],
        widget="Dropdown",
        help="Similarity metric"
    )
    cluster_parser.add_argument(
        "--parallel",
        action="store_true",
        default=True,
        help="Enable parallel computation"
    )

    # Dimensionality reduction
    cluster_parser.add_argument(
        "--embed-method",
        default="umap",
        choices=["umap", "mds", "tsne"],
        widget="Dropdown",
        help="Dimensionality reduction method"
    )
    cluster_parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed"
    )

    # Clustering parameters
    cluster_parser.add_argument(
        "--cluster-method",
        default="meanshift",
        choices=["meanshift", "hdbscan"],
        widget="Dropdown",
        help="Clustering algorithm"
    )
    cluster_parser.add_argument(
        "--k-neighbors",
        type=int,
        default=8,
        help="Number of neighbors for Mean-Shift bandwidth estimation"
    )
    cluster_parser.add_argument(
        "--min-bandwidth",
        type=float,
        default=0.3,
        help="Minimum bandwidth for Mean-Shift"
    )
    cluster_parser.add_argument(
        "--hdbscan-min-cluster-size",
        type=int,
        default=3,
        help="Minimum cluster size for HDBSCAN"
    )
    cluster_parser.add_argument(
        "--auto-fallback",
        action="store_true",
        default=True,
        help="Automatically adjust parameters if clustering fails"
    )

    # Visualization and export
    cluster_parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.3,
        help="Similarity threshold for network edges"
    )
    cluster_parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Do not generate plots (only output table)"
    )
    cluster_parser.add_argument(
        "--wordcloud-dir",
        default="wordclouds",
        help="Word cloud output directory"
    )
    cluster_parser.add_argument(
        "--cytoscape",
        action="store_true",
        default=True,
        help="Export Cytoscape-compatible files"
    )

    # ---------- convert sub-command ----------
    convert_parser = subparsers.add_parser("convert", help="Convert annotation file to standard format")
    convert_parser.add_argument(
        "input",
        widget="FileChooser",
        help="Input file (any Excel format)"
    )
    convert_parser.add_argument(
        "output",
        widget="FileSaver",
        help="Output standard file (.xlsx)"
    )
    convert_parser.add_argument(
        "--obo",
        default="go-basic.obo",
        widget="FileChooser",
        help="obo file (default uses built-in file)"
    )

    args = parser.parse_args()

    # Execute the corresponding sub-command
    if args.command == "enrich":
        run_enrichment(
            go_anno_file=args.go_anno,
            gene_list_file=args.gene_list,
            output_file=args.output,
            ontology=args.ontology,
            use_obo=args.use_obo,
            obo_path=args.obo_path
        )
    elif args.command == "cluster":
        run_clustering(
            input_file=args.input,
            sheet_name=args.sheet,
            output_prefix=args.output_prefix,
            p_cutoff=args.p_cutoff,
            min_genes=args.min_genes,
            logp_median_ratio=args.logp_median_ratio,
            keep_top_n=args.keep_top_n,
            merge_identical=not args.no_merge_identical,
            similarity_metric=args.similarity,
            parallel=args.parallel,
            embed_method=args.embed_method,
            random_state=args.random_state,
            cluster_method=args.cluster_method,
            k_neighbors=args.k_neighbors,
            min_bandwidth=args.min_bandwidth,
            hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
            auto_fallback=args.auto_fallback,
            similarity_threshold=args.similarity_threshold,
            no_plots=args.no_plots,
            wordcloud_dir=args.wordcloud_dir,
            export_cytoscape=args.cytoscape
        )
    elif args.command == "convert":
        convert_any_go_annotation(args.input, args.output, obo_path=args.obo)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()