#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Command-line entry point for the GO enrichment and clustering tool.
Usage examples:
  go-enrich enrich --go-anno go_annotation.xlsx --gene-list id.xlsx --output enrich_result.xlsx
  go-enrich cluster --input enrich_result.xlsx --output-prefix clustering_out
  go-enrich convert any_annotation.xlsx standard.xlsx [--obo go-basic.obo]

  # default: Adjusted P-value
  go-enrich cluster --input enrich.xlsx --p-cutoff 0.01
  # p-value
  go-enrich cluster --input enrich.xlsx --p-cutoff 0.001 --p-value-type raw
"""

import argparse
import sys
import os

try:
    from go_enrich.enrichment import run_enrichment
    from go_enrich.core import run_clustering
    from go_enrich.convert import convert_any_go_annotation
    from go_enrich.ncbi_go import run_download_go, list_species
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from enrichment import run_enrichment
    from core import run_clustering
    from convert import convert_any_go_annotation
    from ncbi_go import run_download_go, list_species


def main():
    parser = argparse.ArgumentParser(
        description="GO enrichment analysis and clustering tool",
        prog="go-enrich"
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="sub-commands")

    # ---------- enrich sub-command ----------
    enrich_parser = subparsers.add_parser("enrich", help="Run GO enrichment analysis")
    enrich_parser.add_argument("--go-anno", required=True, help="GO annotation file (Excel format, with ID and Annotation columns)")
    enrich_parser.add_argument("--gene-list", required=True, help="Target gene list file (Excel format, first column gene IDs)")
    enrich_parser.add_argument("--output", default="go_enrichment_result.xlsx", help="Output Excel file path")
    enrich_parser.add_argument("--ontology", default="BP", choices=["BP", "CC", "MF", "ALL"],
                               help="GO ontology type, default BP")
    enrich_parser.add_argument("--use-obo", action="store_true", default=True,
                               help="Use GO basic ontology file to obtain classification information")
    enrich_parser.add_argument("--obo-path", default="go-basic.obo",
                               help="Local obo file path or download location (default uses built-in file)")

    # ---------- cluster sub-command ----------
    cluster_parser = subparsers.add_parser("cluster", help="Cluster GO enrichment results")
    cluster_parser.add_argument("--input", required=True,
                                help="GO enrichment result Excel file (must contain GO, Description, P-value, Adjusted P-value, Count_all, All_index columns)")
    cluster_parser.add_argument("--sheet", default="Sheet1", help="Excel sheet name")
    cluster_parser.add_argument("--output-prefix", default="go_clustering", help="Output file prefix (multiple files will be generated)")

    # Filtering parameters
    cluster_parser.add_argument("--p-cutoff", type=float, default=0.05, help="Adjusted p-value threshold")
    cluster_parser.add_argument("--min-genes", type=int, default=3, help="Minimum gene count")
    cluster_parser.add_argument("--logp-median-ratio", action="store_true", default=True,
                                help="Whether to use median of negative log p-value for filtering")
    cluster_parser.add_argument("--keep-top-n", type=int, default=15, help="Number of top terms to keep when falling back automatically")
    cluster_parser.add_argument("--no-merge-identical", action="store_true", help="Do not merge terms with identical gene sets")
    # 新增：选择 p-value 类型
    cluster_parser.add_argument("--p-value-type", choices=["adjusted", "raw"], default="adjusted",
                                help="Which p-value to use for filtering: 'adjusted' (Benjamini-Hochberg) or 'raw' (nominal). Default: adjusted")

    # Similarity parameters
    cluster_parser.add_argument("--similarity", default="jaccard", choices=["jaccard", "ochiai"], help="Similarity metric")
    cluster_parser.add_argument("--parallel", action="store_true", default=True, help="Enable parallel computation")

    # Dimensionality reduction parameters
    cluster_parser.add_argument("--embed-method", default="umap", choices=["umap", "mds", "tsne"], help="Dimensionality reduction method")
    cluster_parser.add_argument("--random-state", type=int, default=42, help="Random seed")

    # Clustering parameters
    cluster_parser.add_argument("--cluster-method", default="meanshift", choices=["meanshift", "hdbscan"],
                                help="Clustering algorithm")
    cluster_parser.add_argument("--k-neighbors", type=int, default=8, help="Number of neighbors for Mean-Shift bandwidth estimation")
    cluster_parser.add_argument("--min-bandwidth", type=float, default=0.3, help="Minimum bandwidth for Mean-Shift")
    cluster_parser.add_argument("--hdbscan-min-cluster-size", type=int, default=3, help="Minimum cluster size for HDBSCAN")
    cluster_parser.add_argument("--auto-fallback", action="store_true", default=True, help="Automatically adjust parameters if clustering fails")

    # Visualization and export
    cluster_parser.add_argument("--similarity-threshold", type=float, default=0.3, help="Similarity threshold for network edges")
    cluster_parser.add_argument("--no-plots", action="store_true", help="Do not generate plots (only output table)")
    cluster_parser.add_argument("--wordcloud-dir", default="wordclouds", help="Word cloud output directory")
    cluster_parser.add_argument("--cytoscape", action="store_true", default=True, help="Export Cytoscape-compatible files")

    # ---------- convert sub-command ----------
    convert_parser = subparsers.add_parser("convert", help="Convert any GO annotation file to standard format (using obo file for descriptions)")
    convert_parser.add_argument("input", help="Input file path (supports arbitrary Excel formats)")
    convert_parser.add_argument("output", help="Output file path (standard format .xlsx)")
    convert_parser.add_argument("--obo", default="go-basic.obo",
                                help="GO basic ontology file path (default uses built-in file)")

    # ---------- download-go sub-command ----------
    download_parser = subparsers.add_parser("download-go", help="Download GO annotation (GAF) from NCBI RefSeq and convert to standard format")
    download_parser.add_argument("species", nargs='?',
                                 help="Species Latin name, e.g. Solanum_lycopersicum")
    download_parser.add_argument("--output", default="go_annotation.xlsx",
                                 help="Output standard GO annotation Excel file (default: go_annotation.xlsx)")
    download_parser.add_argument("--obo", default="go-basic.obo",
                                 help="GO basic ontology file path (default uses built-in file)")
    download_parser.add_argument("--assembly",
                                 help="Specific assembly directory name, e.g. GCF_036512215.1-RS_2024_09")
    download_parser.add_argument("--keep-gaf", action="store_true",
                                 help="Keep the downloaded GAF file after conversion")
    download_parser.add_argument("--list-species", action="store_true",
                                 help="List all available plant species on NCBI RefSeq")

    args = parser.parse_args()

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
            export_cytoscape=args.cytoscape,
            p_value_type=args.p_value_type   # 传递
        )
    elif args.command == "convert":
        convert_any_go_annotation(args.input, args.output, obo_path=args.obo)
    elif args.command == "download-go":
        if args.list_species:
            for s in list_species():
                print(s)
        else:
            if not args.species:
                parser.error("Please provide a species name (e.g. Solanum_lycopersicum) or use --list-species")
            run_download_go(
                species=args.species,
                output_file=args.output,
                obo_path=args.obo,
                assembly=args.assembly,
                keep_gaf=args.keep_gaf,
            )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()