#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO富集分析及聚类工具的命令行入口
用法示例:
  go-enrich enrich --go-anno go_annotation.xlsx --gene-list id.xlsx --output enrich_result.xlsx
  go-enrich cluster --input enrich_result.xlsx --output-prefix clustering_out
"""

import argparse
import sys
import os

try:
    from go_enrich.enrichment import run_enrichment
    from go_enrich.core import run_clustering
except ImportError:
    # 如果尚未打包，允许直接运行当前目录下的模块（开发模式）
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from enrichment import run_enrichment
    from core import run_clustering


def main():
    parser = argparse.ArgumentParser(
        description="GO富集分析与聚类工具",
        prog="go-enrich"
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="子命令")

    # ---------- enrich 子命令 ----------
    enrich_parser = subparsers.add_parser("enrich", help="运行GO富集分析")
    enrich_parser.add_argument("--go-anno", required=True, help="GO注释文件 (Excel格式，含ID和Annotation列)")
    enrich_parser.add_argument("--gene-list", required=True, help="目标基因列表文件 (Excel格式，第一列为基因ID)")
    enrich_parser.add_argument("--output", default="go_enrichment_result.xlsx", help="输出Excel文件路径")
    enrich_parser.add_argument("--ontology", default="BP", choices=["BP", "CC", "MF", "ALL"],
                               help="GO本体类型，默认BP")
    enrich_parser.add_argument("--use-obo", action="store_true", default=True,
                               help="使用GO基础本体文件获取分类信息")
    enrich_parser.add_argument("--obo-path", default="go-basic.obo", help="本地obo文件路径或下载保存位置")

    # ---------- cluster 子命令 ----------
    cluster_parser = subparsers.add_parser("cluster", help="对GO富集结果进行聚类分析")
    cluster_parser.add_argument("--input", required=True,
                                help="GO富集结果Excel文件 (包含GO, Description, P-value, Adjusted P-value, Count_all, All_index等列)")
    cluster_parser.add_argument("--sheet", default="Sheet1", help="Excel工作表名")
    cluster_parser.add_argument("--output-prefix", default="go_clustering", help="输出文件前缀 (将生成多个文件)")

    # 筛选参数
    cluster_parser.add_argument("--p-cutoff", type=float, default=0.05, help="校正后p值阈值")
    cluster_parser.add_argument("--min-genes", type=int, default=3, help="最小基因计数")
    cluster_parser.add_argument("--logp-median-ratio", action="store_true", default=True,
                                help="是否使用负对数p值中位数进行筛选")
    cluster_parser.add_argument("--keep-top-n", type=int, default=15, help="自动降级时保留的top术语数")
    cluster_parser.add_argument("--no-merge-identical", action="store_true", help="不合并基因集完全相同的术语")

    # 相似度参数
    cluster_parser.add_argument("--similarity", default="jaccard", choices=["jaccard", "ochiai"], help="相似度度量")
    cluster_parser.add_argument("--parallel", action="store_true", default=True, help="并行计算相似度矩阵")

    # 降维参数
    cluster_parser.add_argument("--embed-method", default="umap", choices=["umap", "mds", "tsne"], help="降维方法")
    cluster_parser.add_argument("--random-state", type=int, default=42, help="随机种子")

    # 聚类参数
    cluster_parser.add_argument("--cluster-method", default="meanshift", choices=["meanshift", "hdbscan"],
                                help="聚类算法")
    cluster_parser.add_argument("--k-neighbors", type=int, default=8, help="Mean-Shift带宽估计的k近邻数")
    cluster_parser.add_argument("--min-bandwidth", type=float, default=0.3, help="Mean-Shift最小带宽")
    cluster_parser.add_argument("--hdbscan-min-cluster-size", type=int, default=3, help="HDBSCAN最小簇大小")
    cluster_parser.add_argument("--auto-fallback", action="store_true", default=True, help="聚类结果异常时自动调整参数")

    # 可视化与导出
    cluster_parser.add_argument("--similarity-threshold", type=float, default=0.3, help="网络图连边相似度阈值")
    cluster_parser.add_argument("--no-plots", action="store_true", help="不生成绘图（仅输出表格）")
    cluster_parser.add_argument("--wordcloud-dir", default="wordclouds", help="词云输出目录")
    cluster_parser.add_argument("--cytoscape", action="store_true", default=True, help="导出Cytoscape兼容文件")

    args = parser.parse_args()

    # 执行相应子命令
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
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()