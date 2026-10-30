#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GUI entry point for go-enrich-tool, built with Gooey.

Wraps the four subcommands (enrich / cluster / convert / download-go) into a
graphical interface. The underlying analysis reuses the same functions as the
command-line interface, so all outputs are identical.
"""

import sys
import os

from gooey import Gooey, GooeyParser

try:
    from .enrichment import run_enrichment
    from .core import run_clustering
    from .convert import convert_any_go_annotation
    from .ncbi_go import run_download_go, list_species
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from enrichment import run_enrichment
    from core import run_clustering
    from convert import convert_any_go_annotation
    from ncbi_go import run_download_go, list_species


def _build_enrich(subs):
    p = subs.add_parser("enrich", help="GO 富集分析")
    g = p.add_argument_group("输入文件")
    g.add_argument("--go-anno", required=True, widget="FileChooser",
                   help="GO 注释文件（Excel，含 ID 和 Annotation 两列）")
    g.add_argument("--gene-list", required=True, widget="FileChooser",
                   help="目标基因列表（Excel，第一列为基因标识）")
    g = p.add_argument_group("输出与选项")
    g.add_argument("--output", widget="FileSaver", default="go_enrichment_result.xlsx",
                   help="结果输出 Excel 文件")
    g.add_argument("--ontology", choices=["BP", "CC", "MF", "ALL"], default="BP",
                   help="GO 本体类型")
    g.add_argument("--use-obo", action="store_true", default=True,
                   help="使用 GO 本体文件获取分类信息")
    g.add_argument("--obo-path", default="go-basic.obo",
                   help="本地 obo 文件路径（默认使用内置文件）")


def _build_cluster(subs):
    p = subs.add_parser("cluster", help="GO 富集结果聚类")
    g = p.add_argument_group("输入与输出")
    g.add_argument("--input", required=True, widget="FileChooser",
                   help="GO 富集结果 Excel 文件")
    g.add_argument("--sheet", default="Sheet1", help="Excel 工作表名")
    g.add_argument("--output-prefix", default="go_clustering",
                   help="输出文件前缀（会生成多个文件）")

    g = p.add_argument_group("过滤参数")
    g.add_argument("--p-cutoff", type=float, default=0.05, widget="DecimalField",
                   help="显著性阈值（校正后 p 值）")
    g.add_argument("--min-genes", type=int, default=3, widget="IntegerField",
                   help="最小基因数")
    g.add_argument("--keep-top-n", type=int, default=15, widget="IntegerField",
                   help="自动回退时保留的 Top 术语数")
    g.add_argument("--p-value-type", choices=["adjusted", "raw"], default="adjusted",
                   help="使用哪种 p 值过滤：adjusted（校正后）或 raw（原始）")
    g.add_argument("--logp-median-ratio", action="store_true", default=True,
                   help="是否用负对数 p 值中位数进行过滤")
    g.add_argument("--no-merge-identical", action="store_true",
                   help="不合并基因集相同的术语")

    g = p.add_argument_group("相似度与降维")
    g.add_argument("--similarity", choices=["jaccard", "ochiai"], default="jaccard",
                   help="相似度度量")
    g.add_argument("--parallel", action="store_true", default=True,
                   help="启用并行计算")
    g.add_argument("--embed-method", choices=["umap", "mds", "tsne"], default="umap",
                   help="降维方法")
    g.add_argument("--random-state", type=int, default=42, widget="IntegerField",
                   help="随机种子")

    g = p.add_argument_group("聚类参数")
    g.add_argument("--cluster-method", choices=["meanshift", "hdbscan"], default="meanshift",
                   help="聚类算法")
    g.add_argument("--k-neighbors", type=int, default=8, widget="IntegerField",
                   help="Mean-Shift 带宽估计的邻居数")
    g.add_argument("--min-bandwidth", type=float, default=0.3, widget="DecimalField",
                   help="Mean-Shift 最小带宽")
    g.add_argument("--hdbscan-min-cluster-size", type=int, default=3, widget="IntegerField",
                   help="HDBSCAN 最小簇大小")
    g.add_argument("--auto-fallback", action="store_true", default=True,
                   help="聚类失败时自动调整参数")

    g = p.add_argument_group("可视化与导出")
    g.add_argument("--similarity-threshold", type=float, default=0.3, widget="DecimalField",
                   help="网络边相似度阈值")
    g.add_argument("--no-plots", action="store_true", help="不生成图（只输出表格）")
    g.add_argument("--wordcloud-dir", default="wordclouds", widget="DirChooser",
                   help="词云输出目录")
    g.add_argument("--cytoscape", action="store_true", default=True,
                   help="导出 Cytoscape 兼容文件")


def _build_convert(subs):
    p = subs.add_parser("convert", help="GO 注释格式转换")
    g = p.add_argument_group("输入与输出")
    g.add_argument("--input", required=True, widget="FileChooser",
                   help="输入文件（任意 Excel 格式的 GO 注释）")
    g.add_argument("--output", required=True, widget="FileSaver",
                   help="输出标准格式文件（.xlsx）")
    g.add_argument("--obo", default="go-basic.obo",
                   help="GO 本体文件路径（默认使用内置文件）")


def _build_download(subs):
    p = subs.add_parser("download-go", help="从 NCBI 下载 GO 注释")
    g = p.add_argument_group("物种与输出")
    g.add_argument("--species", default="",
                   help="物种拉丁名（下划线分隔），例如 Solanum_lycopersicum")
    g.add_argument("--output", widget="FileSaver", default="go_annotation.xlsx",
                   help="输出标准 GO 注释 Excel 文件")
    g.add_argument("--assembly", default="",
                   help="指定组装目录名（可选，例如 GCF_036512215.1-RS_2024_09）")
    g = p.add_argument_group("选项")
    g.add_argument("--obo", default="go-basic.obo",
                   help="GO 本体文件路径（默认使用内置文件）")
    g.add_argument("--keep-gaf", action="store_true",
                   help="转换后保留下载的 GAF 文件")
    g.add_argument("--list-species", action="store_true",
                   help="列出 NCBI RefSeq 上所有可用植物物种（不下载）")


@Gooey(
    program_name="go-enrich-tool",
    program_description="GO 富集分析与聚类工具（图形界面）",
    navigation="TABBED",
    tabbed_groups=True,
    clear_before_run=True,
    show_success_modal=True,
    show_failure_modal=True,
    default_size=(980, 740),
)
def main():
    parser = GooeyParser(description="GO enrichment analysis and clustering tool")
    subs = parser.add_subparsers(dest="command", help="选择功能")
    _build_enrich(subs)
    _build_cluster(subs)
    _build_convert(subs)
    _build_download(subs)
    args = parser.parse_args()

    if args.command == "enrich":
        run_enrichment(
            go_anno_file=args.go_anno,
            gene_list_file=args.gene_list,
            output_file=args.output,
            ontology=args.ontology,
            use_obo=args.use_obo,
            obo_path=args.obo_path,
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
            p_value_type=args.p_value_type,
        )
    elif args.command == "convert":
        convert_any_go_annotation(args.input, args.output, obo_path=args.obo)
    elif args.command == "download-go":
        if args.list_species:
            for s in list_species():
                print(s)
        else:
            if not args.species:
                print("请填写物种拉丁名（例如 Solanum_lycopersicum），或勾选 --list-species 查看可用物种。")
                return
            run_download_go(
                species=args.species,
                output_file=args.output,
                obo_path=args.obo,
                assembly=args.assembly or None,
                keep_gaf=args.keep_gaf,
            )


if __name__ == "__main__":
    main()
