```markdown
# GO Enrichment and Clustering Tool

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**GO Enrichment and Clustering Tool** 是一个用于基因本体（GO）富集分析及结果聚类的一体化命令行工具。它利用超几何检验进行富集分析，并基于基因重叠相似度（Jaccard / Ochiai）对富集得到的 GO 术语进行降维（UMAP / MDS / t‑SNE）和无监督聚类（Mean‑Shift / HDBSCAN），最终生成多种可视化图表及网络文件，帮助研究人员从全局视角挖掘差异基因的功能模块。

## 主要功能

- **GO富集分析**  
  - 支持从本地 GO 注释文件（Excel）自动计算超几何检验 \(p\) 值及 Benjamini‑Hochberg 校正。  
  - 可指定本体类型（BP / CC / MF / ALL），自动下载并解析 GO 基础本体文件（go‑basic.obo）获取 namespace。  
  - 输出包含差异基因列表（`All_index`）等详细信息。

- **动态筛选与冗余合并**  
  - 根据校正后 \(p\) 值、基因计数、负对数中位数自动筛选显著术语。  
  - 自动合并基因集完全相同的冗余 GO 条目，保留 \(p\) 值最小的代表术语。

- **相似度矩阵**  
  - 计算 Jaccard 或 Ochiai 系数，构建术语间基于真实基因重叠的相似性/距离矩阵。  
  - 支持多进程并行加速（术语数 > 100 时自动启用）。

- **降维与聚类**  
  - 降维方法：UMAP（默认）、MDS、t‑SNE。  
  - 聚类算法：Mean‑Shift（自适应带宽）、HDBSCAN。  
  - 轮廓系数评估聚类质量，自动标记边界点。

- **丰富的可视化输出**  
  - 聚类散点图、相似性热图、气泡图、静态/交互式网络图（Pyvis）。  
  - 按簇生成功能词云（WordCloud）。  
  - 导出 Cytoscape 兼容的节点和边文件。

## 安装

### 从源码安装

```bash
# 克隆仓库
git clone https://github.com/SunMintao/go-enrich-tool.git
cd go-enrich-tool

# 推荐使用虚拟环境
python -m venv venv
source venv/bin/activate   # Linux/macOS
venv\Scripts\activate      # Windows

# 安装包
pip install -e .
```

### 依赖

自动安装的核心依赖包括（详见 `pyproject.toml`）：

- `pandas`, `numpy`, `scipy`, `statsmodels`, `requests`
- `matplotlib`, `seaborn`, `wordcloud`, `pyvis`
- `scikit-learn`, `umap-learn`, `hdbscan`
- `openpyxl`

## 快速开始

### 1. 准备输入文件

- **GO注释文件**（Excel，例如 `go_annotation.xlsx`）：  
  至少包含两列：`ID`（基因标识符）和 `Annotation`（GO 注释，格式为 `GO:xxxxxx,描述`；多个条目以分号分隔）。
- **目标基因列表文件**（Excel，例如 `id.xlsx`）：  
  第一列为基因 ID。

### 2. 运行富集分析

```bash
go-enrich enrich \
    --go-anno data/go_annotation.xlsx \
    --gene-list data/id.xlsx \
    --output results/enrichment.xlsx \
    --ontology BP
```

### 3. 对富集结果进行聚类

```bash
go-enrich cluster \
    --input results/enrichment.xlsx \
    --output-prefix results/clustering_output
```

执行后会在 `results/` 目录下生成多个文件（详见下方“输出文件说明”）。

## 命令行参数详情

### `enrich` 子命令

| 参数            | 必选  | 说明                                         |
| ------------- | --- | ------------------------------------------ |
| `--go-anno`   | 是   | GO 注释 Excel 文件路径                           |
| `--gene-list` | 是   | 目标基因列表 Excel 文件路径                          |
| `--output`    | 否   | 富集结果输出文件路径（默认 `go_enrichment_result.xlsx`） |
| `--ontology`  | 否   | 本体类型，可选 `BP`, `CC`, `MF`, `ALL`（默认 `BP`）   |
| `--use-obo`   | 否   | 是否下载并使用 obo 文件获取 namespace（默认启用）           |
| `--obo-path`  | 否   | obo 文件保存/读取路径（默认 `go-basic.obo`）           |

### `cluster` 子命令

| 参数                           | 说明                                                                                                    |
| ---------------------------- | ----------------------------------------------------------------------------------------------------- |
| `--input`                    | 富集结果 Excel 文件路径（必须包含 `GO`, `Description`, `P-value`, `Adjusted P-value`, `Count_all`, `All_index` 等列） |
| `--output-prefix`            | 输出文件前缀（默认 `go_clustering`）                                                                            |
| **筛选参数**                     |                                                                                                       |
| `--p-cutoff`                 | 校正后 \(p\) 值阈值（默认 0.05）                                                                                |
| `--min-genes`                | 最小基因计数（默认 3）                                                                                          |
| `--logp-median-ratio`        | 是否使用负对数 \(p\) 值中位数过滤（默认启用）                                                                            |
| `--keep-top-n`               | 自动降级时保留的 top 术语数量（默认 15）                                                                              |
| `--no-merge-identical`       | 不合并基因集完全相同的术语（默认合并）                                                                                   |
| **相似度参数**                    |                                                                                                       |
| `--similarity`               | 相似度度量，`jaccard` 或 `ochiai`（默认 `jaccard`）                                                              |
| `--parallel`                 | 是否并行计算（默认启用）                                                                                          |
| **降维参数**                     |                                                                                                       |
| `--embed-method`             | 降维方法，`umap`, `mds`, `tsne`（默认 `umap`）                                                                 |
| `--random-state`             | 随机种子（默认 42）                                                                                           |
| **聚类参数**                     |                                                                                                       |
| `--cluster-method`           | 聚类算法，`meanshift` 或 `hdbscan`（默认 `meanshift`）                                                          |
| `--k-neighbors`              | Mean‑Shift 带宽估计的 \(k\) 近邻数（默认 8）                                                                      |
| `--min-bandwidth`            | Mean‑Shift 最小带宽（默认 0.3）                                                                               |
| `--hdbscan-min-cluster-size` | HDBSCAN 最小簇大小（默认 3）                                                                                   |
| `--auto-fallback`            | 聚类异常时自动调整参数（默认启用）                                                                                     |
| **可视化与导出**                   |                                                                                                       |
| `--similarity-threshold`     | 网络图连边相似度阈值（默认 0.3）                                                                                    |
| `--no-plots`                 | 不生成绘图（仅输出表格）                                                                                          |
| `--wordcloud-dir`            | 词云输出目录（默认 `wordclouds`）                                                                               |
| `--cytoscape`                | 是否导出 Cytoscape 文件（默认启用）                                                                               |

## 输出文件说明

运行 `cluster` 子命令后，将在输出前缀所在目录生成以下文件：

| 文件类型         | 文件名示例                          | 说明                                                                                                       |
| ------------ | ------------------------------ | -------------------------------------------------------------------------------------------------------- |
| Excel 结果     | `{prefix}_results.xlsx`        | 包含 4 个工作表：`Clustering_Results`（每个术语的详细结果）、`Cluster_Summary`（簇汇总）、`Parameters`（运行参数）、`Merged_Terms`（合并信息） |
| 散点图          | `{prefix}_scatter.png`         | 二维降维坐标及簇着色                                                                                               |
| 相似性热图        | `{prefix}_heatmap.png`         | 按簇排序的相似度矩阵热图                                                                                             |
| 气泡图          | `{prefix}_bubble.png`          | 术语描述 vs \(-\log_{10}(p)\)，气泡大小表示基因计数                                                                     |
| 静态网络图        | `{prefix}_network.png`         | 相似度超过阈值的术语连接网络                                                                                           |
| 交互式网络        | `{prefix}_network.html`        | 基于 Pyvis 的 HTML 交互网络，支持拖拽、悬停查看信息                                                                         |
| 词云目录         | `{prefix}_wordclouds/`         | 每个簇单独的词云图（PNG）                                                                                           |
| Cytoscape 节点 | `{prefix}_cytoscape_nodes.txt` | 节点信息表（制表符分隔）                                                                                             |
| Cytoscape 边  | `{prefix}_cytoscape_edges.txt` | 边列表（GO_A, GO_B, Similarity）                                                                              |

## 使用示例数据运行富集分析

```
go-enrich enrich --go-anno examples/data/go_annotation.xlsx --gene-list examples/data/id.xlsx --output enrich.xlsx
```

## 完整工作流

```bash
# 1. 富集分析
go-enrich enrich --go-anno tests/go_annotation.xlsx --gene-list tests/id.xlsx --output enrich.xlsx

# 2. 聚类（使用 Ochiai 相似度、t-SNE 降维、HDBSCAN 聚类）
go-enrich cluster --input enrich.xlsx --output-prefix my_clustering \
    --similarity ochiai --embed-method tsne --cluster-method hdbscan

# 3. 仅输出表格，不绘图
go-enrich cluster --input enrich.xlsx --output-prefix my_clustering_no_plots --no-plots
```

## 常见问题

**Q1: 提示缺少模块 `umap` 或 `hdbscan`？**  
A: 请确保已安装所有依赖：`pip install -e .` 会自动安装；或手动安装 `pip install umap-learn hdbscan`。

**Q2: 为什么聚类后所有术语被分到一簇？**  
A: 可能原因是术语数量过少或相似度普遍较高。请尝试使用 `--cluster-method hdbscan` 或调整 `--min-bandwidth` / `--k-neighbors` 参数。

**Q3: 富集分析中某些 GO 术语没有 `GO_level_1` 分类？**  
A: 请确认 `--use-obo` 已启用且网络通畅，工具会自动下载 `go-basic.obo`。若仍失败，可手动下载该文件并指定 `--obo-path`。

## 引用

如果您在研究中使用了本工具，请引用：

> [Sun Mintao]. *GO Enrichment and Clustering Tool*. GitHub repository. https://github.com/SunMintao/go-enrich-tool

## 许可证

本项目采用 [MIT 许可证](LICENSE)，欢迎自由使用和修改。
