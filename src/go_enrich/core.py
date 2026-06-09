#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO 富集结果聚类分析（强制依赖版）
- 必须安装所有依赖包，否则报错
- 包含动态筛选、冗余合并、Jaccard/Ochiai相似度、UMAP降维、Mean‑Shift/HDBSCAN聚类、
  轮廓系数评估、静态/交互式网络图、词云、Cytoscape导出等全部功能
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS, TSNE
from sklearn.cluster import MeanShift
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.neighbors import NearestNeighbors
import umap
import hdbscan
from pyvis.network import Network
from wordcloud import WordCloud
import warnings
import re
from collections import Counter
import os
import itertools

plt.rcParams['font.family'] = 'Arial'
warnings.filterwarnings('ignore')


class GOTermClusteringAnalyzer:
    def __init__(self, file_path, sheet_name='Sheet1'):
        self.file_path = file_path
        self.sheet_name = sheet_name
        self.df = None
        self.gene_sets = None
        self.similarity_matrix = None
        self.distance_matrix = None
        self.coords_2d = None
        self.cluster_labels = None
        self.cluster_centers = None
        self.results_df = None
        self.cluster_summary = None
        self.silhouette_avg = None
        self.parameters = {}
        self.merged_terms_info = []

    # ---------- 1. 数据加载与预处理 ----------
    def load_data(self):
        print(f"加载数据: {self.file_path}")
        self.df = pd.read_excel(self.file_path, sheet_name=self.sheet_name)

        def parse_genes(gene_string):
            if pd.isna(gene_string) or gene_string == '':
                return frozenset()
            return frozenset(gene_string.strip().split(';'))

        self.df['Gene_Set'] = self.df['All_index'].apply(parse_genes)
        self.df = self.df[self.df['Gene_Set'].apply(len) > 0].reset_index(drop=True)
        print(f"原始 GO term 数量: {len(self.df)}")

        empty_ratio = (self.df['All_index'].isna() | (self.df['All_index'] == '')).mean()
        if empty_ratio > 0.3:
            print(f"警告: {empty_ratio:.1%} 的 All_index 为空，请检查输入文件。")

        mismatch = (self.df['Count_all'] != self.df['Gene_Set'].apply(len)).sum()
        if mismatch > 0:
            print(f"警告: {mismatch} 行的 Count_all 与基因集实际大小不一致，已自动校准。")
            self.df['Count_all'] = self.df['Gene_Set'].apply(len)

        if self.df['GO'].duplicated().any():
            self.df = self.df.sort_values('P-value').drop_duplicates('GO').reset_index(drop=True)
            print("已去除重复的 GO ID。")
        return self

    def filter_terms_dynamic(self, p_cutoff=0.05, min_genes=5, logp_median_ratio=True,
                             keep_top_n=15, auto_fallback=True):
        print("动态筛选 GO 术语...")
        mask = (self.df['Adjusted P-value'] < p_cutoff) & (self.df['Count_all'] >= min_genes)
        if logp_median_ratio:
            neg_logp = -np.log10(self.df['P-value'].clip(lower=1e-300))
            median_logp = neg_logp.median()
            mask = mask & (neg_logp > median_logp)
        keep_idx = self.df[mask].index
        if len(keep_idx) < 3 and auto_fallback:
            mask2 = self.df['Adjusted P-value'] < p_cutoff
            keep_idx = self.df[mask2].index
            if len(keep_idx) < 3:
                keep_idx = self.df.nsmallest(keep_top_n, 'P-value').index
                print(f"显著术语不足，按 P-value 取 top {keep_top_n}")
        self.df = self.df.loc[keep_idx].reset_index(drop=True)
        print(f"筛选后保留 {len(self.df)} 个术语")
        self.parameters['filter'] = {
            'p_cutoff': p_cutoff, 'min_genes': min_genes,
            'logp_median_ratio': logp_median_ratio,
            'keep_top_n': keep_top_n, 'final_count': len(self.df)
        }
        return self

    def merge_identical_sets(self):
        print("检查并合并完全相同基因集...")
        set_to_indices = {}
        for idx, row in self.df.iterrows():
            fs = row['Gene_Set']
            set_to_indices.setdefault(fs, []).append(idx)
        to_keep = []
        self.merged_terms_info = []
        for fs, indices in set_to_indices.items():
            if len(indices) == 1:
                to_keep.append(indices[0])
            else:
                sub_df = self.df.loc[indices]
                best_idx = sub_df['P-value'].idxmin()
                to_keep.append(best_idx)
                others = [i for i in indices if i != best_idx]
                merged_ids = sub_df.loc[others, 'GO'].tolist()
                self.merged_terms_info.append({
                    'representative_GO': sub_df.loc[best_idx, 'GO'],
                    'merged_GOs': ','.join(merged_ids),
                    'description': sub_df.loc[best_idx, 'Description']
                })
                print(f"  合并 {len(others)} 个术语: {sub_df.loc[best_idx, 'GO']} 保留了 {merged_ids}")
        if to_keep:
            self.df = self.df.loc[to_keep].reset_index(drop=True)
        self.df['Merged_GOs'] = ''
        for info in self.merged_terms_info:
            rep_go = info['representative_GO']
            self.df.loc[self.df['GO'] == rep_go, 'Merged_GOs'] = info['merged_GOs']
        print(f"合并后剩余 {len(self.df)} 个非冗余术语")
        self.parameters['merged'] = {'n_merged_groups': len(self.merged_terms_info)}
        return self

    # ---------- 2. 相似度计算 ----------
    def calculate_similarity(self, metric='jaccard', parallel=True):
        n = len(self.df)
        if n == 0:
            raise ValueError("无有效 GO term")
        print(f"计算 {metric} 相似度矩阵 (n={n})...")
        self.gene_sets = [row['Gene_Set'] for _, row in self.df.iterrows()]

        def sim_pair(i, j):
            si = self.gene_sets[i]
            sj = self.gene_sets[j]
            inter = len(si & sj)
            if metric == 'jaccard':
                union = len(si | sj)
                return inter / union if union > 0 else 0.0
            elif metric == 'ochiai':
                denom = np.sqrt(len(si) * len(sj))
                return inter / denom if denom > 0 else 0.0
            else:
                raise ValueError("metric must be 'jaccard' or 'ochiai'")

        self.similarity_matrix = np.eye(n)
        if parallel and n > 100:
            try:
                from multiprocessing import Pool
                pairs = list(itertools.combinations(range(n), 2))
                with Pool() as pool:
                    results = pool.starmap(sim_pair, pairs)
                for (i, j), val in zip(pairs, results):
                    self.similarity_matrix[i, j] = val
                    self.similarity_matrix[j, i] = val
                print("并行计算完成")
            except Exception as e:
                print(f"并行失败，回退串行: {e}")
                self._calc_similarity_serial(sim_pair, n)
        else:
            self._calc_similarity_serial(sim_pair, n)

        self.distance_matrix = 1 - self.similarity_matrix
        self.parameters['similarity'] = {'metric': metric, 'n_terms': n}
        return self

    def _calc_similarity_serial(self, sim_pair_func, n):
        for i in range(n):
            for j in range(i+1, n):
                sim = sim_pair_func(i, j)
                self.similarity_matrix[i, j] = sim
                self.similarity_matrix[j, i] = sim

    # ---------- 3. 降维 ----------
    def embed_coords(self, method='umap', random_state=42, **kwargs):
        print(f"降维方法: {method}")
        n = len(self.df)
        if method.lower() == 'mds':
            mds = MDS(n_components=2, dissimilarity='precomputed',
                      random_state=random_state, max_iter=3000, eps=1e-9)
            self.coords_2d = mds.fit_transform(self.distance_matrix)
        elif method.lower() == 'tsne':
            tsne = TSNE(n_components=2, metric='precomputed',
                        random_state=random_state, perplexity=min(30, n-1))
            self.coords_2d = tsne.fit_transform(self.distance_matrix)
        elif method.lower() == 'umap':
            reducer = umap.UMAP(n_components=2, metric='precomputed',
                                random_state=random_state, **kwargs)
            self.coords_2d = reducer.fit_transform(self.distance_matrix)
        else:
            raise ValueError("method 可选 mds, tsne, umap")
        self.parameters['embedding'] = {'method': method, 'kwargs': kwargs}
        return self

    # ---------- 4. 聚类 ----------
    def estimate_bandwidth_knn(self, k=8):
        nbrs = NearestNeighbors(n_neighbors=k).fit(self.coords_2d)
        distances, _ = nbrs.kneighbors(self.coords_2d)
        bandwidth = np.mean(distances[:, -1])
        return bandwidth

    def perform_clustering(self, cluster_method='meanshift', k_neighbors=8,
                           min_bandwidth=0.3, hdbscan_min_cluster_size=3,
                           auto_fallback=True):
        if self.coords_2d is None:
            raise ValueError("请先执行 embed_coords()")
        print(f"聚类方法: {cluster_method}")
        if cluster_method == 'meanshift':
            bandwidth = self.estimate_bandwidth_knn(k=k_neighbors)
            if bandwidth < min_bandwidth:
                bandwidth = min_bandwidth
            print(f"自适应带宽: {bandwidth:.4f}")
            ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=True)
            self.cluster_labels = ms.fit_predict(self.coords_2d)
            self.cluster_centers = ms.cluster_centers_
            unique = np.unique(self.cluster_labels)
            print(f"Mean-Shift 完成，{len(unique)} 个簇")
            if auto_fallback and len(unique) == 1 and len(self.df) > 5:
                print("警告：所有点为一簇，缩小带宽重试...")
                ms2 = MeanShift(bandwidth=bandwidth*0.6, bin_seeding=True)
                self.cluster_labels = ms2.fit_predict(self.coords_2d)
                self.cluster_centers = ms2.cluster_centers_
                print(f"重试后簇数: {len(np.unique(self.cluster_labels))}")
        elif cluster_method == 'hdbscan':
            clusterer = hdbscan.HDBSCAN(min_cluster_size=hdbscan_min_cluster_size,
                                        metric='euclidean')
            self.cluster_labels = clusterer.fit_predict(self.coords_2d)
            noise_mask = self.cluster_labels == -1
            if noise_mask.any():
                max_label = self.cluster_labels.max()
                self.cluster_labels[noise_mask] = max_label + 1
                print(f"将 {noise_mask.sum()} 个噪声点转为新簇 {max_label+1}")
            self.cluster_centers = None
            print(f"HDBSCAN 完成，{len(np.unique(self.cluster_labels))} 个簇")
        else:
            raise ValueError("cluster_method 须为 'meanshift' 或 'hdbscan'")
        self.parameters['clustering'] = {'method': cluster_method,
                                         'n_clusters': len(np.unique(self.cluster_labels))}
        return self

    # ---------- 5. 质量评估 ----------
    def evaluate_clustering(self):
        if self.cluster_labels is None:
            raise ValueError("请先执行聚类")
        if self.results_df is None:
            raise ValueError("请先执行 compile_results() 以创建 results_df")
        if len(np.unique(self.cluster_labels)) < 2:
            self.silhouette_avg = -1
            print("无法计算轮廓系数（簇数 < 2）")
            return self
        self.silhouette_avg = silhouette_score(self.distance_matrix, self.cluster_labels,
                                               metric='precomputed')
        print(f"整体轮廓系数: {self.silhouette_avg:.4f}")
        sil_vals = silhouette_samples(self.distance_matrix, self.cluster_labels,
                                      metric='precomputed')
        self.results_df['Silhouette'] = sil_vals
        self.results_df['Is_Boundary'] = sil_vals < 0.1
        self.parameters['quality'] = {'silhouette_score': self.silhouette_avg}
        return self

    # ---------- 6. 汇总结果 ----------
    def compile_results(self):
        self.results_df = self.df.copy()
        self.results_df['Cluster'] = self.cluster_labels
        self.results_df['X_coord'] = self.coords_2d[:, 0]
        self.results_df['Y_coord'] = self.coords_2d[:, 1]

        cluster_summary_list = []
        for cl in np.unique(self.cluster_labels):
            idx = self.results_df[self.results_df['Cluster'] == cl].index.tolist()
            size = len(idx)
            if size > 1:
                sub = self.similarity_matrix[np.ix_(idx, idx)]
                intra_sim = np.mean(sub[np.triu_indices_from(sub, k=1)])
            else:
                intra_sim = 1.0
            mean_p = self.results_df.loc[idx, 'Adjusted P-value'].mean()
            cluster_summary_list.append({
                'Cluster': cl,
                'Size': size,
                'Avg_Intra_Similarity': intra_sim,
                'Mean_Adjusted_P': mean_p
            })
        self.cluster_summary = pd.DataFrame(cluster_summary_list)
        return self

    # ---------- 7. 可视化 ----------
    def plot_cluster_scatter(self, save_path='cluster_scatter.png'):
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(self.coords_2d[:, 0], self.coords_2d[:, 1],
                   c=self.cluster_labels, cmap='tab20', s=100,
                   alpha=0.8, edgecolors='black')
        if self.cluster_centers is not None:
            ax.scatter(self.cluster_centers[:, 0], self.cluster_centers[:, 1],
                       marker='x', s=200, c='red', label='Centers')
        ax.set_xlabel('Dimension 1')
        ax.set_ylabel('Dimension 2')
        ax.set_title('GO term clustering (Mean‑Shift)')
        ax.legend()
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"散点图: {save_path}")

    def plot_similarity_heatmap(self, save_path='similarity_heatmap.png'):
        order = np.argsort(self.cluster_labels)
        sorted_mat = self.similarity_matrix[order][:, order]
        fig, ax = plt.subplots(figsize=(12, 10))
        im = ax.imshow(sorted_mat, cmap='YlOrRd', vmin=0, vmax=1)
        pos = 0
        for cl in np.unique(self.cluster_labels[order]):
            cnt = np.sum(self.cluster_labels[order] == cl)
            if cnt > 0:
                ax.add_patch(plt.Rectangle((pos-0.5, pos-0.5), cnt, cnt,
                                           linewidth=2, edgecolor='blue', facecolor='none'))
                pos += cnt
        plt.colorbar(im, label='Similarity')
        ax.set_title('Similarity matrix (sorted by cluster)')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"热图: {save_path}")

    def plot_bubble_chart(self, save_path='bubble_chart.png', x_cap=10):
        plot_df = self.results_df.copy()
        plot_df['neg_log_p'] = -np.log10(plot_df['Adjusted P-value'].clip(lower=1e-300))
        plot_df['neg_log_p_capped'] = plot_df['neg_log_p'].clip(upper=x_cap)
        plot_df['capped'] = plot_df['neg_log_p'] > x_cap
        plot_df['Short_Desc'] = plot_df['Description'].apply(lambda s: s[:60]+'...' if len(s)>63 else s)
        plot_df = plot_df.sort_values('Adjusted P-value')
        n = len(plot_df)
        fig_height = max(4, n*0.4+1)
        fig, ax = plt.subplots(figsize=(6, fig_height))
        clusters = np.unique(plot_df['Cluster'])
        colors = plt.cm.tab20(np.linspace(0,1,len(clusters)))
        color_dict = {cl: colors[i] for i,cl in enumerate(clusters)}
        for i, (_, row) in enumerate(plot_df.iterrows()):
            size = max(row['Count_all']*3, 5)
            ax.scatter(row['neg_log_p_capped'], i, s=size, alpha=0.7,
                       color=color_dict[row['Cluster']], edgecolors='black')
            if row['capped']:
                ax.annotate(f'({row["neg_log_p"]:.1f})', (row['neg_log_p_capped'], i),
                            xytext=(3,0), textcoords='offset points', fontsize=7)
        ax.set_yticks(range(n))
        ax.set_yticklabels(plot_df['Short_Desc'])
        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.set_title('Bubble chart')
        ax.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"气泡图: {save_path}")

    def plot_network_view(self, similarity_threshold=0.3, save_path='network_view.png'):
        fig, ax = plt.subplots(figsize=(14, 10))
        uniq = np.unique(self.cluster_labels)
        colors = plt.cm.tab20(np.linspace(0,1,len(uniq)))
        for cl, col in zip(uniq, colors):
            idx = self.cluster_labels == cl
            ax.scatter(self.coords_2d[idx,0], self.coords_2d[idx,1],
                       s=120, color=col, label=f'Cluster {cl}', alpha=0.8)
        n = len(self.df)
        for i in range(n):
            for j in range(i+1, n):
                if self.similarity_matrix[i,j] > similarity_threshold:
                    ax.plot([self.coords_2d[i,0], self.coords_2d[j,0]],
                            [self.coords_2d[i,1], self.coords_2d[j,1]],
                            'gray', alpha=0.3, linewidth=0.8)
        for i, (x,y) in enumerate(self.coords_2d):
            ax.annotate(self.results_df.iloc[i]['Description'][:40], (x,y),
                        fontsize=7, ha='center', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.6))
        ax.set_xlabel('Dimension 1')
        ax.set_ylabel('Dimension 2')
        ax.set_title(f'Network (sim > {similarity_threshold})')
        ax.legend(bbox_to_anchor=(1.05,1))
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"网络图: {save_path}")

    def plot_interactive_network(self, similarity_threshold=0.3, save_path='network_interactive.html'):
        net = Network(height="800px", width="100%", bgcolor="#ffffff")
        for idx, row in self.results_df.iterrows():
            title = f"GO:{row['GO']}<br>{row['Description']}<br>P={row['P-value']:.2e}<br>Cluster={row['Cluster']}"
            net.add_node(idx, label=row['Description'][:30], title=title,
                         color=self._get_cluster_color(row['Cluster']), size=20)
        for i in range(len(self.df)):
            for j in range(i + 1, len(self.df)):
                if self.similarity_matrix[i, j] > similarity_threshold:
                    net.add_edge(i, j, value=self.similarity_matrix[i, j],
                                 title=f"Sim: {self.similarity_matrix[i, j]:.2f}")
        # 使用合法的 JSON 字符串设置选项
        net.set_options('{"physics": {"enabled": true, "stabilization": {"iterations": 100}}}')
        net.save_graph(save_path)
        print(f"交互网络: {save_path}")

    def _get_cluster_color(self, cluster_id):
        palette = plt.cm.tab20(np.linspace(0,1,20))
        return '#%02x%02x%02x' % (int(palette[cluster_id%20][0]*255),
                                  int(palette[cluster_id%20][1]*255),
                                  int(palette[cluster_id%20][2]*255))

    def plot_wordclouds(self, output_dir='wordclouds'):
        os.makedirs(output_dir, exist_ok=True)
        stopwords = set(['the','of','and','to','in','for','on','by','with','is','are','a','an','be','this','that',
                         'activity','process','response','regulation','positive','negative','involved','mediated'])
        for cl in np.unique(self.cluster_labels):
            texts = self.results_df[self.results_df['Cluster'] == cl]['Description'].tolist()
            corpus = ' '.join(texts).lower()
            corpus = re.sub(r'[^\w\s]', ' ', corpus)
            wc = WordCloud(width=800, height=400, background_color='white',
                           stopwords=stopwords, max_words=50).generate(corpus)
            plt.figure(figsize=(10,5))
            plt.imshow(wc, interpolation='bilinear')
            plt.axis('off')
            plt.title(f'Cluster {cl} Word Cloud')
            plt.savefig(f"{output_dir}/cluster_{cl}_wordcloud.png", dpi=300)
            plt.close()
        print(f"词云保存至 {output_dir}")

    # ---------- 8. 输出结果 ----------
    def save_results(self, output_path='clustering_results.xlsx'):
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            self.results_df[['GO','Description','P-value','Adjusted P-value','Cluster','Silhouette','Is_Boundary',
                             'X_coord','Y_coord','Count_all','Merged_GOs','Gene_Set','All_index']].to_excel(
                writer, sheet_name='Clustering_Results', index=False)
            self.cluster_summary.to_excel(writer, sheet_name='Cluster_Summary', index=False)
            param_df = pd.DataFrame([(k,str(v)) for k,v in self.parameters.items()], columns=['Parameter','Value'])
            param_df.to_excel(writer, sheet_name='Parameters', index=False)
            if self.merged_terms_info:
                merged_df = pd.DataFrame(self.merged_terms_info)
                merged_df.to_excel(writer, sheet_name='Merged_Terms', index=False)
        print(f"结果保存: {output_path}")

    def export_cytoscape(self, similarity_threshold=0.3,
                         node_file='cytoscape_nodes.txt', edge_file='cytoscape_edges.txt'):
        nodes = self.results_df[['GO','Description','Cluster','Count_all','Adjusted P-value']].copy()
        nodes['neg_log10p'] = -np.log10(nodes['Adjusted P-value'].clip(lower=1e-300))
        nodes.to_csv(node_file, sep='\t', index=False)
        edges = []
        n = len(self.df)
        for i in range(n):
            for j in range(i+1, n):
                sim = self.similarity_matrix[i,j]
                if sim > similarity_threshold:
                    edges.append([self.df.iloc[i]['GO'], self.df.iloc[j]['GO'], sim])
        pd.DataFrame(edges, columns=['GO_A','GO_B','Similarity']).to_csv(edge_file, sep='\t', index=False)
        print(f"Cytoscape 文件: {node_file}, {edge_file}")


def run_clustering(input_file, sheet_name='Sheet1', output_prefix='go_clustering',
                   p_cutoff=0.05, min_genes=3, logp_median_ratio=True, keep_top_n=15,
                   merge_identical=True, similarity_metric='jaccard', parallel=True,
                   embed_method='umap', random_state=42, cluster_method='meanshift',
                   k_neighbors=8, min_bandwidth=0.3, hdbscan_min_cluster_size=3,
                   auto_fallback=True, similarity_threshold=0.3, no_plots=False,
                   wordcloud_dir='wordclouds', export_cytoscape=True):
    """对GO富集结果进行聚类分析，生成多种输出文件"""
    # 实例化 GOTermClusteringAnalyzer，按参数调用各方法
    # 示例框架：
    analyzer = GOTermClusteringAnalyzer(input_file, sheet_name=sheet_name)
    analyzer.load_data()
    analyzer.filter_terms_dynamic(p_cutoff=p_cutoff, min_genes=min_genes,
                                  logp_median_ratio=logp_median_ratio,
                                  keep_top_n=keep_top_n)
    if merge_identical:
        analyzer.merge_identical_sets()
    analyzer.calculate_similarity(metric=similarity_metric, parallel=parallel)
    analyzer.embed_coords(method=embed_method, random_state=random_state)
    analyzer.perform_clustering(cluster_method=cluster_method,
                                k_neighbors=k_neighbors,
                                min_bandwidth=min_bandwidth,
                                hdbscan_min_cluster_size=hdbscan_min_cluster_size,
                                auto_fallback=auto_fallback)
    analyzer.compile_results()
    analyzer.evaluate_clustering()
    if not no_plots:
        analyzer.plot_cluster_scatter(f"{output_prefix}_scatter.png")
        analyzer.plot_similarity_heatmap(f"{output_prefix}_heatmap.png")
        analyzer.plot_bubble_chart(f"{output_prefix}_bubble.png")
        analyzer.plot_network_view(similarity_threshold=similarity_threshold,
                                   save_path=f"{output_prefix}_network.png")
        analyzer.plot_interactive_network(similarity_threshold=similarity_threshold,
                                          save_path=f"{output_prefix}_network.html")
        analyzer.plot_wordclouds(output_dir=f"{output_prefix}_{wordcloud_dir}")
    analyzer.save_results(f"{output_prefix}_results.xlsx")
    if export_cytoscape:
        analyzer.export_cytoscape(similarity_threshold=similarity_threshold,
                                  node_file=f"{output_prefix}_cytoscape_nodes.txt",
                                  edge_file=f"{output_prefix}_cytoscape_edges.txt")