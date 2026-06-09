import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import requests
import os
import re
from collections import defaultdict

# --------------------------- 1. 读取 GO 注释 ---------------------------
def read_go_annotation(file_path):
    """
    读取 go_annotation.xlsx，返回:
    gene2gos: dict {gene: set of GO ids}
    all_genes: list of all genes (背景基因集)
    go_descriptions: dict {GO_id: description}
    """
    df = pd.read_excel(file_path)
    # 假设列名是 'ID' 和 'Annotation'
    gene2gos = {}
    all_genes = []
    go_descriptions = {}
    for idx, row in df.iterrows():
        gene = str(row['ID']).strip()
        ann = str(row['Annotation'])
        if pd.isna(ann) or ann == 'nan':
            continue
        # 按分号分割每个 GO 条目
        go_terms = ann.split(';')
        go_ids = set()
        for term in go_terms:
            term = term.strip()
            if ',' in term:
                go_id, desc = term.split(',', 1)
                go_id = go_id.strip()
                desc = desc.strip()
                go_ids.add(go_id)
                if go_id not in go_descriptions:
                    go_descriptions[go_id] = desc
            else:
                # 如果没有逗号，可能只是 GO ID？跳过
                continue
        gene2gos[gene] = go_ids
        all_genes.append(gene)
    return gene2gos, all_genes, go_descriptions

def read_gene_list(file_path):
    """读取 id.xlsx，返回目标基因列表（第一列）"""
    df = pd.read_excel(file_path, skiprows=1, header=None)
    genes = df[0].dropna().astype(str).str.strip().tolist()
    return genes

# --------------------------- 2. 获取 GO 分类（namespace） ---------------------------
def download_obo(url="http://purl.obolibrary.org/obo/go/go-basic.obo", local_path="go-basic.obo"):
    """如果本地没有 obo 文件，则从网络下载"""
    if not os.path.exists(local_path):
        print("Downloading go-basic.obo ...")
        r = requests.get(url)
        with open(local_path, 'wb') as f:
            f.write(r.content)
    return local_path

def parse_obo_namespace(obo_path):
    """解析 obo 文件，返回 dict {GO_id: namespace} 其中 namespace 为 'biological_process', 'cellular_component', 'molecular_function'"""
    go_namespace = {}
    with open(obo_path, 'r') as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split(": ")[1]
            elif line.startswith("namespace: ") and current_id:
                ns = line.split(": ")[1]
                go_namespace[current_id] = ns
                current_id = None  # 重置，避免下一个 term 误用
    return go_namespace

# --------------------------- 3. 富集分析 ---------------------------
def enrichment_analysis(gene2gos, go_descriptions, all_genes, target_genes, go_namespace=None, ontology='BP'):
    """
    对指定 ontology (BP/CC/MF) 或全部 GO term 进行富集分析
    返回 DataFrame
    """
    # 背景基因总数
    N = len(all_genes)
    n = len(target_genes)
    # 目标基因集合
    target_set = set(target_genes)

    # 收集所有出现的 GO term
    all_go_terms = set()
    for gos in gene2gos.values():
        all_go_terms.update(gos)

    # 过滤出指定 ont 的 GO term（如果提供了 namespace 且 ontology 不是 'ALL'）
    if go_namespace and ontology != 'ALL':
        # 映射字符串简化
        onto_map = {'BP': 'biological_process', 'CC': 'cellular_component', 'MF': 'molecular_function'}
        wanted_ns = onto_map.get(ontology, ontology)
        go_terms = [go for go in all_go_terms if go_namespace.get(go) == wanted_ns]
    else:
        go_terms = list(all_go_terms)

    # 统计背景中每个 GO term 的基因数 (K)
    bg_counts = defaultdict(int)
    for gene, gos in gene2gos.items():
        for go in gos:
            if go in go_terms:
                bg_counts[go] += 1

    # 统计目标中每个 GO term 的基因数 (k) 和具体的基因列表
    target_counts = defaultdict(int)
    target_genes_list = defaultdict(list)
    for gene in target_set:
        if gene in gene2gos:
            for go in gene2gos[gene]:
                if go in go_terms:
                    target_counts[go] += 1
                    target_genes_list[go].append(gene)

    # 准备结果列表
    results = []
    for go in go_terms:
        if go not in bg_counts or bg_counts[go] == 0:
            continue
        K = bg_counts[go]
        k = target_counts.get(go, 0)
        if k == 0:
            # 可以跳过 p-value 非常大的项，但保留以便输出所有？通常只输出有差异的，但和示例一致，可以保留
            pass
        # 超几何检验：从 N 个基因中抽取 n 个，其中该 GO 背景有 K 个，观测到 k 个的概率 (双侧？通常用右尾)
        # 使用 hypergeom.sf(k-1, N, K, n) 得到 P(X >= k)
        pval = hypergeom.sf(k-1, N, K, n) if k > 0 else 1.0
        # 计算比例
        diff_ratio = f"{k}/{n} {k/n*100:.2f}%" if n>0 else "0/0 0%"
        bg_ratio = f"{K}/{N} {K/N*100:.2f}%" if N>0 else "0/0 0%"
        all_index = ';'.join(target_genes_list.get(go, []))
        results.append({
            'GO': go,
            'Description': go_descriptions.get(go, ''),
            'DiffRatio': diff_ratio,
            'BgRatio': bg_ratio,
            'P-value': pval,
            'All_index': all_index,
            'Count_all': k
        })

    df_res = pd.DataFrame(results)
    if df_res.empty:
        return df_res

    # 多重假设校正 (Benjamini-Hochberg)
    _, pvals_adj, _, _ = multipletests(df_res['P-value'].values, method='fdr_bh')
    df_res['Adjusted P-value'] = pvals_adj

    # 按 p-value 排序
    df_res = df_res.sort_values('P-value').reset_index(drop=True)

    # 添加 GO_level_1 列: 根据 namespace 填入 (如果存在)
    if go_namespace:
        df_res['GO_level_1'] = df_res['GO'].apply(lambda x: go_namespace.get(x, ''))
        # 转换为可读格式
        df_res['GO_level_1'] = df_res['GO_level_1'].map({
            'biological_process': 'Biological process',
            'cellular_component': 'Cellular component',
            'molecular_function': 'Molecular function'
        }).fillna('')
    else:
        df_res['GO_level_1'] = ''

    # 调整列顺序以匹配示例
    cols = ['GO_level_1', 'GO', 'Description', 'DiffRatio', 'BgRatio', 'P-value', 'Adjusted P-value', 'All_index', 'Count_all']
    df_res = df_res[cols]
    return df_res

# --------------------------- 主程序 ---------------------------
def run_enrichment(go_anno_file, gene_list_file, output_file, ontology='BP', use_obo=True, obo_path='go-basic.obo'):
    # 读取数据
    gene2gos, all_genes, go_descriptions = read_go_annotation(go_anno_file)
    target_genes = read_gene_list(gene_list_file)

    # 确保目标基因都在背景中
    missing = set(target_genes) - set(all_genes)
    if missing:
        print(f"Warning: {len(missing)} target genes not found in annotation file: {missing}")

    # 获取 GO 分类
    go_namespace = None
    if use_obo:
        try:
            obo_path = download_obo(local_path=obo_path)
            go_namespace = parse_obo_namespace(obo_path)
        except Exception as e:
            print(f"Failed to load obo file: {e}. Proceeding without GO classification.")
            go_namespace = None

    # 执行富集分析
    df_res = enrichment_analysis(gene2gos, go_descriptions, all_genes, target_genes, go_namespace, ontology)

    # 输出 Excel
    if not df_res.empty:
        df_res.to_excel(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("No significant GO terms found.")