import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import requests
import os
import re
from collections import defaultdict

# --------------------------- 1. Read GO annotation ---------------------------
def read_go_annotation(file_path):
    """
    Read go_annotation.xlsx and return:
    gene2gos: dict {gene: set of GO ids}
    all_genes: list of all genes (background gene set)
    go_descriptions: dict {GO_id: description}
    """
    df = pd.read_excel(file_path)
    # Assume columns are 'ID' and 'Annotation'
    gene2gos = {}
    all_genes = []
    go_descriptions = {}
    for idx, row in df.iterrows():
        gene = str(row['ID']).strip()
        ann = str(row['Annotation'])
        if pd.isna(ann) or ann == 'nan':
            continue
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
                # If no comma, maybe only GO ID? skip
                continue
        gene2gos[gene] = go_ids
        all_genes.append(gene)
    return gene2gos, all_genes, go_descriptions

def read_gene_list(file_path):
    """Read id.xlsx and return target gene list (first column)"""
    df = pd.read_excel(file_path, skiprows=1, header=None)
    genes = df[0].dropna().astype(str).str.strip().tolist()
    return genes

# --------------------------- 2. Obtain GO classification (namespace) ---------------------------
def download_obo(url="http://purl.obolibrary.org/obo/go/go-basic.obo", local_path="go-basic.obo"):
    """Download obo file from network if not present locally"""
    if not os.path.exists(local_path):
        print("Downloading go-basic.obo ...")
        r = requests.get(url)
        with open(local_path, 'wb') as f:
            f.write(r.content)
    return local_path

def parse_obo_namespace(obo_path):
    """
    Parse obo file and return dict {GO_id: namespace}
    where namespace is 'biological_process', 'cellular_component', or 'molecular_function'
    """
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
                current_id = None  # reset to avoid misuse for next term
    return go_namespace

# --------------------------- 3. Enrichment analysis ---------------------------
def enrichment_analysis(gene2gos, go_descriptions, all_genes, target_genes, go_namespace=None, ontology='BP'):
    """
    Perform enrichment analysis for specified ontology (BP/CC/MF) or all GO terms.
    Return DataFrame.
    """
    N = len(all_genes)
    n = len(target_genes)
    target_set = set(target_genes)

    # Collect all GO terms appearing
    all_go_terms = set()
    for gos in gene2gos.values():
        all_go_terms.update(gos)

    # Filter by ontology if namespace is provided and ontology is not 'ALL'
    if go_namespace and ontology != 'ALL':
        onto_map = {'BP': 'biological_process', 'CC': 'cellular_component', 'MF': 'molecular_function'}
        wanted_ns = onto_map.get(ontology, ontology)
        go_terms = [go for go in all_go_terms if go_namespace.get(go) == wanted_ns]
    else:
        go_terms = list(all_go_terms)

    # Count background genes for each GO term (K)
    bg_counts = defaultdict(int)
    for gene, gos in gene2gos.items():
        for go in gos:
            if go in go_terms:
                bg_counts[go] += 1

    # Count target genes for each GO term (k) and list of genes
    target_counts = defaultdict(int)
    target_genes_list = defaultdict(list)
    for gene in target_set:
        if gene in gene2gos:
            for go in gene2gos[gene]:
                if go in go_terms:
                    target_counts[go] += 1
                    target_genes_list[go].append(gene)

    # Prepare results list
    results = []
    for go in go_terms:
        if go not in bg_counts or bg_counts[go] == 0:
            continue
        K = bg_counts[go]
        k = target_counts.get(go, 0)
        # p-value using hypergeometric test (right-tailed)
        pval = hypergeom.sf(k-1, N, K, n) if k > 0 else 1.0
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

    # Multiple testing correction (Benjamini-Hochberg)
    _, pvals_adj, _, _ = multipletests(df_res['P-value'].values, method='fdr_bh')
    df_res['Adjusted P-value'] = pvals_adj

    # Sort by p-value
    df_res = df_res.sort_values('P-value').reset_index(drop=True)

    # Add GO_level_1 column based on namespace if available
    if go_namespace:
        df_res['GO_level_1'] = df_res['GO'].apply(lambda x: go_namespace.get(x, ''))
        df_res['GO_level_1'] = df_res['GO_level_1'].map({
            'biological_process': 'Biological process',
            'cellular_component': 'Cellular component',
            'molecular_function': 'Molecular function'
        }).fillna('')
    else:
        df_res['GO_level_1'] = ''

    # Reorder columns to match expected format
    cols = ['GO_level_1', 'GO', 'Description', 'DiffRatio', 'BgRatio', 'P-value', 'Adjusted P-value', 'All_index', 'Count_all']
    df_res = df_res[cols]
    return df_res

# --------------------------- Main program ---------------------------
def run_enrichment(go_anno_file, gene_list_file, output_file, ontology='BP', use_obo=True, obo_path='go-basic.obo'):
    # Read data
    gene2gos, all_genes, go_descriptions = read_go_annotation(go_anno_file)
    target_genes = read_gene_list(gene_list_file)

    # Ensure target genes are in background
    missing = set(target_genes) - set(all_genes)
    if missing:
        print(f"Warning: {len(missing)} target genes not found in annotation file: {missing}")

    # Obtain GO classification
    go_namespace = None
    if use_obo:
        try:
            obo_path = download_obo(local_path=obo_path)
            go_namespace = parse_obo_namespace(obo_path)
        except Exception as e:
            print(f"Failed to load obo file: {e}. Proceeding without GO classification.")
            go_namespace = None

    # Perform enrichment analysis
    df_res = enrichment_analysis(gene2gos, go_descriptions, all_genes, target_genes, go_namespace, ontology)

    # Output Excel
    if not df_res.empty:
        df_res.to_excel(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("No significant GO terms found.")