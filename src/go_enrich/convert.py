#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO annotation file format conversion module (obtain descriptions from obo file)
Function: Automatically identify the gene ID column in any Excel file, extract all GO IDs,
          and then retrieve the corresponding names from the GO basic ontology file (.obo)
          as descriptions. Output a standard two-column format (ID, Annotation).
Usage (command line):
    go-enrich convert input.xlsx output.xlsx [--obo go-basic.obo]
Usage (Python):
    from go_enrich.convert import convert_any_go_annotation
    convert_any_go_annotation('input.xlsx', 'output.xlsx', obo_path='go-basic.obo')
"""

import pandas as pd
import re
from collections import defaultdict
import os


def detect_id_column(df):
    """
    Automatically identify the ID column.
    Prefer columns whose names contain 'id' or 'gene' (case-insensitive);
    if none found, default to the first column.
    """
    lower_cols = {col: col.lower() for col in df.columns}
    for col, lc in lower_cols.items():
        if 'id' in lc or 'gene' in lc:
            return col
    return df.columns[0]


def load_go_names(obo_path):
    """
    Parse the GO basic ontology file (.obo) and return a dict {GO_id: name}.
    """
    if not os.path.exists(obo_path):
        raise FileNotFoundError(f"GO ontology file not found: {obo_path}")

    go_names = {}
    current_id = None
    current_name = None
    with open(obo_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[Term]'):
                # New Term starts, reset
                current_id = None
                current_name = None
            elif line.startswith('id: GO:'):
                current_id = line.split('id: ')[1].strip()
            elif line.startswith('name: '):
                current_name = line.split('name: ')[1].strip()
            elif line == '' and current_id and current_name:
                # Term ends, save
                go_names[current_id] = current_name
                current_id = None
                current_name = None
    # Handle the last Term if no trailing blank line
    if current_id and current_name:
        go_names[current_id] = current_name

    print(f"Loaded {len(go_names)} GO terms from {obo_path}")
    return go_names


def extract_go_ids_from_text(text):
    """
    Extract all GO IDs (GO:xxxxxx format) from arbitrary text, return list of GO IDs.
    """
    if pd.isna(text):
        return []
    text = str(text)
    # Match GO: followed by digits (optional dot and version)
    return re.findall(r'GO:\d+(?:\.\d+)?', text)


def extract_all_go_entries(df, id_col):
    """
    Traverse all columns except the ID column, extract all GO IDs, aggregate by gene ID.
    Return dict: {gene_id: set of go_ids}
    """
    gene_go_map = defaultdict(set)
    other_cols = [col for col in df.columns if col != id_col]
    for _, row in df.iterrows():
        gene_id = str(row[id_col]).strip()
        if not gene_id:
            continue
        combined_text = ''
        for col in other_cols:
            val = row[col]
            if pd.notna(val):
                combined_text += str(val) + '; '
        if combined_text:
            go_ids = extract_go_ids_from_text(combined_text)
            for go in go_ids:
                if go:
                    gene_go_map[gene_id].add(go)
    return gene_go_map


def convert_any_go_annotation(input_path, output_path, obo_path='go-basic.obo'):
    """
    Main conversion function: read arbitrary Excel file, extract GO IDs,
    use obo file to get descriptions, output standard format.
    """
    # Read input file
    try:
        df = pd.read_excel(input_path, engine='openpyxl')
    except Exception:
        df = pd.read_excel(input_path)  # fallback default engine
    print(f"Reading file: {input_path}, {len(df)} rows, columns: {list(df.columns)}")

    # Identify ID column
    id_col = detect_id_column(df)
    print(f"Identified ID column: '{id_col}'")

    # Load obo mapping
    go_names = load_go_names(obo_path)

    # Extract all GO IDs
    gene_go_map = extract_all_go_entries(df, id_col)
    print(f"Extracted GO annotations for {len(gene_go_map)} genes")

    rows = []
    for gene, go_set in gene_go_map.items():
        sorted_go = sorted(go_set)  # sort by GO ID
        annotation_parts = []
        for go in sorted_go:
            desc = go_names.get(go, '')  # empty if not found
            if desc:
                annotation_parts.append(f"{go},{desc}")
            else:
                annotation_parts.append(go)  # only GO ID if no description
        annotation = '; '.join(annotation_parts)
        rows.append({'ID': gene, 'Annotation': annotation})

    result_df = pd.DataFrame(rows)
    result_df.to_excel(output_path, index=False)
    print(f"Conversion complete! Saved to {output_path}, {len(result_df)} genes")


# If run directly, act as a standalone tool
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python convert.py input.xlsx output.xlsx [--obo go-basic.obo]")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    obo_path = "go-basic.obo"
    if len(sys.argv) >= 4 and sys.argv[3].startswith('--obo'):
        obo_path = sys.argv[4] if len(sys.argv) > 4 else "go-basic.obo"
    convert_any_go_annotation(input_file, output_file, obo_path)