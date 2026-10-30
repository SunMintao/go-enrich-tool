#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download GO annotation (GAF) from NCBI RefSeq and convert it to the standard format.

Supported source: NCBI RefSeq plant genomes
    https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/{species}/
      annotation_releases/{assembly}/{assembly}_gene_ontology.gaf.gz

The converted standard file uses NCBI gene IDs (the GAF "DB Object ID" column),
so it only supports gene IDs in the NCBI RefSeq namespace.
"""

import os
import re
import gzip
from urllib.parse import unquote

import pandas as pd
import requests

from .utils import get_default_obo_path, parse_obo_file

NCBI_REFSEQ_PLANT_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/"


def _list_remote_dir(url):
    """List a remote NCBI FTP directory over HTTPS; return (dirs, files)."""
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    dirs, files = [], []
    for m in re.finditer(r'href="([^"]+)"', resp.text):
        href = m.group(1)
        if href in ('../', './') or href.startswith('?'):
            continue
        name = unquote(href).rstrip('/')
        if not name:
            continue
        if href.endswith('/'):
            if name not in dirs:
                dirs.append(name)
        else:
            if name not in files:
                files.append(name)
    return dirs, files


def list_species():
    """Return the available plant species names (Latin names) on NCBI RefSeq."""
    dirs, _ = _list_remote_dir(NCBI_REFSEQ_PLANT_BASE)
    return sorted(dirs)


def resolve_gaf_url(species, assembly=None):
    """Resolve the GAF file URL for a species. Return (gaf_url, assembly_dir)."""
    base = NCBI_REFSEQ_PLANT_BASE.rstrip('/') + '/' + species + '/'
    dirs, _ = _list_remote_dir(base)
    if 'annotation_releases' in dirs:
        release_dir = 'annotation_releases'
    elif 'all_assembly_versions' in dirs:
        release_dir = 'all_assembly_versions'
    else:
        raise ValueError(f"No annotation directory found for species '{species}' under {base}")

    rel_url = base + release_dir + '/'
    asm_dirs, _ = _list_remote_dir(rel_url)
    gcf_dirs = sorted(d for d in asm_dirs if d.startswith('GCF_'))
    if not gcf_dirs:
        raise ValueError(f"No GCF_ assembly found under {rel_url}")

    if assembly:
        if assembly not in asm_dirs:
            raise ValueError(f"Assembly '{assembly}' not found under {rel_url}; available: {asm_dirs}")
        asm_dir = assembly
    else:
        asm_dir = gcf_dirs[-1]  # latest assembly by name

    _, files = _list_remote_dir(rel_url + asm_dir + '/')
    gaf_files = sorted(f for f in files if f.endswith('_gene_ontology.gaf.gz'))
    if not gaf_files:
        raise ValueError(f"No *_gene_ontology.gaf.gz found under {rel_url}{asm_dir}/")
    gaf_url = rel_url + asm_dir + '/' + gaf_files[0]
    return gaf_url, asm_dir


def download_gaf(gaf_url, output_dir='.'):
    """Download a (gzipped) GAF file into output_dir; return the local path."""
    os.makedirs(output_dir, exist_ok=True)
    fname = gaf_url.rstrip('/').split('/')[-1]
    local = os.path.join(output_dir, fname)
    resp = requests.get(gaf_url, stream=True, timeout=600)
    resp.raise_for_status()
    with open(local, 'wb') as f:
        for chunk in resp.iter_content(1024 * 1024):
            if chunk:
                f.write(chunk)
    return local


def parse_gaf(gaf_path):
    """Parse a GAF 2.x file into {gene_id: set(go_ids)}, skipping NOT annotations."""
    gene2go = {}
    opener = gzip.open if str(gaf_path).endswith('.gz') else open
    with opener(gaf_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                continue
            gene_id = cols[1].strip()
            go_id = cols[4].strip()
            qualifiers = [q.strip() for q in cols[3].split('|')] if cols[3].strip() else []
            if not gene_id or not go_id.startswith('GO:'):
                continue
            if 'NOT' in qualifiers:
                continue
            gene2go.setdefault(gene_id, set()).add(go_id)
    return gene2go


def gaf_to_standard(gaf_path, output_path, obo_path=None):
    """Convert a GAF file to the standard two-column format (ID, Annotation)."""
    if obo_path and os.path.exists(obo_path):
        pass
    elif obo_path == 'go-basic.obo' or not obo_path:
        default = get_default_obo_path()
        if default:
            obo_path = default
            print(f"Using built-in obo file: {obo_path}")
    go_names = parse_obo_file(obo_path)
    gene2go = parse_gaf(gaf_path)

    rows = []
    for gene in sorted(gene2go):
        parts = []
        for go in sorted(gene2go[gene]):
            desc = go_names.get(go, '')
            parts.append(f"{go},{desc}" if desc else go)
        rows.append({'ID': gene, 'Annotation': '; '.join(parts)})
    df = pd.DataFrame(rows)
    df.to_excel(output_path, index=False, engine='openpyxl')
    print(f"Converted {len(df)} genes from {gaf_path} to {output_path}")
    return df


def run_download_go(species, output_file, obo_path=None, assembly=None, keep_gaf=False):
    """Download GO annotation for a species and convert it to the standard format."""
    print(f"Resolving GO annotation for species '{species}' from NCBI RefSeq ...")
    gaf_url, asm_dir = resolve_gaf_url(species, assembly=assembly)
    print(f"Found GAF file: {gaf_url}")
    local_gaf = download_gaf(gaf_url)
    print(f"Downloaded: {local_gaf}")
    gaf_to_standard(local_gaf, output_file, obo_path=obo_path)
    if not keep_gaf:
        os.remove(local_gaf)
        print(f"Removed temporary GAF file: {local_gaf}")
    return output_file
