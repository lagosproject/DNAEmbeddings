#!/usr/bin/env python3
import numpy as np
import pandas as pd
import umap
import json
import glob
import os
import re

# Ensure we work from project root (one level up from scripts/)
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

def natural_key(string):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', string)]

def main():
    print("=" * 60)
    print("ALPHAGenome UMAP Explorer Extension: Gene & Pathway Projections")
    print("=" * 60)

    # 1. Load genomic bin embeddings and metadata
    emb_files = glob.glob("res/*_embeddings.npy")
    chromosomes = []
    for f in emb_files:
        basename = os.path.basename(f)
        match = re.match(r"^([a-zA-Z0-9]+)_embeddings\.npy$", basename)
        if match:
            chromosomes.append(match.group(1))
    chromosomes = sorted(chromosomes, key=natural_key)
    
    print(f"Loading embeddings for chromosomes: {chromosomes}...")
    all_embeddings = []
    all_metadata = []
    for chrom in chromosomes:
        emb_path = f"res/{chrom}_embeddings.npy"
        meta_path = f"res/{chrom}_metadata.csv"
        if os.path.exists(emb_path) and os.path.exists(meta_path):
            emb = np.load(emb_path)
            meta = pd.read_csv(meta_path)
            if len(emb) != len(meta):
                min_len = min(len(emb), len(meta))
                emb = emb[:min_len]
                meta = meta.iloc[:min_len]
            all_embeddings.append(emb)
            all_metadata.append(meta)
    
    if not all_embeddings:
        print("Error: No embedding files found in res/")
        return
        
    X = np.vstack(all_embeddings)
    df = pd.concat(all_metadata, ignore_index=True)
    print(f"Loaded {len(X)} genomic bins with embedding dimension {X.shape[1]}.")

    # 2. Load annotations
    annotations_path = "res/gene_annotations.json"
    if not os.path.exists(annotations_path):
        print(f"Error: {annotations_path} not found! Please run download_annotations.py first.")
        return
        
    with open(annotations_path, 'r') as f:
        annotations = json.load(f)
    bin_genes = annotations["bin_genes"]
    gene_data = annotations["gene_data"]

    # 3. Build Gene-level Embeddings
    print("\nProcessing Gene-level Embeddings...")
    gene_to_bins = {}
    for i, genes in enumerate(bin_genes):
        for g in genes:
            gene_to_bins.setdefault(g, []).append(i)
            
    gene_symbols = sorted(list(gene_to_bins.keys()))
    print(f"Found {len(gene_symbols)} unique genes with at least one genomic bin mapping.")
    
    gene_embeddings = []
    gene_chroms = []
    gene_starts = []
    gene_ends = []
    
    for gene in gene_symbols:
        bin_indices = gene_to_bins[gene]
        # Average the embeddings of the bins this gene overlaps with
        emb_subset = X[bin_indices]
        gene_embeddings.append(np.mean(emb_subset, axis=0))
        
        # Get genomic location coordinates based on first bin (or overall span)
        chroms = [df.iloc[idx]['chrom'] for idx in bin_indices]
        starts = [int(df.iloc[idx]['start']) for idx in bin_indices]
        ends = [int(df.iloc[idx]['end']) for idx in bin_indices]
        
        # Take primary chromosome (most common)
        primary_chrom = max(set(chroms), key=chroms.count)
        gene_chroms.append(primary_chrom)
        gene_starts.append(min(starts))
        gene_ends.append(max(ends))
        
    gene_X = np.array(gene_embeddings)
    print(f"Gene embedding matrix shape: {gene_X.shape}")

    # 4. Run Gene-level UMAP
    print("Running UMAP on Gene Embeddings (this may take a minute)...")
    reducer_gene = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        metric='cosine',
        random_state=42,
        verbose=True
    )
    gene_coords = reducer_gene.fit_transform(gene_X)
    
    gene_results = {
        "symbol": gene_symbols,
        "x": gene_coords[:, 0].astype(float).tolist(),
        "y": gene_coords[:, 1].astype(float).tolist(),
        "chrom": gene_chroms,
        "start": gene_starts,
        "end": gene_ends
    }
    
    gene_output_path = "res/gene_umap_results.json"
    print(f"Saving Gene UMAP results to {gene_output_path}...")
    with open(gene_output_path, 'w') as f:
        json.dump(gene_results, f)

    # 5. Build Pathway-level Embeddings
    print("\nProcessing Pathway-level Embeddings...")
    pathway_to_genes = {}
    pathway_meta = {} # Maps pathway ID to (name, source)
    
    # Map gene symbol back to index in gene_symbols list to retrieve its gene embedding
    gene_symbol_to_idx = {sym: idx for idx, sym in enumerate(gene_symbols)}
    
    for g, info in gene_data.items():
        g_upper = g.upper()
        if g_upper not in gene_symbol_to_idx:
            continue
        for p in info.get("pathways", []):
            p_id = p["id"]
            pathway_to_genes.setdefault(p_id, []).append(g_upper)
            pathway_meta[p_id] = (p["name"], p["source"])
            
    pathway_ids = sorted(list(pathway_to_genes.keys()))
    print(f"Found {len(pathway_ids)} unique pathways with at least one gene mapping.")
    
    pathway_embeddings = []
    pathway_names = []
    pathway_sources = []
    
    for p_id in pathway_ids:
        mapped_genes = pathway_to_genes[p_id]
        gene_indices = [gene_symbol_to_idx[g] for g in mapped_genes]
        # Average the embeddings of the genes belonging to this pathway
        pathway_emb_subset = gene_X[gene_indices]
        pathway_embeddings.append(np.mean(pathway_emb_subset, axis=0))
        
        name, source = pathway_meta[p_id]
        pathway_names.append(name)
        pathway_sources.append(source)
        
    pathway_X = np.array(pathway_embeddings)
    print(f"Pathway embedding matrix shape: {pathway_X.shape}")

    # 6. Run Pathway-level UMAP
    print("Running UMAP on Pathway Embeddings...")
    reducer_pathway = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        metric='cosine',
        random_state=42,
        verbose=True
    )
    pathway_coords = reducer_pathway.fit_transform(pathway_X)
    
    pathway_results = {
        "id": pathway_ids,
        "name": pathway_names,
        "source": pathway_sources,
        "x": pathway_coords[:, 0].astype(float).tolist(),
        "y": pathway_coords[:, 1].astype(float).tolist()
    }
    
    pathway_output_path = "res/pathway_umap_results.json"
    print(f"Saving Pathway UMAP results to {pathway_output_path}...")
    with open(pathway_output_path, 'w') as f:
        json.dump(pathway_results, f)

    print("\nSUCCESS: Calculated both Gene-level and Pathway-level UMAP coordinates!")
    print("=" * 60)

if __name__ == "__main__":
    main()
