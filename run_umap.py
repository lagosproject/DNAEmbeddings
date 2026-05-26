import numpy as np
import pandas as pd
import umap
import json
import glob
import os
import re

def natural_key(string):
    """Key function for natural sorting of alphanumeric strings (e.g. chr1, chr2, chr10)."""
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', string)]

def main():
    # Dynamically find all embedding files in the res/ directory
    emb_files = glob.glob("res/*_embeddings.npy")
    
    chromosomes = []
    for f in emb_files:
        basename = os.path.basename(f)
        match = re.match(r"^([a-zA-Z0-9]+)_embeddings\.npy$", basename)
        if match:
            chromosomes.append(match.group(1))
            
    chromosomes = sorted(chromosomes, key=natural_key)
    
    print(f"Dynamically discovered chromosomes: {chromosomes}")
    
    all_embeddings = []
    all_metadata = []

    print("Loading embeddings and metadata...")
    for chrom in chromosomes:
        emb_path = f"res/{chrom}_embeddings.npy"
        meta_path = f"res/{chrom}_metadata.csv"
        
        if os.path.exists(emb_path) and os.path.exists(meta_path):
            emb = np.load(emb_path)
            meta = pd.read_csv(meta_path)
            
            # Verify shapes
            if len(emb) != len(meta):
                print(f"Warning: Size mismatch for {chrom}. Embeddings: {len(emb)}, Metadata: {len(meta)}. Aligning...")
                min_len = min(len(emb), len(meta))
                emb = emb[:min_len]
                meta = meta.iloc[:min_len]
            
            all_embeddings.append(emb)
            all_metadata.append(meta)
            print(f"Loaded {chrom}: {len(emb)} bins.")
        else:
            print(f"Could not find matching metadata file for {chrom}")

    if not all_embeddings:
        print("No embedding files found!")
        return

    # Concatenate all embeddings and metadata
    X = np.vstack(all_embeddings)
    df = pd.concat(all_metadata, ignore_index=True)
    
    print(f"\nCombined shape: {X.shape}")

    # 1. Run Global UMAP (Genome-wide)
    print("\nRunning global UMAP (Genome-wide)...")
    reducer_global = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        metric='cosine',
        random_state=42,
        verbose=True
    )
    global_coords = reducer_global.fit_transform(X)
    df['global_x'] = global_coords[:, 0].astype(float)
    df['global_y'] = global_coords[:, 1].astype(float)

    # 2. Run Local UMAP (Chromosome-specific)
    print("\nRunning local UMAPs per chromosome...")
    local_x = np.zeros(len(df))
    local_y = np.zeros(len(df))

    for chrom in chromosomes:
        idx = df[df['chrom'] == chrom].index
        if len(idx) > 0:
            print(f"Processing UMAP for {chrom}...")
            X_chrom = X[idx]
            
            n_neighbors = min(15, len(X_chrom) - 1)
            if n_neighbors < 2:
                coords = np.zeros((len(X_chrom), 2))
            else:
                reducer_local = umap.UMAP(
                    n_neighbors=n_neighbors,
                    min_dist=0.1,
                    metric='cosine',
                    random_state=42,
                    verbose=False
                )
                coords = reducer_local.fit_transform(X_chrom)
            
            local_x[idx] = coords[:, 0]
            local_y[idx] = coords[:, 1]

    df['local_x'] = local_x
    df['local_y'] = local_y

    # Save to JSON
    output_path = "res/umap_results.json"
    print(f"\nSaving results to {output_path}...")
    result_dict = df.to_dict(orient='list')
    with open(output_path, 'w') as f:
        json.dump(result_dict, f)
    
    print("Done! UMAP calculations completed successfully.")

if __name__ == "__main__":
    main()
