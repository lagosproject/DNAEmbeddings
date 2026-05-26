#!/usr/bin/env python3
import os
import json
import urllib.request
import pandas as pd
import requests

# Ensure we work from project root (one level up from scripts/)
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

def main():
    print("=" * 60)
    print("ALPHAGenome UMAP Explorer Annotation Downloader")
    print("=" * 60)
    
    # 1. Load UMAP results to get genomic bins
    umap_results_path = "res/umap_results.json"
    output_path = "res/gene_annotations.json"
    if not os.path.exists(umap_results_path):
        print(f"Error: {umap_results_path} not found! Please run run_umap.py first.")
        return
        
    print("Loading UMAP results...")
    with open(umap_results_path, 'r') as f:
        umap_data = json.load(f)
        
    chrom_list = umap_data['chrom']
    start_list = umap_data['start']
    end_list = umap_data['end']
    n_bins = len(chrom_list)
    print(f"Loaded {n_bins} genomic bins.")
    
    # 2. Load or Download RefSeq Gene annotations from UCSC
    refgene_local_path = "res/refGene.txt.gz"
    refgene_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz"
    
    if os.path.exists(refgene_local_path):
        print(f"Loading RefSeq gene annotations from local cache: {refgene_local_path}...")
        source_path = refgene_local_path
    else:
        print(f"Downloading RefSeq gene annotations from UCSC...")
        print(f"URL: {refgene_url}")
        source_path = refgene_url
        
    try:
        refgene_df = pd.read_csv(
            source_path, 
            sep='\t', 
            header=None, 
            usecols=[2, 4, 5, 12],
            names=['chrom', 'txStart', 'txEnd', 'symbol']
        )
        print("RefSeq database loaded successfully.")
        
        # Cache locally if downloaded from URL
        if source_path == refgene_url:
            print(f"Caching RefSeq database locally to {refgene_local_path}...")
            try:
                urllib.request.urlretrieve(refgene_url, refgene_local_path)
                print("Successfully cached RefSeq database.")
            except Exception as cache_err:
                print(f"Warning: Could not save local cache: {cache_err}")
    except Exception as e:
        print(f"Error loading RefSeq database: {e}")
        return
        
    # Filter chromosomes to match those in our UMAP results
    unique_umap_chroms = set(chrom_list)
    refgene_df = refgene_df[refgene_df['chrom'].isin(unique_umap_chroms)]
    print(f"Filtered to {len(refgene_df)} genes matching target chromosomes.")
    
    # 3. Intersect genomic bins with genes
    print("Mapping genes to genomic bins (intersection)...")
    
    # Group genes by chromosome and sort by txStart
    genes_by_chrom = {}
    for chrom, group in refgene_df.groupby('chrom'):
        # Keep unique records to reduce redundant calculations
        unique_group = group.drop_duplicates(subset=['txStart', 'txEnd', 'symbol'])
        genes_by_chrom[chrom] = sorted(
            list(zip(unique_group['txStart'], unique_group['txEnd'], unique_group['symbol'])),
            key=lambda x: x[0]
        )
        
    # Group bins by chromosome to intersect chromosome-by-chromosome
    bins_by_chrom = {}
    for i in range(n_bins):
        chrom = chrom_list[i]
        if chrom not in bins_by_chrom:
            bins_by_chrom[chrom] = []
        bins_by_chrom[chrom].append((i, start_list[i], end_list[i]))
        
    # List of list of gene symbols for each bin index
    bin_genes = [[] for _ in range(n_bins)]
    
    for chrom, chrom_bins in bins_by_chrom.items():
        chrom_genes = genes_by_chrom.get(chrom, [])
        if not chrom_genes:
            continue
            
        print(f"  Processing {chrom} ({len(chrom_bins)} bins, {len(chrom_genes)} genes)...")
        
        # Intersection using sorted sweep-line
        for bin_idx, b_start, b_end in chrom_bins:
            overlapping = set()
            for g_start, g_end, symbol in chrom_genes:
                if g_start >= b_end:
                    break
                if g_end > b_start:
                    # Strip any version numbers or whitespace from gene symbols
                    clean_symbol = str(symbol).strip().upper()
                    if clean_symbol:
                        overlapping.add(clean_symbol)
            bin_genes[bin_idx] = sorted(list(overlapping))
            
    # Calculate unique genes
    all_unique_genes = sorted(list(set(g for genes in bin_genes for g in genes)))
    print(f"Mapped genes to all bins. Found {len(all_unique_genes)} unique gene symbols.")
    
    # 4. Fetch pathways and summaries for unique genes from MyGene.info in batches
    print("Fetching gene pathways and summaries from MyGene.info...")
    
    # Load existing annotations if available to use as a cache
    existing_gene_data = {}
    if os.path.exists(output_path):
        try:
            print(f"Loading existing annotations from {output_path} to reuse cached gene information...")
            with open(output_path, 'r') as f:
                existing_data = json.load(f)
                existing_gene_data = existing_data.get("gene_data", {})
            print(f"Found cached details for {len(existing_gene_data)} genes.")
        except Exception as cache_err:
            print(f"Warning: Could not load existing annotations cache: {cache_err}")

    gene_data = {}
    genes_to_fetch = []
    
    # Identify which genes are already in the cache vs. need to be fetched
    for sym in all_unique_genes:
        sym_upper = sym.upper()
        if sym_upper in existing_gene_data:
            gene_data[sym_upper] = existing_gene_data[sym_upper]
        else:
            genes_to_fetch.append(sym)
            
    n_to_fetch = len(genes_to_fetch)
    print(f"Need to fetch MyGene.info annotations for {n_to_fetch} new genes (skipped {len(all_unique_genes) - n_to_fetch} already cached).")
    
    batch_size = 500
    
    for i in range(0, n_to_fetch, batch_size):
        batch = genes_to_fetch[i:i+batch_size]
        print(f"  Fetching batch {i//batch_size + 1}/{(n_to_fetch-1)//batch_size + 1} ({len(batch)} genes)...")
        
        try:
            r = requests.post(
                "https://mygene.info/v3/query",
                data={
                    "q": ",".join(batch),
                    "scopes": "symbol",
                    "fields": "symbol,pathway,summary",
                    "species": "human"
                },
                timeout=30
            )
            r.raise_for_status()
            data = r.json()
            
            for hit in data:
                query_symbol = hit.get("query")
                if not query_symbol:
                    continue
                
                symbol = query_symbol.strip().upper()
                if "notfound" in hit:
                    # Keep empty structure
                    if symbol not in gene_data:
                        gene_data[symbol] = {"summary": "", "pathways": []}
                    continue
                    
                if symbol not in gene_data:
                    gene_data[symbol] = {
                        "summary": hit.get("summary", ""),
                        "pathways": []
                    }
                else:
                    if not gene_data[symbol]["summary"] and hit.get("summary"):
                        gene_data[symbol]["summary"] = hit.get("summary")
                        
                # Extract KEGG and Reactome pathways
                pathway_obj = hit.get("pathway")
                if pathway_obj:
                    # KEGG
                    kegg = pathway_obj.get("kegg")
                    if kegg:
                        kegg_list = kegg if isinstance(kegg, list) else [kegg]
                        for p in kegg_list:
                            gene_data[symbol]["pathways"].append({
                                "source": "KEGG",
                                "id": p.get("id"),
                                "name": p.get("name")
                            })
                    # Reactome
                    reactome = pathway_obj.get("reactome")
                    if reactome:
                        react_list = reactome if isinstance(reactome, list) else [reactome]
                        for p in react_list:
                            gene_data[symbol]["pathways"].append({
                                "source": "Reactome",
                                "id": p.get("id"),
                                "name": p.get("name")
                            })
                            
        except Exception as e:
            print(f"  Error fetching batch: {e}")
            # Ensure we still have empty entries for these so they don't break the app
            for sym in batch:
                sym_upper = sym.upper()
                if sym_upper not in gene_data:
                    gene_data[sym_upper] = {"summary": "", "pathways": []}
                    
    # 5. Save results to a single JSON file
    print(f"Saving compiled annotations to {output_path}...")
    compiled_data = {
        "bin_genes": bin_genes,
        "gene_data": gene_data
    }
    
    with open(output_path, 'w') as f:
        json.dump(compiled_data, f, indent=2)
        
    print("\nSUCCESS: All annotations downloaded and cached successfully!")
    print(f"Output file size: {os.path.getsize(output_path) / (1024 * 1024):.2f} MB")
    print("=" * 60)

if __name__ == "__main__":
    main()
