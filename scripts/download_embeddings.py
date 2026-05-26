"""
Download ALPHAGenome embeddings from HuggingFace into res/.

Usage:
    python3 scripts/download_embeddings.py           # all chromosomes
    python3 scripts/download_embeddings.py chr1 chr2 # specific chromosomes

Requires: pip install huggingface_hub
"""

import sys
import os
from pathlib import Path
from huggingface_hub import hf_hub_download, snapshot_download

REPO_ID  = "lagosproject/ALPHAGenome-Embeddings"
RES_DIR  = Path(__file__).parent.parent / "res"
ALL_CHRS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

RES_DIR.mkdir(exist_ok=True)

chroms = sys.argv[1:] if len(sys.argv) > 1 else ALL_CHRS

if chroms == ALL_CHRS:
    print(f"Downloading all chromosomes (~312 MB) from {REPO_ID} ...")
    snapshot_download(
        REPO_ID,
        repo_type="dataset",
        local_dir=str(RES_DIR),
        allow_patterns=["data/chr*"],
    )
    # Move files from res/data/ to res/
    data_subdir = RES_DIR / "data"
    if data_subdir.exists():
        for f in data_subdir.iterdir():
            f.rename(RES_DIR / f.name)
        data_subdir.rmdir()
else:
    for chrom in chroms:
        for suffix in ("_embeddings.npy", "_metadata.csv"):
            filename = f"{chrom}{suffix}"
            dest = RES_DIR / filename
            if dest.exists():
                print(f"  {filename} already present, skipping.")
                continue
            print(f"  Downloading {filename} ...")
            path = hf_hub_download(
                REPO_ID,
                filename=f"data/{filename}",
                repo_type="dataset",
                local_dir=str(RES_DIR),
            )
            # hf_hub_download may put it in a subdir — move if needed
            downloaded = Path(path)
            if downloaded != dest:
                downloaded.rename(dest)

print(f"\nFiles saved to {RES_DIR}")
