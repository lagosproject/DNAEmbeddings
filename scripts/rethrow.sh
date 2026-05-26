#!/usr/bin/env bash
# Run from project root
cd "$(dirname "$0")/.." || exit 1

python3 scripts/run_umap.py
python3 scripts/download_annotations.py
python3 scripts/run_umap_extensions.py
