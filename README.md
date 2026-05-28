# ALPHAGenome DNA Embeddings — UMAP Explorer

[![GitHub Pages](https://img.shields.io/badge/🧬_Live_Demo-GitHub_Pages-6366f1?style=for-the-badge)](https://lagosproject.github.io/DNAEmbeddings/)
[![HuggingFace Dataset](https://img.shields.io/badge/🤗_Dataset-HuggingFace-FFD21E?style=for-the-badge)](https://huggingface.co/datasets/lagosproject/ALPHAGenome-Embeddings)

An interactive visualization of latent relationships between genomic regions, built from **ALPHAGenome** DNA foundation model embeddings.

<p align="center">
  <img src="docs/preview.png" alt="UMAP Explorer Preview" width="900">
</p>

---

## 🌐 Live Demo

**→ [Open the UMAP Explorer](https://lagosproject.github.io/DNAEmbeddings/)**

The web app runs entirely in the browser — no server needed. All pre-computed UMAP coordinates and gene/pathway annotations are bundled as static JSON files.

---

## 📊 Data & Dataset

The raw embeddings are available on HuggingFace:

**→ [lagosproject/ALPHAGenome-Embeddings](https://huggingface.co/datasets/lagosproject/ALPHAGenome-Embeddings)**

This dataset contains the pre-computed 3,072-dimensional embeddings for all genomic bins across the human genome (hg38), along with metadata and reference annotations.

---

## 🧬 What This Visualizes

The human genome (hg38 / GRCh38) is divided into **~12,400 non-overlapping 131 KB bins**. Each bin's DNA sequence is embedded into a 3,072-dimensional latent space using the ALPHAGenome foundation model. These embeddings reveal:

- **Chromatin structure clusters** — Regions with similar regulatory profiles cluster together
- **Gene-level patterns** — Averaged bin embeddings per gene, projected to a Gene UMAP
- **Pathway-level patterns** — Averaged gene embeddings per KEGG/Reactome pathway

### Three Mapping Levels

| Level | Points | Description |
|-------|--------|-------------|
| **Genomic Bins** | ~12,400 | Raw 131 KB windows, colored by chromosome or position |
| **Individual Genes** | ~25,000 | Average embedding per RefSeq gene |
| **Biological Pathways** | ~1,200 | Average embedding per KEGG/Reactome pathway |

---

## 📂 Repository Structure

```
DNAEmbeddings/
├── index.html                        # Main web application (single-file, self-contained)
├── res/                              # Pre-computed data (JSON files tracked in git)
│   ├── umap_results.json             # Bin-level UMAP coordinates (global + per-chromosome)
│   ├── gene_umap_results.json        # Gene-level UMAP coordinates
│   ├── pathway_umap_results.json     # Pathway-level UMAP coordinates
│   └── gene_annotations.json         # Gene/pathway annotations (RefSeq + MyGene.info)
├── notebooks/
│   └── generate_embeddings.ipynb     # Colab notebook: run AlphaGenome on hg38 → produces chr*_embeddings.npy
├── scripts/                          # Pipeline and utility scripts
│   ├── download_embeddings.py        # Download pre-computed embeddings from HuggingFace
│   ├── run_umap.py                   # Compute UMAP from raw embeddings
│   ├── run_umap_extensions.py        # Compute gene & pathway UMAP projections
│   ├── download_annotations.py       # Download gene/pathway annotations
│   ├── serve.py                      # Local development server
│   └── rethrow.sh                    # Re-run full pipeline
└── docs/                             # Documentation assets
```

### Files NOT in this repo (too large for GitHub)

The raw embeddings live on HuggingFace: **[lagosproject/ALPHAGenome-Embeddings](https://huggingface.co/datasets/lagosproject/ALPHAGenome-Embeddings)**

| File pattern | Size | Purpose |
|---|---|---|
| `res/chr*_embeddings.npy` | ~250 MB total | Raw 3,072-dim ALPHAGenome embeddings per chromosome |
| `res/chr*_metadata.csv` | ~0.5 MB total | Bin coordinate metadata (chrom, start, end) |
| `res/refGene.txt.gz` | ~8 MB | UCSC RefSeq gene database cache |

Download them with the provided script (see below). These are only needed to **re-run** the UMAP computation pipeline. The web app works without them.

---


## 🛠️ Local Development

### Quick Start (just the viewer)

```bash
python3 scripts/serve.py
```

This opens the interactive explorer at `http://localhost:8000/index.html`.

### Full Pipeline (re-generate everything from embeddings)

> **Want to re-run the model?** Open [`notebooks/generate_embeddings.ipynb`](notebooks/generate_embeddings.ipynb) in Google Colab (GPU recommended). It downloads hg38, runs AlphaGenome on every 131 KB bin, and outputs embeddings.

```bash
# 0. Download embeddings from HuggingFace (~312 MB)
python3 scripts/download_embeddings.py

# Or download only specific chromosomes:
python3 scripts/download_embeddings.py chr1 chr22

# 1. Run UMAP dimensionality reduction
python3 scripts/run_umap.py

# 2. Download gene/pathway annotations
python3 scripts/download_annotations.py

# 3. Compute gene & pathway UMAP projections
python3 scripts/run_umap_extensions.py

# 4. Launch the web viewer
python3 scripts/serve.py

# Or run the full pipeline at once:
bash scripts/rethrow.sh
```

#### Python dependencies:
```
numpy
pandas
umap-learn
requests
huggingface_hub
```

---

## 🎨 Features

- **Interactive UMAP Scatter Plot** — Powered by Plotly.js with zoom, pan, hover inspection
- **Three Mapping Levels** — Toggle between genomic bins, genes, and biological pathways
- **Chromosome Selection** — View all chromosomes or isolate specific ones (single-click toggle, double-click solo)
- **Advanced Highlighting** — Search and highlight by gene symbol, biological pathway, or unannotated "desert" regions
- **Inspector Panel** — Click any point to see coordinates, overlapping genes, functional summaries, and KEGG/Reactome pathways
- **Physical Karyotype** — Bottom panel shows physical genome positions of visible UMAP points
- **Genomic Search** — Search by coordinate range (e.g. `chr3:10,000,000-20,000,000`), gene symbol, or pathway name
- **Fully Offline** — All data is pre-bundled; works without any API calls

---

## 📊 Data Sources

| Source | Usage |
|--------|-------|
| **ALPHAGenome** | DNA sequence embeddings (3,072-dim per 131 KB bin) |
| **UCSC RefSeq** (hg38) | Gene coordinate annotations |
| **MyGene.info** | Gene functional summaries, KEGG & Reactome pathway mappings |
| **UMAP** | Dimensionality reduction (cosine metric, n_neighbors=15, min_dist=0.1) |

---

## 📜 License

This project is provided for research and educational purposes.

---

## 🙏 Acknowledgments

Built with [Plotly.js](https://plotly.com/javascript/), [UMAP](https://umap-learn.readthedocs.io/), and genomic data from [UCSC Genome Browser](https://genome.ucsc.edu/) and [MyGene.info](https://mygene.info/).
