# HantaVec: Latent-Space Audit of Orthohantavirus Glycoprotein Embeddings

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Status: MVP](https://img.shields.io/badge/Status-MVP-orange.svg)]()

## Overview

**Question:** Do pre-trained protein language models (ESM-2) preserve species-level structure in Orthohantavirus glycoproteins despite compressing pairwise similarities?

**Answer:** Yes. ESM-2 embeddings (398 sequences, 480-dim) correlate strongly with sequence identity (ρ=0.689, p<10⁻¹⁴) and recover species clusters in UMAP space, despite 19.6× compression in similarity range.

## Key Results

| Finding | Value | Implication |
|---------|-------|-------------|
| ESM-2 / Identity Correlation | ρ = 0.689 | Strong preservation of evolutionary distance |
| Similarity Compression | 78.3 → 4.0 point range (19.6×) | High-resolution within narrow band |
| Species Clustering | Puumala, Hantaan, Seoul well-separated | Structure recoverable despite compression |
| Conservation | Mean 0.281 ± 0.072 | Moderate variability (expected for RNA virus) |

---

## Quick Start

```bash
# Install
git clone https://github.com/cmorregof/hanta_vec.git
cd hanta_vec && pip install -r requirements.txt

# Run pipeline (4 phases)
python scripts/02_build_dataset.py          # Phase 1: Fetch + QC → 398 sequences
python scripts/03_compute_embeddings.py     # Phase 2: ESM-2 embeddings (cached)
python scripts/04_mvp_figures.py            # Phase 3: PCA/UMAP + figures
jupyter notebook notebooks/05_structure.ipynb  # Phase 4: Conservation + 3D
```

**Load embeddings:**
```python
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

emb = np.load('results/embeddings/embeddings_level1.npy')  # (398, 480)
sim = cosine_similarity(emb)  # pairwise similarities
```

---

## Data & Methods

### Dataset
- **Source:** NCBI (3 parallel sources: protein DB, M-segment, RefSeq)
- **Curation:** 2,115 candidates → 398 after QC (18.8% pass rate)
- **Species:** Puumala (84), Hantaan (57), Seoul (46), Dobrava (34), others (177)
- **Length:** 200–677 aa (mean 437 ± 69)

### Pipeline
1. **ESM-2 (35M):** Mean pooling over residues (exclude [BOS]/[EOS]), SHA256 caching
2. **PCA:** 480 → 50 → 2 dims (42.4% variance in PC1+PC2)
3. **UMAP:** 50D → 2D (n_neighbors=15, min_dist=0.1, cosine metric)
4. **Baselines:** Sequence identity (Bio.Align), conservation scoring (Shannon entropy)

---

## Results

### F1: Dataset Distribution
![Dataset overview](results/figures/small/F1_dataset_overview.png)
Species, lengths, and extraction methods across 398 sequences.

### F2: PCA by Species
![PCA scatter](results/figures/small/F2_pca_species.png)
Clear species separation along PC1; scree inset shows cumulative variance.

### F3/F4: UMAP Clustering
![UMAP by species](results/figures/small/F3_umap_species.png)
Puumala, Hantaan, Seoul form tight clusters; Old/New World partially separated.

### F5: Similarity Heatmap
![Similarity heatmap](results/figures/small/F5_similarity_heatmap.png)
80 stratified sequences; strong within-species blocks despite 0.960–1.000 range.

### F6: ESM-2 vs Sequence Identity
![Correlation](results/figures/small/F6_esm2_vs_identity.png)
**Spearman ρ = 0.689** (N=100, seed=42, p<10⁻¹⁴). Strong monotonic relationship.

### F6b: Compression Analysis
![Compression](results/figures/small/F6b_compression_analysis.png)
Identity range 0.194–0.977 vs ESM-2 range 0.960–1.000 (19.6× compression).

### F7: Conservation Scores
![Conservation](results/figures/small/F7_conservation_scores.png)
MSA of 79 sequences (top 5 species); mean conservation 0.281 (genus-level divergence).

---

## What We Claim / Don't Claim

✅ ESM-2 preserves species structure despite compression  
✅ Correlation with sequence identity quantified (ρ=0.689)  
✅ Fully reproducible (config-driven, seed=42, cached embeddings)  

❌ Pathogenicity, transmissibility, or fitness predictions  
❌ Vaccine design or therapeutic claims  
❌ Fine-tuned model (pre-trained only)  

---

## Reproducibility

- **Config:** All hyperparameters in `config/config.yaml`
- **Manifests:** `results/manifests/{dataset_manifest, qc_report, cache_index}.{tsv,json}`
- **Caching:** SHA256-based, auto-recomputed if missing
- **Seeds:** Python/NumPy/Sklearn all use seed=42

---

## Structure

```
scripts/
  ├── 02_build_dataset.py       # Phase 1: NCBI fetch + QC
  ├── 03_compute_embeddings.py  # Phase 2: ESM-2 inference
  └── 04_mvp_figures.py         # Phase 3: Visualization

src/
  ├── data/fetch.py             # 3-source parallel NCBI fetching
  ├── embeddings/esm2.py        # Model loading + mean pooling
  ├── reduction/{pca,umap}.py   # Dimensionality reduction
  └── baselines/identity.py     # Sequence identity baseline

data/processed/
  ├── gn_sequences_level1.fasta  # 398 curated sequences
  └── metadata_level1.tsv        # Species, country, method

results/
  ├── embeddings/embeddings_level1.npy  # (398, 480)
  └── figures/small/{F1–F7}.png
```

---

## Requirements

```
Python 3.9+, torch, transformers, scikit-learn, umap-learn, biopython, 
matplotlib, seaborn, pandas, jupyter
```

**Hardware:** CPU only (all phases). MPS/GPU optional for Phase 2.  
**Runtime:** ~15 min (Phase 1 NCBI), ~1 min (Phase 2), ~20 sec (Phases 3–4).

---

## Limitations

- Pre-trained only (no fine-tuning on viral sequences)
- Gn domain only (~480 aa, not full GPC)
- Geographic metadata sparse
- Pilot study (not peer-reviewed)
- No cross-genus validation

---

## Citation

```bibtex
@software{hantavec2026,
  title={HantaVec: Latent-Space Audit of Orthohantavirus Glycoprotein Embeddings},
  year={2026},
  url={https://github.com/cmorregof/hanta_vec},
  note={MVP, pre-peer-review}
}
```

**License:** MIT  
**Status:** MVP (Phase 5 ready: fine-tuning + downstream clustering)
