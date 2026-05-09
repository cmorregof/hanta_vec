# HantaVec

A scalable, reproducible latent-space and structure-aware explorer for Orthohantavirus glycoproteins.

## Overview

HantaVec builds a curated dataset of Orthohantavirus Gn/Gc/GPC sequences from NCBI GenBank, generates protein language model embeddings (ESM-2), and evaluates whether these embeddings organize sequences by known evolutionary and geographic structure.

**Scope:** This is a descriptive representation audit. We do not:
- Generate or optimize sequences
- Predict pathogenicity, transmissibility, or immune escape
- Design vaccines or therapeutics
- Make causal or predictive claims

## Quick Start

```bash
# 1. Clone and setup
git clone <repo>
cd HantaVec
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# 2. Install dependencies
pip install -r requirements.txt

# 3. Configure (NCBI access)
# Edit config/config.yaml: add your NCBI email and API key

# 4. Run environment check
python scripts/00_check_environment.py

# 5. Probe NCBI for available data
python scripts/01_ncbi_probe.py

# 6. Review Phase 0 report
cat reports/phase0_summary.md
```

## Project Structure

```
HantaVec/
├── config/
│   └── config.yaml          # Single source of truth for all parameters
├── data/
│   ├── raw/                 # Original GenBank/FASTA files (gitignored)
│   ├── processed/           # Curated, QC'd sequences
│   └── structures/          # PDB files
├── src/
│   ├── data/                # Fetch, QC, metadata, deduplication
│   ├── embeddings/          # ESM-2 inference + caching
│   ├── baselines/           # k-mer, identity, alignment baselines
│   ├── reduction/           # PCA, UMAP
│   └── visualization/       # Figures
├── scripts/                 # Standalone executables for each phase
├── notebooks/               # Jupyter notebooks for exploration
├── results/
│   ├── manifests/           # Checksums, QC reports (in Git)
│   ├── embeddings/          # .npy embeddings (gitignored)
│   ├── baselines/           # Baseline matrices (gitignored)
│   └── figures/
│       ├── small/           # PNG figures ≤1MB (in Git)
│       └── large/           # HTML interactive figures (gitignored)
└── tests/                   # Unit tests
```

## Phases

**Phase 0:** Environment check + NCBI probe  
**Phase 1:** Data fetch, QC, deduplication  
**Phase 2:** ESM-2 embeddings + baselines  
**Phase 3:** MVP figures (F1–F6)  
**Phase 4:** 3D structure demo  

Each phase produces a summary report and must pass sanity checks before proceeding.

## Configuration

All parameters live in `config/config.yaml`. Key sections:
- `ncbi`: NCBI API configuration (email, key, batch size)
- `taxa`: NCBI taxonomy IDs for Orthohantavirus species
- `qc`: Quality control thresholds (length, ambiguous AAs, etc.)
- `embeddings`: ESM-2 model selection, pooling strategy, device
- `reduction`: PCA/UMAP hyperparameters
- `paths`: All input/output directories

**Never commit:** NCBI API keys, email addresses. Use environment variables or `.env`.

## Dataset

- **Level 0:** 10–20 sequences (smoke test)
- **Level 1:** 100–500 sequences (MVP, if NCBI has them)
- **Level 2:** Full dataset (paper extension, GPU-enabled)

## Embeddings

- **Model:** ESM-2 35M (CPU-feasible, reproducible)
- **Pooling:** Mean over all residue tokens (excluding special tokens)
- **Caching:** SHA256(sequence + model) → `.npy`
- **Device:** Auto-detect CUDA; fall back to CPU

## Baselines

1. **k-mer composition** (k=2, k=3)
2. **Pairwise sequence identity**
3. **MSA-based alignment distance** (if MAFFT available)

**Figure F6 (key):** ESM-2 cosine similarity vs. pairwise identity scatter. If ESM-2 captures structure beyond identity, points diverge from diagonal.

## Scientific Claims

**Permitted:**
- "ESM-2 embeddings cluster sequences by known species"
- "Embedding similarity shows correlation (ρ = X) with sequence identity"
- "Old World and New World sequences occupy distinct embedding regions"

**Prohibited:**
- Predictive claims (transmissibility, pathogenicity, emergence)
- Causal claims (fitness, evolution)
- Clinical/therapeutic claims
- "ESM-2 is superior to identity" (without rigorous comparison)

## Reproducibility

The pipeline is deterministic given:
1. `config.yaml` (parameters)
2. NCBI database snapshot (via manifests)
3. Random seeds (in config)

To reproduce:
```bash
python scripts/run_pipeline.sh  # (Phase 1+)
```

To audit data:
```bash
cat results/manifests/dataset_manifest_level0.tsv  # SHA256, accession, metadata
```

## Requirements

- Python ≥ 3.9
- 4–8 GB RAM (CPU mode); GPU optional but recommended for Phase 2+
- ~5 GB disk (raw data + results)
- NCBI account (free, for API key)

## License

MIT

## Authors

Implemented as a reproducible bioinformatics project.

## Contributing

This is a pilot project. Issues and PRs welcome.
