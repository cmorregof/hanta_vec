# Data Directory

## Structure

```
data/
├── raw/              # NCBI raw downloads (gitignored)
├── processed/        # Curated sequences (gitignored)
├── structures/       # PDB files (gitignored)
├── manifests/        # IN REPO: checksums & metadata
└── README.md         # This file
```

## Directory Details

### `raw/` — NCBI Downloads
- **Contents:** Raw GenBank/FASTA files from NCBI Virus, nucleotide, RefSeq
- **Status:** Gitignored (large, re-fetchable)
- **Re-create:** Run `python scripts/02_build_dataset.py`

### `processed/` — Curated Sequences
- **Contents:** 
  - `gn_sequences_level1.fasta` — 398 cleaned, deduplicated Gn sequences
  - `metadata_level1.tsv` — Species, country, extraction method
- **Status:** Gitignored (auto-generated from Phase 1)
- **Re-create:** Run `python scripts/02_build_dataset.py`

### `structures/` — PDB Files
- **Contents:** 6Y6P.pdb, 5LK2.pdb, 6YRQ.pdb, 6Y6Q.pdb (Orthohantavirus structures)
- **Status:** Gitignored (downloaded from RCSB PDB)
- **Re-create:** Run Phase 4 notebook or download from https://www.rcsb.org/

### `manifests/` — Reproducibility Checksums (IN REPO)
- **Contents:**
  - `dataset_manifest_level1.tsv` — Accessions, sequences, metadata checksums
  - `qc_report.json` — QC statistics (pass/fail counts, dedupe info)
  - `cache_index.json` — ESM-2 embedding cache registry
  - `embedding_report.json` — Embedding sanity checks
- **Status:** Committed to Git
- **Purpose:** Audit trail for reproducibility

## Workflow

1. **Phase 1:** `02_build_dataset.py` fetches from NCBI → `raw/`, extracts → `processed/`, writes manifest → `manifests/`
2. **Phase 2–4:** Downstream scripts read from `processed/` and `structures/`, write results to `results/`

## Reproducibility Notes

- To re-fetch all data: `python scripts/02_build_dataset.py` (requires NCBI API key in `config/config.yaml`)
- To verify: Compare SHA256 hashes in `manifests/dataset_manifest_level1.tsv`
- All gitignored directories are auto-regenerated on first run
