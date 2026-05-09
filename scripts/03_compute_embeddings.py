#!/usr/bin/env python3
"""
Phase 2: Compute ESM-2 embeddings for Level 1 dataset.

Features:
- Batch processing with configurable batch size
- Resume-capable (progress tracking)
- SHA256-based caching
- Sanity checks and statistics
"""

import json
import sys
import logging
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import load_config
from embeddings.esm2 import load_model, get_or_compute_embedding
from embeddings.cache import build_cache_index


def setup_logging(log_dir: Path) -> logging.Logger:
    """Setup logging."""
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "phase2.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def load_progress(progress_file: Path) -> Dict[str, str]:
    """Load embedding progress status."""
    if progress_file.exists():
        with open(progress_file) as f:
            return json.load(f)
    return {}


def save_progress(progress_file: Path, progress: Dict[str, str]):
    """Save embedding progress status."""
    progress_file.parent.mkdir(parents=True, exist_ok=True)
    with open(progress_file, "w") as f:
        json.dump(progress, f, indent=2)


def main():
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    log_dir = Path(__file__).parent.parent / "logs"
    logger = setup_logging(log_dir)

    logger.info("="*80)
    logger.info("PHASE 2: COMPUTE ESM-2 EMBEDDINGS")
    logger.info("="*80)

    config = load_config(config_path)

    # Load data
    logger.info("\n1. LOAD DATA")
    fasta_path = Path(__file__).parent.parent / "data" / "processed" / "gn_sequences_level1.fasta"
    metadata_path = Path(__file__).parent.parent / "data" / "processed" / "metadata_level1.tsv"

    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    metadata = pd.read_csv(metadata_path, sep="\t")

    logger.info(f"  Loaded {len(sequences)} sequences")
    logger.info(f"  Loaded {len(metadata)} metadata entries")

    # Load model
    logger.info("\n2. LOAD MODEL")
    model, tokenizer, device = load_model(config)
    logger.info(f"  Model: {config['embeddings']['model']}")
    logger.info(f"  Hidden dim: {model.config.hidden_size}")

    # Initialize progress tracking
    logger.info("\n3. INITIALIZE PROGRESS TRACKING")
    progress_file = Path(__file__).parent.parent / "results" / "manifests" / "embedding_progress.json"
    progress = load_progress(progress_file)

    pending = [acc for acc in sequences.keys() if progress.get(acc) != "done"]
    logger.info(f"  Pending: {len(pending)}/{len(sequences)}")
    logger.info(f"  Already done: {len(sequences) - len(pending)}")

    # Compute embeddings
    logger.info("\n4. COMPUTE EMBEDDINGS")
    embeddings = {}
    batch_size = config["embeddings"]["batch_size_cpu"]

    with tqdm(total=len(pending), desc="Embeddings") as pbar:
        for i, accession in enumerate(pending):
            try:
                seq = sequences[accession]
                emb = get_or_compute_embedding(
                    seq,
                    accession,
                    model,
                    tokenizer,
                    device,
                    config,
                )
                embeddings[accession] = emb
                progress[accession] = "done"
            except Exception as e:
                logger.warning(f"Failed to embed {accession}: {e}")
                progress[accession] = "failed"

            # Save progress every batch
            if (i + 1) % batch_size == 0:
                save_progress(progress_file, progress)

            pbar.update(1)

    # Final progress save
    save_progress(progress_file, progress)

    logger.info(f"  ✓ Computed {len(embeddings)} embeddings")

    # Reconstruct full embedding matrix (in order of accessions)
    logger.info("\n5. ASSEMBLE EMBEDDING MATRIX")
    accessions = list(sequences.keys())
    embedding_matrix = np.array([
        embeddings.get(acc, np.zeros(model.config.hidden_size))
        for acc in accessions
    ])

    logger.info(f"  Shape: {embedding_matrix.shape}")

    # Save embeddings
    logger.info("\n6. SAVE EMBEDDINGS")
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    embeddings_dir.mkdir(parents=True, exist_ok=True)

    np.save(embeddings_dir / "embeddings_level1.npy", embedding_matrix)
    with open(embeddings_dir / "accessions_level1.txt", "w") as f:
        f.write("\n".join(accessions))

    logger.info(f"  ✓ Saved embeddings to {embeddings_dir / 'embeddings_level1.npy'}")
    logger.info(f"  ✓ Saved accessions to {embeddings_dir / 'accessions_level1.txt'}")

    # Sanity checks
    logger.info("\n7. SANITY CHECKS")
    logger.info(f"  Shape: {embedding_matrix.shape}")
    logger.info(f"  Any NaN: {np.isnan(embedding_matrix).any()}")

    norms = np.linalg.norm(embedding_matrix, axis=1)
    logger.info(f"  Embedding norms: {norms.mean():.3f} ± {norms.std():.3f}")
    logger.info(f"  Norm range: [{norms.min():.3f}, {norms.max():.3f}]")

    # Similarity checks
    logger.info("\n8. SIMILARITY SANITY CHECKS")
    from sklearn.metrics.pairwise import cosine_similarity

    sim_matrix = cosine_similarity(embedding_matrix[:min(20, len(embedding_matrix))])
    species_list = metadata.iloc[:20]["species_clean"].values

    # Intra-species similarity
    sp_counts = metadata["species_clean"].value_counts()
    for sp in sp_counts.head(3).index:
        idx = [i for i, s in enumerate(metadata["species_clean"]) if s == sp]
        if len(idx) > 1:
            intra_sim = sim_matrix[np.ix_(
                [i for i, s in enumerate(species_list) if s == sp],
                [i for i, s in enumerate(species_list) if s == sp],
            )]
            if intra_sim.size > 1:
                np.fill_diagonal(intra_sim, np.nan)
                logger.info(f"  {sp[:30]}: intra-sim={np.nanmean(intra_sim):.3f}")

    # Generate report
    logger.info("\n9. GENERATE REPORT")
    report = {
        "total_sequences": len(sequences),
        "embeddings_computed": len(embeddings),
        "embeddings_failed": sum(1 for v in progress.values() if v == "failed"),
        "embedding_shape": list(embedding_matrix.shape),
        "model": config["embeddings"]["model"],
        "has_nan": bool(np.isnan(embedding_matrix).any()),
        "norm_mean": float(norms.mean()),
        "norm_std": float(norms.std()),
        "norm_min": float(norms.min()),
        "norm_max": float(norms.max()),
    }

    manifests_dir = Path(__file__).parent.parent / "results" / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    with open(manifests_dir / "embedding_report.json", "w") as f:
        json.dump(report, f, indent=2)

    logger.info(json.dumps(report, indent=2))

    # Build cache index
    logger.info("\n10. BUILD CACHE INDEX")
    cache_index = build_cache_index(config, manifests_dir / "cache_index.json")
    logger.info(f"  ✓ Cached {cache_index['total_cached']} embeddings")

    logger.info("\n" + "="*80)
    logger.info("✓ PHASE 2 COMPLETE")
    logger.info("="*80)

    if report["embeddings_failed"] > 0:
        logger.warning(f"⚠ {report['embeddings_failed']} embeddings failed")
        return 1

    logger.info("✓ Ready for Phase 3 (reduction & visualization)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
