#!/usr/bin/env python3
"""
Phase 2: Compute ESM-2 embeddings for three-segment Level 1 dataset.

Three-segment strategy:
- S-segment: Nucleocapsid protein (~430 aa, no truncation needed)
- M-segment: Glycoprotein precursor (~900-1150 aa, some truncation at 1022)
- L-segment: RNA polymerase (~2100 aa, all truncated to 1022)

Features:
- Separate embedding files per segment
- Truncation logging and tracking
- Resume-capable per segment
- SHA256-based caching
- Comprehensive per-segment statistics
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


def process_segment_embeddings(
    segment: str,
    sequences: dict,
    metadata: pd.DataFrame,
    model, tokenizer, device, config,
    logger, embeddings_dir: Path
) -> dict:
    """Process embeddings for a single segment."""

    logger.info(f"\n--- {segment.upper()}-SEGMENT EMBEDDINGS ---")

    # Progress tracking per segment
    progress_file = embeddings_dir / f"embedding_progress_{segment}.json"
    progress = load_progress(progress_file)

    pending = [acc for acc in sequences.keys() if progress.get(acc) != "done"]
    logger.info(f"  Pending: {len(pending)}/{len(sequences)}")
    logger.info(f"  Already done: {len(sequences) - len(pending)}")

    embeddings = {}
    truncation_log = {}
    max_len = config["embeddings"]["max_length"]  # 1022

    # Compute embeddings
    batch_size = config["embeddings"]["batch_size_cpu"]
    failed_count = 0

    with tqdm(total=len(pending), desc=f"{segment}-segment") as pbar:
        for i, accession in enumerate(pending):
            try:
                seq = sequences[accession]
                original_length = len(seq)

                # Track truncation
                if original_length > max_len:
                    truncation_log[accession] = {
                        'original_length': original_length,
                        'truncated_length': max_len
                    }
                    logger.debug(f"  Truncating {accession}: {original_length} → {max_len}")

                emb = get_or_compute_embedding(
                    seq, accession, model, tokenizer, device, config
                )
                embeddings[accession] = emb
                progress[accession] = "done"

            except Exception as e:
                logger.warning(f"Failed to embed {accession}: {e}")
                progress[accession] = "failed"
                failed_count += 1

            # Save progress every batch
            if (i + 1) % batch_size == 0:
                save_progress(progress_file, progress)

            pbar.update(1)

    # Final progress save
    save_progress(progress_file, progress)

    # Assemble embedding matrix
    accessions = list(sequences.keys())
    embedding_matrix = np.array([
        embeddings.get(acc, np.zeros(model.config.hidden_size))
        for acc in accessions
    ])

    # Save segment embeddings
    np.save(embeddings_dir / f"embeddings_{segment}.npy", embedding_matrix)

    # Save metadata subset for this segment
    segment_metadata = metadata[metadata['segment'] == segment.upper()].copy()
    segment_metadata.to_csv(embeddings_dir / f"metadata_{segment}.csv", index=False)

    # Save accession order
    with open(embeddings_dir / f"accessions_{segment}.txt", "w") as f:
        f.write("\n".join(accessions))

    logger.info(f"  ✓ Computed {len(embeddings)} embeddings")
    logger.info(f"  ✓ Failed: {failed_count}")
    logger.info(f"  ✓ Truncated: {len(truncation_log)} sequences")
    logger.info(f"  ✓ Saved to embeddings_{segment}.npy")

    # Segment statistics
    norms = np.linalg.norm(embedding_matrix, axis=1)

    segment_report = {
        'segment': segment.upper(),
        'total_sequences': len(sequences),
        'embeddings_computed': len(embeddings),
        'embeddings_failed': failed_count,
        'truncated_count': len(truncation_log),
        'mean_original_length': float(np.mean([len(seq) for seq in sequences.values()])),
        'mean_final_length': float(np.mean([min(len(seq), max_len) for seq in sequences.values()])),
        'embedding_shape': list(embedding_matrix.shape),
        'has_nan': bool(np.isnan(embedding_matrix).any()),
        'norm_mean': float(norms.mean()),
        'norm_std': float(norms.std()),
        'truncation_log': truncation_log
    }

    return segment_report


def main():
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    log_dir = Path(__file__).parent.parent / "logs"
    logger = setup_logging(log_dir)

    logger.info("="*80)
    logger.info("PHASE 2: COMPUTE THREE-SEGMENT ESM-2 EMBEDDINGS")
    logger.info("="*80)

    config = load_config(config_path)

    # Load three-segment data
    logger.info("\n1. LOAD THREE-SEGMENT DATA")
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load sequences from three FASTA files
    segment_sequences = {}
    for segment in ['S', 'M', 'L']:
        fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
        sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
        segment_sequences[segment] = sequences
        logger.info(f"  {segment}-segment: {len(sequences)} sequences")

    # Load combined metadata
    metadata_path = processed_dir / "metadata_level1.tsv"
    metadata = pd.read_csv(metadata_path, sep="\t")
    logger.info(f"  Combined metadata: {len(metadata)} records")

    # Load model
    logger.info("\n2. LOAD ESM-2 MODEL")
    model, tokenizer, device = load_model(config)
    logger.info(f"  Model: {config['embeddings']['model']}")
    logger.info(f"  Hidden dim: {model.config.hidden_size}")
    logger.info(f"  Max length: {config['embeddings']['max_length']}")

    # Setup output directory
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    embeddings_dir.mkdir(parents=True, exist_ok=True)

    # Process each segment
    logger.info("\n3. COMPUTE EMBEDDINGS PER SEGMENT")
    segment_reports = {}

    for segment in ['S', 'M', 'L']:
        segment_report = process_segment_embeddings(
            segment,
            segment_sequences[segment],
            metadata,
            model, tokenizer, device, config,
            logger, embeddings_dir
        )
        segment_reports[segment] = segment_report

    # Generate comprehensive report
    logger.info("\n4. GENERATE COMPREHENSIVE REPORT")

    total_sequences = sum(len(seqs) for seqs in segment_sequences.values())
    total_computed = sum(report['embeddings_computed'] for report in segment_reports.values())
    total_failed = sum(report['embeddings_failed'] for report in segment_reports.values())
    total_truncated = sum(report['truncated_count'] for report in segment_reports.values())

    comprehensive_report = {
        'model': config['embeddings']['model'],
        'max_length': config['embeddings']['max_length'],
        'total_sequences': total_sequences,
        'total_computed': total_computed,
        'total_failed': total_failed,
        'total_truncated': total_truncated,
        'per_segment': segment_reports,
        'truncation_by_segment': {
            segment: {
                'count': report['truncated_count'],
                'percentage': report['truncated_count'] / report['total_sequences'] * 100
            }
            for segment, report in segment_reports.items()
        }
    }

    # Save comprehensive report
    manifests_dir = Path(__file__).parent.parent / "results" / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    with open(manifests_dir / "embedding_report.json", "w") as f:
        json.dump(comprehensive_report, f, indent=2)

    logger.info("\n" + "="*50)
    logger.info("EMBEDDING REPORT SUMMARY")
    logger.info("="*50)
    logger.info(f"Total sequences processed: {total_sequences}")
    logger.info(f"Total embeddings computed: {total_computed}")
    logger.info(f"Total failed: {total_failed}")
    logger.info(f"Total truncated: {total_truncated}")

    for segment, report in segment_reports.items():
        logger.info(f"\n{segment}-segment:")
        logger.info(f"  Computed: {report['embeddings_computed']}")
        logger.info(f"  Failed: {report['embeddings_failed']}")
        logger.info(f"  Truncated: {report['truncated_count']} ({report['truncated_count']/report['total_sequences']*100:.1f}%)")
        logger.info(f"  Mean length: {report['mean_original_length']:.1f} → {report['mean_final_length']:.1f}")

    # Build cache index
    logger.info("\n5. BUILD CACHE INDEX")
    cache_index = build_cache_index(config, manifests_dir / "cache_index.json")
    logger.info(f"  ✓ Cached {cache_index['total_cached']} embeddings")

    logger.info("\n" + "="*80)
    logger.info("✓ PHASE 2 COMPLETE - THREE SEGMENTS")
    logger.info("="*80)

    if total_failed > 0:
        logger.warning(f"⚠ {total_failed} total embeddings failed")
        return 1

    logger.info("✓ Ready for Phase 3 (reduction & visualization)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
