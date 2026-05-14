#!/usr/bin/env python3
"""
Step 3: Fixed embedding computation script.

FIXED STRATEGY:
- Compute embeddings into a list, save once at end
- Never use np.zeros pre-allocation with incremental fill
- Only save embeddings that were successfully computed
- Handle failures explicitly without zero padding
"""

import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import load_config
from embeddings.esm2 import load_model, get_or_compute_embedding

def load_missing_sequences(embeddings_dir, processed_dir, segment):
    """Load sequences that need embedding computation for a segment."""

    print(f"  📂 Loading missing sequences for {segment}...")

    # Load recomputation plan
    plan_file = embeddings_dir / "recomputation_plan.json"
    if not plan_file.exists():
        raise FileNotFoundError(f"Recomputation plan not found: {plan_file}")

    with open(plan_file, 'r') as f:
        plan = json.load(f)

    missing_accessions = plan['missing_by_segment'][segment]['missing_accessions']
    print(f"    Missing accessions for {segment}: {len(missing_accessions)}")

    # Load sequences from FASTA
    fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
    all_sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}

    # Filter to only missing sequences
    missing_sequences = {acc: all_sequences[acc] for acc in missing_accessions if acc in all_sequences}

    print(f"    Sequences found in FASTA: {len(missing_sequences)}")
    print(f"    Sample missing: {list(missing_sequences.keys())[:3]}")

    return missing_sequences

def compute_embeddings_fixed(segment, missing_sequences, model, tokenizer, device, config):
    """Compute embeddings using FIXED strategy - no zero padding."""

    print(f"  🧮 Computing {segment} embeddings (FIXED STRATEGY)...")

    # Strategy: Collect successful embeddings in lists, no pre-allocation
    successful_embeddings = []
    successful_accessions = []
    failed_accessions = []

    max_len = config["embeddings"]["max_length"]
    truncation_log = {}

    with tqdm(total=len(missing_sequences), desc=f"{segment} embeddings") as pbar:
        for i, (accession, seq) in enumerate(missing_sequences.items()):
            try:
                original_length = len(seq)

                # Track truncation
                if original_length > max_len:
                    truncation_log[accession] = {
                        'original_length': original_length,
                        'truncated_length': max_len
                    }

                # Compute embedding
                embedding = get_or_compute_embedding(
                    seq, accession, model, tokenizer, device, config
                )

                # Validate embedding (no zero vectors allowed)
                norm = np.linalg.norm(embedding)
                if norm > 1e-6:  # Valid embedding
                    successful_embeddings.append(embedding)
                    successful_accessions.append(accession)
                else:
                    print(f"    ⚠️  Zero embedding detected for {accession}, skipping")
                    failed_accessions.append(accession)

                # Progress reporting
                if (i + 1) % 100 == 0:
                    print(f"    Progress: {i+1}/{len(missing_sequences)}, "
                          f"successful: {len(successful_embeddings)}, "
                          f"failed: {len(failed_accessions)}")

            except Exception as e:
                print(f"    ❌ Failed to embed {accession}: {e}")
                failed_accessions.append(accession)

            pbar.update(1)

    print(f"  📊 {segment} Results:")
    print(f"    Successful embeddings: {len(successful_embeddings)}")
    print(f"    Failed embeddings: {len(failed_accessions)}")
    print(f"    Truncated sequences: {len(truncation_log)}")

    if len(successful_embeddings) == 0:
        print(f"    ❌ No valid embeddings computed for {segment}")
        return None, None, None

    # Convert to numpy array ONLY for successful embeddings
    embedding_matrix = np.array(successful_embeddings)
    print(f"    Final embedding matrix shape: {embedding_matrix.shape}")
    print(f"    Mean norm: {np.linalg.norm(embedding_matrix, axis=1).mean():.6f}")

    return embedding_matrix, successful_accessions, {
        'truncation_log': truncation_log,
        'failed_accessions': failed_accessions,
        'successful_count': len(successful_embeddings),
        'failed_count': len(failed_accessions)
    }

def save_embeddings_fixed(embeddings_dir, segment, embedding_matrix, accessions, report):
    """Save embeddings using FIXED strategy - no zero padding."""

    print(f"  💾 Saving {segment} embeddings...")

    # Save new embeddings
    new_file = embeddings_dir / f"embeddings_{segment}_new.npy"
    np.save(new_file, embedding_matrix)

    # Save accession mapping
    mapping_file = embeddings_dir / f"embeddings_{segment}_new_index.json"
    mapping_data = {
        'accession_ids': accessions,
        'embedding_shape': list(embedding_matrix.shape),
        'computation_report': report
    }

    with open(mapping_file, 'w') as f:
        json.dump(mapping_data, f, indent=2)

    print(f"    ✓ Saved embeddings: {new_file}")
    print(f"    ✓ Saved index: {mapping_file}")

def merge_with_clean_embeddings(embeddings_dir, segment):
    """Merge new embeddings with existing clean embeddings."""

    print(f"  🔗 Merging {segment} embeddings...")

    # Load existing clean embeddings
    clean_file = embeddings_dir / f"embeddings_{segment}_clean.npy"
    clean_index_file = embeddings_dir / f"embeddings_{segment}_index.json"

    if not clean_file.exists():
        print(f"    ⚠️  No clean embeddings found, using only new embeddings")
        # Copy new to final
        new_file = embeddings_dir / f"embeddings_{segment}_new.npy"
        new_index_file = embeddings_dir / f"embeddings_{segment}_new_index.json"
        final_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        final_index_file = embeddings_dir / f"embeddings_{segment}_final_index.json"

        np.save(final_file, np.load(new_file))
        with open(new_index_file, 'r') as f:
            data = json.load(f)
        with open(final_index_file, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"    ✓ Created final: {final_file}")
        return

    # Load clean embeddings
    clean_embeddings = np.load(clean_file)
    with open(clean_index_file, 'r') as f:
        clean_index = json.load(f)
    clean_accessions = clean_index['accession_ids']

    # Load new embeddings
    new_file = embeddings_dir / f"embeddings_{segment}_new.npy"
    new_index_file = embeddings_dir / f"embeddings_{segment}_new_index.json"

    if new_file.exists():
        new_embeddings = np.load(new_file)
        with open(new_index_file, 'r') as f:
            new_index = json.load(f)
        new_accessions = new_index['accession_ids']
    else:
        new_embeddings = np.empty((0, clean_embeddings.shape[1]))
        new_accessions = []

    # Combine embeddings (clean + new)
    combined_embeddings = np.vstack([clean_embeddings, new_embeddings])
    combined_accessions = clean_accessions + new_accessions

    print(f"    Clean: {len(clean_accessions)} embeddings")
    print(f"    New: {len(new_accessions)} embeddings")
    print(f"    Combined: {len(combined_accessions)} embeddings")

    # Save final combined embeddings
    final_file = embeddings_dir / f"embeddings_{segment}_final.npy"
    final_index_file = embeddings_dir / f"embeddings_{segment}_final_index.json"

    np.save(final_file, combined_embeddings)

    final_index_data = {
        'accession_ids': combined_accessions,
        'embedding_shape': list(combined_embeddings.shape),
        'clean_count': len(clean_accessions),
        'new_count': len(new_accessions),
        'total_count': len(combined_accessions)
    }

    with open(final_index_file, 'w') as f:
        json.dump(final_index_data, f, indent=2)

    print(f"    ✓ Saved final: {final_file}")
    print(f"    ✓ Saved final index: {final_index_file}")

def main():
    print("=" * 70)
    print("STEP 3: FIXED EMBEDDING COMPUTATION (NO ZERO PADDING)")
    print("=" * 70)

    # Setup paths
    project_root = Path(__file__).parent.parent
    embeddings_dir = project_root / "results" / "embeddings"
    processed_dir = project_root / "data" / "processed"
    config_path = project_root / "config" / "config.yaml"

    # Load config
    config = load_config(config_path)

    # Load model
    print(f"\n🤖 LOADING ESM-2 MODEL")
    print("=" * 50)
    model, tokenizer, device = load_model(config)
    print(f"Model: {config['embeddings']['model']}")
    print(f"Hidden dim: {model.config.hidden_size}")
    print(f"Max length: {config['embeddings']['max_length']}")

    # Process each segment
    print(f"\n⚡ COMPUTING EMBEDDINGS (FIXED STRATEGY)")
    print("=" * 50)

    segment_reports = {}

    for segment in ['S', 'M', 'L']:
        print(f"\n🧬 Processing {segment} segment...")

        # Load missing sequences
        missing_sequences = load_missing_sequences(embeddings_dir, processed_dir, segment)

        if len(missing_sequences) == 0:
            print(f"  ✓ No missing sequences for {segment}")
            segment_reports[segment] = {'status': 'no_missing'}
            continue

        # Compute embeddings using FIXED strategy
        embedding_matrix, accessions, report = compute_embeddings_fixed(
            segment, missing_sequences, model, tokenizer, device, config
        )

        if embedding_matrix is not None:
            # Save new embeddings
            save_embeddings_fixed(embeddings_dir, segment, embedding_matrix, accessions, report)

            # Merge with clean embeddings
            merge_with_clean_embeddings(embeddings_dir, segment)

            segment_reports[segment] = {
                'status': 'completed',
                'computed': report['successful_count'],
                'failed': report['failed_count'],
                'shape': list(embedding_matrix.shape)
            }
        else:
            segment_reports[segment] = {
                'status': 'failed',
                'computed': 0,
                'failed': len(missing_sequences)
            }

    # Final summary
    print(f"\n📋 EMBEDDING COMPUTATION SUMMARY")
    print("=" * 50)

    total_computed = 0
    total_failed = 0

    for segment, report in segment_reports.items():
        status = report['status']
        if status == 'completed':
            computed = report['computed']
            failed = report['failed']
            total_computed += computed
            total_failed += failed
            print(f"{segment}: {computed} computed, {failed} failed")
        elif status == 'no_missing':
            print(f"{segment}: no missing sequences")
        else:
            print(f"{segment}: computation failed")

    print(f"\nTotal: {total_computed} computed, {total_failed} failed")
    print(f"Success rate: {total_computed/(total_computed+total_failed)*100:.1f}%")

    print(f"\n✅ Step 3 complete! Fixed embedding computation finished.")
    print(f"Next: Run Step 4 validation to ensure zero vectors = 0%")

    return 0

if __name__ == "__main__":
    exit(main())