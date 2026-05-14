#!/usr/bin/env python3
"""
Step 1: Strip zero padding from embedding files.
Extract only valid embeddings and save as clean files.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

def strip_zero_padding(segment, embeddings_dir, metadata):
    """Strip zero padding from a segment embedding file."""

    print(f"\n🧹 STRIPPING ZERO PADDING: {segment} segment")
    print("-" * 50)

    # Load embedding file
    emb_file = embeddings_dir / f"embeddings_{segment}.npy"
    if not emb_file.exists():
        print(f"❌ File not found: {emb_file}")
        return None, None

    embeddings = np.load(emb_file)
    print(f"Original shape: {embeddings.shape}")

    # Find non-zero embeddings
    norms = np.linalg.norm(embeddings, axis=1)
    valid_indices = np.where(norms > 1e-6)[0]
    zero_indices = np.where(norms <= 1e-6)[0]

    print(f"Valid embeddings: {len(valid_indices)}/{len(embeddings)} ({len(valid_indices)/len(embeddings)*100:.1f}%)")
    print(f"Zero embeddings: {len(zero_indices)}/{len(embeddings)} ({len(zero_indices)/len(embeddings)*100:.1f}%)")

    if len(valid_indices) == 0:
        print(f"❌ No valid embeddings found for {segment}")
        return None, None

    # Extract valid embeddings
    clean_embeddings = embeddings[valid_indices]
    print(f"Clean embeddings shape: {clean_embeddings.shape}")

    # Get corresponding metadata for this segment
    seg_metadata = metadata[metadata['segment_used'] == segment].reset_index(drop=True)
    print(f"Metadata entries for {segment}: {len(seg_metadata)}")

    # Extract accession IDs for valid embeddings
    accession_ids = []
    for idx in valid_indices:
        if idx < len(seg_metadata):
            accession_id = seg_metadata.iloc[idx]['genbank_id']
            accession_ids.append(accession_id)
            if len(accession_ids) <= 5:  # Show first few
                print(f"  Valid embedding {idx} → {accession_id}")
        else:
            print(f"  Warning: embedding index {idx} out of bounds for metadata")
            accession_ids.append(f"UNKNOWN_{idx}")

    print(f"Valid accession IDs extracted: {len(accession_ids)}")

    # Save clean embeddings
    clean_file = embeddings_dir / f"embeddings_{segment}_clean.npy"
    np.save(clean_file, clean_embeddings)
    print(f"✓ Saved clean embeddings: {clean_file}")

    # Save accession ID index
    index_file = embeddings_dir / f"embeddings_{segment}_index.json"
    index_data = {
        'accession_ids': accession_ids,
        'valid_indices': valid_indices.tolist(),
        'original_shape': embeddings.shape,
        'clean_shape': clean_embeddings.shape,
        'mean_norm': float(np.mean(norms[valid_indices]))
    }

    with open(index_file, 'w') as f:
        json.dump(index_data, f, indent=2)
    print(f"✓ Saved accession index: {index_file}")

    return clean_embeddings, accession_ids

def main():
    print("=" * 70)
    print("STEP 1: STRIPPING ZERO PADDING FROM EMBEDDING FILES")
    print("=" * 70)

    # Setup paths
    project_root = Path(__file__).parent.parent
    embeddings_dir = project_root / "results" / "embeddings"
    processed_dir = project_root / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    # Load metadata
    if not metadata_path.exists():
        print(f"❌ Metadata not found: {metadata_path}")
        return 1

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded metadata: {len(metadata)} sequences")

    # Report segment distribution
    segment_counts = metadata['segment_used'].value_counts()
    print(f"Metadata segment distribution:")
    for seg, count in segment_counts.items():
        print(f"  {seg}: {count}")

    # Process each segment
    results = {}
    total_valid = 0

    for segment in ['S', 'M', 'L']:
        clean_embs, acc_ids = strip_zero_padding(segment, embeddings_dir, metadata)

        if clean_embs is not None:
            results[segment] = {
                'valid_count': len(acc_ids),
                'accession_ids': acc_ids
            }
            total_valid += len(acc_ids)
        else:
            results[segment] = {
                'valid_count': 0,
                'accession_ids': []
            }

    # Summary report
    print(f"\n📊 ZERO PADDING STRIPPING SUMMARY")
    print("=" * 50)
    print(f"Valid embeddings per segment after stripping:")
    for segment in ['S', 'M', 'L']:
        count = results[segment]['valid_count']
        print(f"  {segment}: {count} valid embeddings")

    print(f"\nTotal valid embeddings: {total_valid}")
    print(f"Original total: {len(metadata)} sequences")
    print(f"Reduction: {len(metadata) - total_valid} sequences need recomputation")

    # Check for overlaps (shouldn't be any)
    all_ids = []
    for segment in ['S', 'M', 'L']:
        all_ids.extend(results[segment]['accession_ids'])

    unique_ids = set(all_ids)
    if len(all_ids) != len(unique_ids):
        print(f"⚠️  Warning: {len(all_ids) - len(unique_ids)} duplicate accession IDs across segments")
    else:
        print(f"✓ No duplicate accession IDs across segments")

    print(f"\n✅ Zero padding stripped! Clean embedding files created.")
    print(f"Next: Run Step 2 to identify what needs recomputation")

    return 0

if __name__ == "__main__":
    exit(main())