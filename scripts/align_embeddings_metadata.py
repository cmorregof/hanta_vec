#!/usr/bin/env python3
"""
Align embeddings with filtered metadata (remove embeddings for excluded sequences).
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path

def main():
    print("=" * 70)
    print("ALIGNING EMBEDDINGS WITH FILTERED METADATA")
    print("=" * 70)

    # Load combined metadata
    metadata_file = Path("data/processed/metadata_level1_with_embeddings.tsv")
    metadata = pd.read_csv(metadata_file, sep="\t")
    print(f"📋 Metadata: {len(metadata)} sequences")

    # Load segment embeddings and indices
    embeddings_dir = Path("results/embeddings")

    aligned_embeddings = []
    total_kept = 0

    for segment in ['S', 'M', 'L']:
        print(f"\n📋 Processing {segment} segment...")

        # Load embeddings
        embedding_file = embeddings_dir / f"embeddings_{segment}.npy"
        embeddings = np.load(embedding_file)

        # Load index
        index_file = embeddings_dir / f"embeddings_{segment}_index.json"
        with open(index_file, 'r') as f:
            index = json.load(f)
        embedding_accessions = index['accession_ids']

        print(f"  Embeddings shape: {embeddings.shape}")
        print(f"  Index accessions: {len(embedding_accessions)}")

        # Get metadata for this segment
        segment_metadata = metadata[metadata['segment'] == segment]
        metadata_accessions = segment_metadata['genbank_id'].tolist()

        print(f"  Metadata accessions: {len(metadata_accessions)}")

        # Find which embeddings to keep
        keep_indices = []
        for i, acc in enumerate(embedding_accessions):
            if acc in metadata_accessions:
                keep_indices.append(i)

        # Filter embeddings
        filtered_embeddings = embeddings[keep_indices]
        aligned_embeddings.append(filtered_embeddings)

        print(f"  Kept: {len(keep_indices)}/{len(embeddings)} embeddings")
        total_kept += len(keep_indices)

    # Combine aligned embeddings
    final_embeddings = np.vstack(aligned_embeddings)
    print(f"\n✅ Final embeddings shape: {final_embeddings.shape}")
    print(f"✅ Metadata shape: {len(metadata)}")

    if len(final_embeddings) == len(metadata):
        print(f"✅ Perfect alignment: {len(final_embeddings)} sequences")
    else:
        print(f"❌ Still misaligned: {len(final_embeddings)} embeddings vs {len(metadata)} metadata")

    # Save aligned files
    aligned_embeddings_path = embeddings_dir / "embeddings_level1.npy"
    np.save(aligned_embeddings_path, final_embeddings)

    print(f"\n💾 Saved aligned embeddings: {aligned_embeddings_path}")
    print(f"   Shape: {final_embeddings.shape}")

    # Verify no zero vectors
    zero_count = sum(1 for emb in final_embeddings if np.linalg.norm(emb) < 1e-6)
    print(f"   Zero vectors: {zero_count}")

    return 0

if __name__ == "__main__":
    exit(main())