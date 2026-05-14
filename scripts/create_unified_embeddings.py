#!/usr/bin/env python3
"""
Create unified embeddings and metadata files for figure generation.
"""

import numpy as np
import pandas as pd
from pathlib import Path

def main():
    print("=" * 70)
    print("CREATING UNIFIED EMBEDDINGS AND METADATA")
    print("=" * 70)

    # Load segment embeddings and metadata
    embeddings_dir = Path("results/embeddings")
    processed_dir = Path("data/processed")

    all_embeddings = []
    all_metadata = []

    for segment in ['S', 'M', 'L']:
        print(f"\n📋 Loading {segment} segment...")

        # Load embeddings
        embedding_file = embeddings_dir / f"embeddings_{segment}.npy"
        embeddings = np.load(embedding_file)
        print(f"  Embeddings: {embeddings.shape}")

        # Load metadata
        metadata_file = embeddings_dir / f"metadata_{segment}.csv"
        metadata = pd.read_csv(metadata_file)
        print(f"  Metadata: {len(metadata)} rows")

        # Append to combined arrays
        all_embeddings.append(embeddings)
        all_metadata.append(metadata)

    # Combine embeddings
    combined_embeddings = np.vstack(all_embeddings)
    print(f"\n✅ Combined embeddings shape: {combined_embeddings.shape}")

    # Combine metadata
    combined_metadata = pd.concat(all_metadata, ignore_index=True)
    print(f"✅ Combined metadata: {len(combined_metadata)} rows")

    # Save unified files
    unified_embeddings_path = embeddings_dir / "embeddings_level1.npy"
    unified_metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    np.save(unified_embeddings_path, combined_embeddings)
    combined_metadata.to_csv(unified_metadata_path, sep="\t", index=False)

    print(f"\n💾 Saved unified files:")
    print(f"  Embeddings: {unified_embeddings_path}")
    print(f"  Metadata: {unified_metadata_path}")

    # Verify dimensions match
    if len(combined_embeddings) == len(combined_metadata):
        print(f"\n✅ Dimensions match: {len(combined_embeddings)} sequences")
    else:
        print(f"\n❌ Dimension mismatch: {len(combined_embeddings)} embeddings vs {len(combined_metadata)} metadata")

    # Show species breakdown
    print(f"\n📊 Species distribution:")
    species_counts = combined_metadata['species'].value_counts()
    for species, count in species_counts.head(10).items():
        print(f"  {species}: {count}")

    return 0

if __name__ == "__main__":
    exit(main())