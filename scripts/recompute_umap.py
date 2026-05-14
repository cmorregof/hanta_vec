#!/usr/bin/env python3
"""
Recompute UMAP coordinates for filtered dataset (767 sequences).
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
import umap

def main():
    print("=" * 70)
    print("RECOMPUTING UMAP FOR FILTERED DATASET")
    print("=" * 70)

    # Load aligned embeddings and metadata
    embeddings_dir = Path("results/embeddings")

    embeddings_path = embeddings_dir / "embeddings_level1.npy"
    metadata_path = Path("data/processed/metadata_level1_with_embeddings.tsv")

    embeddings = np.load(embeddings_path)
    metadata = pd.read_csv(metadata_path, sep="\t")

    print(f"📋 Loaded data:")
    print(f"  Embeddings: {embeddings.shape}")
    print(f"  Metadata: {len(metadata)} rows")

    if len(embeddings) != len(metadata):
        print(f"❌ Dimension mismatch!")
        return 1

    # Compute PCA (50 dimensions)
    print(f"\n🔄 Computing PCA (50 dimensions)...")
    pca = PCA(n_components=50, random_state=42)
    embeddings_pca = pca.fit_transform(embeddings)

    variance_explained = np.sum(pca.explained_variance_ratio_)
    print(f"  PCA shape: {embeddings_pca.shape}")
    print(f"  Variance explained: {variance_explained:.3f}")

    # Compute UMAP
    print(f"\n🔄 Computing UMAP (2D)...")
    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        metric='cosine'
    )

    umap_coords = reducer.fit_transform(embeddings_pca)
    print(f"  UMAP shape: {umap_coords.shape}")

    # Save results
    pca_file = embeddings_dir / "pca_50d_level1.npy"
    umap_file = embeddings_dir / "umap_coords_level1.npy"

    np.save(pca_file, embeddings_pca)
    np.save(umap_file, umap_coords)

    print(f"\n💾 Saved dimensionality reduction:")
    print(f"  PCA: {pca_file}")
    print(f"  UMAP: {umap_file}")

    # Show species distribution in UMAP space
    print(f"\n📊 Species distribution in UMAP space:")
    species_counts = metadata['species'].value_counts()
    for species, count in species_counts.head(8).items():
        species_data = metadata[metadata['species'] == species]
        if len(species_data) > 0:
            indices = species_data.index.tolist()
            species_coords = umap_coords[indices]
            x_range = f"{species_coords[:, 0].min():.2f} to {species_coords[:, 0].max():.2f}"
            y_range = f"{species_coords[:, 1].min():.2f} to {species_coords[:, 1].max():.2f}"
            print(f"  {species} ({count}): X({x_range}), Y({y_range})")

    print(f"\n✅ UMAP recomputation complete!")
    print(f"Ready for F3 generation with {len(umap_coords)} coordinates")

    return 0

if __name__ == "__main__":
    exit(main())