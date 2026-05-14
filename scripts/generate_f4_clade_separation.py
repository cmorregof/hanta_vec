#!/usr/bin/env python3
"""
Generate F4: UMAP Old World vs New World clade separation.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def classify_clade(species):
    """Classify species into Old World vs New World clades."""
    old_world = ['Seoul', 'Hantaan', 'Puumala', 'Dobrava', 'Prospect Hill']
    new_world = ['Andes', 'Sin Nombre', 'Choclo', 'Bayou', 'Caño Delgadito', 'Black Creek Canal', 'Muleshoe']

    if species in old_world:
        return 'Old World'
    elif species in new_world:
        return 'New World'
    else:
        return 'Unknown'

def main():
    print("=" * 70)
    print("GENERATING F4: OLD WORLD vs NEW WORLD CLADE SEPARATION")
    print("=" * 70)

    # Load data
    umap_coords = np.load("results/embeddings/umap_coords_level1.npy")
    metadata = pd.read_csv("data/processed/metadata_level1_with_embeddings.tsv", sep="\t")

    print(f"📋 Loaded:")
    print(f"  UMAP coordinates: {umap_coords.shape}")
    print(f"  Metadata: {len(metadata)} rows")

    # Add clade classification
    metadata['clade'] = metadata['species'].apply(classify_clade)

    # Count clades
    clade_counts = metadata['clade'].value_counts()
    print(f"\n📊 Clade distribution:")
    for clade, count in clade_counts.items():
        print(f"  {clade}: {count}")

    # Create UMAP plot by clade
    plt.figure(figsize=(12, 8))

    # Define colors
    clade_colors = {
        'Old World': '#4e79a7',  # Blue
        'New World': '#e15759',  # Red
        'Unknown': '#cccccc'     # Gray
    }

    # Plot each clade
    for clade in ['Old World', 'New World', 'Unknown']:
        if clade in clade_counts and clade_counts[clade] > 0:
            clade_mask = metadata['clade'] == clade
            clade_coords = umap_coords[clade_mask]

            plt.scatter(
                clade_coords[:, 0],
                clade_coords[:, 1],
                c=clade_colors[clade],
                label=f'{clade} (n={clade_counts[clade]})',
                alpha=0.7,
                s=30,
                edgecolors='white',
                linewidth=0.5
            )

    plt.xlabel('UMAP 1', fontsize=12)
    plt.ylabel('UMAP 2', fontsize=12)
    plt.title('UMAP of ESM-2 Embeddings — Orthohantavirus Clade Separation\n' +
              f'Old World vs New World Clades (N={len(metadata)})', fontsize=14)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save figure
    output_dir = Path("results/figures/small")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / "F4_umap_clade_separation.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✅ Saved F4: {output_path}")

    # Generate summary statistics
    print(f"\n📈 Clade separation analysis:")

    for clade in ['Old World', 'New World']:
        if clade in clade_counts:
            clade_mask = metadata['clade'] == clade
            clade_coords = umap_coords[clade_mask]

            x_range = clade_coords[:, 0].max() - clade_coords[:, 0].min()
            y_range = clade_coords[:, 1].max() - clade_coords[:, 1].min()
            centroid_x = clade_coords[:, 0].mean()
            centroid_y = clade_coords[:, 1].mean()

            print(f"  {clade}:")
            print(f"    Centroid: ({centroid_x:.2f}, {centroid_y:.2f})")
            print(f"    Range: X={x_range:.2f}, Y={y_range:.2f}")

    # Calculate inter-clade distance
    if 'Old World' in clade_counts and 'New World' in clade_counts:
        old_world_mask = metadata['clade'] == 'Old World'
        new_world_mask = metadata['clade'] == 'New World'

        old_world_coords = umap_coords[old_world_mask]
        new_world_coords = umap_coords[new_world_mask]

        old_centroid = np.mean(old_world_coords, axis=0)
        new_centroid = np.mean(new_world_coords, axis=0)

        inter_clade_distance = np.linalg.norm(old_centroid - new_centroid)
        print(f"\n  Inter-clade distance: {inter_clade_distance:.2f}")

    print(f"\n✅ F4 clade separation analysis complete!")

    return 0

if __name__ == "__main__":
    exit(main())