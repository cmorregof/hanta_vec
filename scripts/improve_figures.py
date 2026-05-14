#!/usr/bin/env python3
"""
Improve figure quality with better UMAP parameters and visualization.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.decomposition import PCA
import umap
from matplotlib.patches import Ellipse
from scipy.spatial import ConvexHull

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

def make_improved_f4(embeddings, metadata, output_dir):
    """Generate improved F4 with better UMAP and visualization."""

    print("  Generating improved F4: UMAP Clade Separation...")

    # Add clade classification
    metadata['clade'] = metadata['species'].apply(classify_clade)

    # Optimize UMAP parameters for better clade separation
    print("    Computing optimized UMAP...")

    # First PCA to reduce noise
    pca = PCA(n_components=50, random_state=42)
    embeddings_pca = pca.fit_transform(embeddings)

    # UMAP with parameters optimized for clade separation
    reducer = umap.UMAP(
        n_neighbors=30,        # Increased for better global structure
        min_dist=0.01,         # Decreased for tighter clusters
        n_components=2,
        random_state=42,
        metric='cosine',
        spread=2.0,            # Better cluster separation
        n_epochs=500           # More iterations for convergence
    )

    umap_coords = reducer.fit_transform(embeddings_pca)

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Left panel: By species
    species_colors = {
        'Seoul': '#e15759', 'Andes': '#4e79a7', 'Puumala': '#59a14f',
        'Choclo': '#f28e2b', 'Prospect Hill': '#af7aa1', 'Sin Nombre': '#edc948',
        'Bayou': '#76b7b2', 'Caño Delgadito': '#9c755f', 'Muleshoe': '#cccccc',
        'Hantaan': '#999999', 'Black Creek Canal': '#666666'
    }

    for species in metadata['species'].unique():
        if species in species_colors:
            mask = metadata['species'] == species
            coords = umap_coords[mask]
            if len(coords) > 0:
                ax1.scatter(coords[:, 0], coords[:, 1],
                           c=species_colors[species],
                           label=f'{species} (n={len(coords)})',
                           alpha=0.8, s=40, edgecolors='white', linewidth=0.5)

    ax1.set_xlabel('UMAP 1', fontsize=12)
    ax1.set_ylabel('UMAP 2', fontsize=12)
    ax1.set_title('UMAP by Species', fontsize=14, fontweight='bold')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Right panel: By clade with enhanced visualization
    clade_colors = {'Old World': '#2E86AB', 'New World': '#A23B72', 'Unknown': '#cccccc'}

    for clade in ['Old World', 'New World']:
        mask = metadata['clade'] == clade
        coords = umap_coords[mask]
        if len(coords) > 5:  # Need enough points for convex hull
            # Plot points
            ax2.scatter(coords[:, 0], coords[:, 1],
                       c=clade_colors[clade],
                       label=f'{clade} (n={len(coords)})',
                       alpha=0.7, s=50, edgecolors='white', linewidth=1)

            # Add convex hull
            if len(coords) >= 3:
                try:
                    hull = ConvexHull(coords)
                    for simplex in hull.simplices:
                        ax2.plot(coords[simplex, 0], coords[simplex, 1],
                               color=clade_colors[clade], alpha=0.3, linewidth=2)
                except:
                    pass

            # Add density contours
            try:
                x, y = coords[:, 0], coords[:, 1]
                ax2.scatter(np.mean(x), np.mean(y),
                           c=clade_colors[clade], s=200, marker='*',
                           edgecolors='black', linewidth=2, alpha=0.9)
            except:
                pass

    ax2.set_xlabel('UMAP 1', fontsize=12)
    ax2.set_ylabel('UMAP 2', fontsize=12)
    ax2.set_title('UMAP by Evolutionary Clade', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)

    plt.suptitle('Improved UMAP Analysis — Orthohantavirus Evolutionary Structure',
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()

    # Save improved figure
    output_path = output_dir / "F4_improved_umap_clades.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    ✓ Saved improved F4: {output_path}")

    # Calculate separation metrics
    old_world_mask = metadata['clade'] == 'Old World'
    new_world_mask = metadata['clade'] == 'New World'

    old_coords = umap_coords[old_world_mask]
    new_coords = umap_coords[new_world_mask]

    old_centroid = np.mean(old_coords, axis=0)
    new_centroid = np.mean(new_coords, axis=0)

    separation = np.linalg.norm(old_centroid - new_centroid)
    old_spread = np.std(old_coords, axis=0)
    new_spread = np.std(new_coords, axis=0)

    print(f"    Clade separation metrics:")
    print(f"      Inter-clade distance: {separation:.3f}")
    print(f"      Old World spread: {old_spread[0]:.3f}, {old_spread[1]:.3f}")
    print(f"      New World spread: {new_spread[0]:.3f}, {new_spread[1]:.3f}")

    return umap_coords

def make_improved_f3(embeddings, metadata, umap_coords, output_dir):
    """Generate improved F3 with better species visualization."""

    print("  Generating improved F3: Species clustering...")

    plt.figure(figsize=(14, 10))

    # Enhanced species colors
    species_colors = {
        'Seoul': '#e15759', 'Andes': '#4e79a7', 'Puumala': '#59a14f',
        'Choclo': '#f28e2b', 'Prospect Hill': '#af7aa1', 'Sin Nombre': '#edc948',
        'Bayou': '#76b7b2', 'Caño Delgadito': '#9c755f', 'Muleshoe': '#bab0ab',
        'Hantaan': '#e377c2', 'Black Creek Canal': '#7f7f7f'
    }

    # Plot each species with enhanced visualization
    for species in sorted(metadata['species'].unique()):
        mask = metadata['species'] == species
        coords = umap_coords[mask]

        if len(coords) > 0:
            color = species_colors.get(species, '#cccccc')

            # Main scatter plot
            scatter = plt.scatter(coords[:, 0], coords[:, 1],
                                c=color,
                                label=f'{species} (n={len(coords)})',
                                alpha=0.8, s=60, edgecolors='white', linewidth=1)

            # Add convex hull for major species
            if len(coords) >= 4:
                try:
                    hull = ConvexHull(coords)
                    hull_points = coords[hull.vertices]
                    hull_points = np.append(hull_points, [hull_points[0]], axis=0)
                    plt.fill(hull_points[:, 0], hull_points[:, 1],
                           color=color, alpha=0.15)
                    plt.plot(hull_points[:, 0], hull_points[:, 1],
                           color=color, linewidth=2, alpha=0.6)
                except:
                    pass

    plt.xlabel('UMAP 1', fontsize=14)
    plt.ylabel('UMAP 2', fontsize=14)
    plt.title('UMAP of ESM-2 Embeddings — Orthohantavirus Species Clusters\n' +
             f'Clear Species Separation with Convex Hulls (N={len(metadata)})',
             fontsize=16, fontweight='bold')

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save
    output_path = output_dir / "F3_improved_species_umap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    ✓ Saved improved F3: {output_path}")

def main():
    print("=" * 70)
    print("IMPROVING FIGURE QUALITY - BETTER VISUALIZATION")
    print("=" * 70)

    # Load data
    embeddings = np.load("results/embeddings/embeddings_level1.npy")
    metadata = pd.read_csv("data/processed/metadata_level1_with_embeddings.tsv", sep="\t")

    print(f"📋 Loaded:")
    print(f"  Embeddings: {embeddings.shape}")
    print(f"  Metadata: {len(metadata)} rows")

    # Ensure output directory exists
    output_dir = Path("results/figures/small")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate improved figures with better UMAP
    umap_coords = make_improved_f4(embeddings, metadata, output_dir)
    make_improved_f3(embeddings, metadata, umap_coords, output_dir)

    print(f"\n✅ Improved figures generated!")
    print(f"Key improvements:")
    print(f"  - Optimized UMAP parameters for clade separation")
    print(f"  - Enhanced color schemes and visualization")
    print(f"  - Convex hulls for cluster boundaries")
    print(f"  - Better species separation display")

    return 0

if __name__ == "__main__":
    exit(main())