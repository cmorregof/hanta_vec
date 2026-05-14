#!/usr/bin/env python3
"""
Generate missing F7 (3D visualization) and F8 (conservation analysis) figures.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
from sklearn.manifold import TSNE
import seaborn as sns

def make_f7_3d_visualization(coords_pca_50d, metadata, output_dir):
    """Generate F7: 3D PCA + t-SNE visualization."""

    # 3D PCA
    pca_3d = coords_pca_50d[:, :3]

    # 3D t-SNE
    print("  Computing 3D t-SNE...")
    tsne_3d = TSNE(n_components=3, random_state=42, perplexity=30, max_iter=1000)
    tsne_coords = tsne_3d.fit_transform(coords_pca_50d)

    # Create combined 3D visualization
    fig = plt.figure(figsize=(20, 8))

    # PCA 3D subplot
    ax1 = fig.add_subplot(121, projection='3d')

    species_colors = {}
    colors = plt.cm.tab10(np.linspace(0, 1, len(metadata['species'].unique())))
    for i, species in enumerate(metadata['species'].unique()):
        species_colors[species] = colors[i]

    for species in metadata['species'].unique():
        mask = metadata['species'] == species
        ax1.scatter(
            pca_3d[mask, 0], pca_3d[mask, 1], pca_3d[mask, 2],
            c=[species_colors[species]], label=species, s=50, alpha=0.7
        )

    ax1.set_xlabel('PC1')
    ax1.set_ylabel('PC2')
    ax1.set_zlabel('PC3')
    ax1.set_title('3D PCA by Species')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # t-SNE 3D subplot
    ax2 = fig.add_subplot(122, projection='3d')

    for species in metadata['species'].unique():
        mask = metadata['species'] == species
        ax2.scatter(
            tsne_coords[mask, 0], tsne_coords[mask, 1], tsne_coords[mask, 2],
            c=[species_colors[species]], label=species, s=50, alpha=0.7
        )

    ax2.set_xlabel('t-SNE 1')
    ax2.set_ylabel('t-SNE 2')
    ax2.set_zlabel('t-SNE 3')
    ax2.set_title('3D t-SNE by Species')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Save static version
    output_path = output_dir / "small" / "F7_3d_visualization.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved {output_path}")

    # Create interactive 3D plot with plotly
    df_plot = pd.DataFrame({
        'PC1': pca_3d[:, 0],
        'PC2': pca_3d[:, 1],
        'PC3': pca_3d[:, 2],
        'species': metadata['species'].values,
        'clade': metadata['clade'].values,
        'accession': metadata['genbank_id'].values
    })

    fig_interactive = px.scatter_3d(
        df_plot, x='PC1', y='PC2', z='PC3',
        color='species',
        hover_data=['accession', 'clade'],
        title='Interactive 3D PCA Visualization by Species',
        width=900, height=700
    )

    # Save interactive version
    output_path_html = output_dir / "large" / "F7_3d_visualization.html"
    fig_interactive.write_html(output_path_html)
    print(f"✓ Saved {output_path_html}")

def make_f8_conservation_analysis(embeddings, metadata, sequences, output_dir):
    """Generate F8: Conservation and variability analysis."""

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Panel A: Per-position embedding variance across species
    print("  Computing per-position variance...")

    species_embeddings = {}
    for species in metadata['species'].unique():
        if species != 'Unknown':
            mask = metadata['species'] == species
            species_embeddings[species] = embeddings[mask]

    # Compute mean embedding per species
    species_means = {}
    for species, embs in species_embeddings.items():
        if len(embs) > 0:
            species_means[species] = np.mean(embs, axis=0)

    # Variance across species means
    if len(species_means) > 1:
        all_means = np.array(list(species_means.values()))
        position_variance = np.var(all_means, axis=0)

        ax1.plot(position_variance, color='darkblue', linewidth=1)
        ax1.set_xlabel('Embedding Dimension')
        ax1.set_ylabel('Variance')
        ax1.set_title('Per-Dimension Variance Across Species')
        ax1.grid(True, alpha=0.3)

    # Panel B: Intra vs Inter-species distances
    print("  Computing species distance distributions...")

    intra_distances = []
    inter_distances = []

    # Sample for computational efficiency
    max_samples = 50
    for species in species_embeddings.keys():
        species_embs = species_embeddings[species]
        if len(species_embs) > 1:
            # Sample intra-species distances
            indices = np.random.choice(len(species_embs), min(max_samples, len(species_embs)), replace=False)
            sample_embs = species_embs[indices]

            for i in range(len(sample_embs)):
                for j in range(i+1, len(sample_embs)):
                    dist = np.linalg.norm(sample_embs[i] - sample_embs[j])
                    intra_distances.append(dist)

    # Sample inter-species distances
    species_list = list(species_means.keys())
    for i in range(len(species_list)):
        for j in range(i+1, len(species_list)):
            dist = np.linalg.norm(species_means[species_list[i]] - species_means[species_list[j]])
            inter_distances.append(dist)

    if intra_distances and inter_distances:
        ax2.hist(intra_distances, bins=30, alpha=0.7, label='Intra-species', color='lightblue', density=True)
        ax2.hist(inter_distances, bins=30, alpha=0.7, label='Inter-species', color='salmon', density=True)
        ax2.set_xlabel('Euclidean Distance')
        ax2.set_ylabel('Density')
        ax2.set_title('Intra vs Inter-species Distances')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    # Panel C: Clade separation in embedding space
    print("  Analyzing clade separation...")

    old_world_mask = metadata['clade'] == 'Old World'
    new_world_mask = metadata['clade'] == 'New World'

    if np.any(old_world_mask) and np.any(new_world_mask):
        old_world_embs = embeddings[old_world_mask]
        new_world_embs = embeddings[new_world_mask]

        # Project to 2D for visualization
        from sklearn.decomposition import PCA
        pca_2d = PCA(n_components=2)
        all_embs_2d = pca_2d.fit_transform(embeddings)

        ax3.scatter(all_embs_2d[old_world_mask, 0], all_embs_2d[old_world_mask, 1],
                   c='red', alpha=0.6, s=30, label=f'Old World (n={np.sum(old_world_mask)})')
        ax3.scatter(all_embs_2d[new_world_mask, 0], all_embs_2d[new_world_mask, 1],
                   c='blue', alpha=0.6, s=30, label=f'New World (n={np.sum(new_world_mask)})')

        ax3.set_xlabel(f'PC1 ({pca_2d.explained_variance_ratio_[0]:.1%} var)')
        ax3.set_ylabel(f'PC2 ({pca_2d.explained_variance_ratio_[1]:.1%} var)')
        ax3.set_title('Clade Separation in PCA Space')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

    # Panel D: Species clustering dendrogram
    print("  Creating species dendrogram...")

    if len(species_means) > 2:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist, squareform

        species_names = list(species_means.keys())
        species_matrix = np.array(list(species_means.values()))

        # Compute pairwise distances
        distances = pdist(species_matrix, metric='euclidean')
        linkage_matrix = linkage(distances, method='ward')

        dendrogram(linkage_matrix, labels=species_names, ax=ax4, orientation='left')
        ax4.set_title('Species Clustering (Ward Linkage)')
        ax4.set_xlabel('Distance')

    plt.tight_layout()

    # Save figure
    output_path = output_dir / "small" / "F8_conservation_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved {output_path}")

def main():
    print("="*60)
    print("GENERATING F7 AND F8 FIGURES")
    print("="*60)

    # Load data
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    figures_dir = Path(__file__).parent.parent / "results" / "figures"

    # Check if we have aligned data or need to use original
    aligned_meta_path = Path(__file__).parent.parent / "data" / "processed" / "metadata_level1_with_embeddings.tsv"
    if aligned_meta_path.exists():
        print("Using aligned metadata...")
        metadata = pd.read_csv(aligned_meta_path, sep="\t")
        embeddings = np.load(embeddings_dir / "embeddings_level1.npy")

        # Load sequences
        sequences = {}
        fasta_path = Path(__file__).parent.parent / "data" / "processed" / "sequences_level1_for_visualization.fasta"
        if fasta_path.exists():
            from Bio import SeqIO
            sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    else:
        print("❌ No aligned data found - run embedding consolidation first")
        return 1

    print(f"Loaded {len(embeddings)} embeddings and {len(metadata)} metadata entries")

    # Load PCA coordinates
    pca_coords_path = embeddings_dir / "pca_coords_level1.npy"
    if pca_coords_path.exists():
        coords_pca = np.load(pca_coords_path)
        print(f"Loaded PCA coords: {coords_pca.shape}")
    else:
        print("❌ No PCA coordinates found")
        return 1

    # Ensure output directories exist
    (figures_dir / "small").mkdir(parents=True, exist_ok=True)
    (figures_dir / "large").mkdir(parents=True, exist_ok=True)

    # Generate F7
    print("\nGenerating F7: 3D Visualization...")
    # For 50D coords, we need to recompute or use the first 50 dims
    coords_50d = coords_pca if coords_pca.shape[1] >= 50 else np.hstack([coords_pca, np.zeros((len(coords_pca), 50-coords_pca.shape[1]))])
    make_f7_3d_visualization(coords_50d[:, :50], metadata, figures_dir)

    # Generate F8
    print("\nGenerating F8: Conservation Analysis...")
    make_f8_conservation_analysis(embeddings, metadata, sequences, figures_dir)

    print("\n✅ F7 and F8 generation complete!")

    return 0

if __name__ == "__main__":
    exit(main())