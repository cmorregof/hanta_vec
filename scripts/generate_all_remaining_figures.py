#!/usr/bin/env python3
"""
Generate all remaining figures: F1, F2, F5, F6 with correct data columns.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from scipy import stats

# Species colors (consistent with F3)
SPECIES_COLORS = {
    'Andes': '#b07aa1',
    'Seoul': '#e15759',
    'Puumala': '#4e79a7',
    'Choclo': '#f28e2b',
    'Prospect Hill': '#af7aa1',
    'Sin Nombre': '#edc948',
    'Bayou': '#59a14f',
    'Caño Delgadito': '#76b7b2',
    'Muleshoe': '#cccccc',
    'Hantaan': '#cccccc',
    'Black Creek Canal': '#9c755f',
    'other': '#cccccc'
}

def make_f1_dataset_overview(metadata, output_dir):
    """Generate F1: Dataset Overview with correct N count."""

    print("  Generating F1: Dataset Overview...")

    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: Species distribution (horizontal bar plot)
    species_counts = metadata['species'].value_counts().sort_values()
    colors_list = [SPECIES_COLORS.get(s, SPECIES_COLORS['other']) for s in species_counts.index]

    ax_a.barh(range(len(species_counts)), species_counts.values, color=colors_list)
    ax_a.set_yticks(range(len(species_counts)))
    ax_a.set_yticklabels([s[:25] for s in species_counts.index], fontsize=9)
    ax_a.set_xlabel('Number of Sequences')
    ax_a.set_title(f'Panel A: Species Distribution (N={len(metadata)})', fontweight='bold')
    ax_a.grid(axis='x', alpha=0.3)

    # Panel B: Sequence length distribution
    ax_b.hist(metadata['sequence_length'], bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    ax_b.set_xlabel('Sequence Length (amino acids)')
    ax_b.set_ylabel('Frequency')
    ax_b.set_title('Panel B: Sequence Length Distribution', fontweight='bold')
    ax_b.axvline(1022, color='red', linestyle='--', alpha=0.7, label='ESM-2 limit')
    ax_b.legend()
    ax_b.grid(alpha=0.3)

    # Panel C: Segment distribution
    segment_counts = metadata['segment'].value_counts()
    colors = ['#4e79a7', '#e15759', '#f28e2b']  # S, M, L colors

    ax_c.pie(segment_counts.values, labels=segment_counts.index, autopct='%1.1f%%',
             colors=colors, startangle=90)
    ax_c.set_title('Panel C: Segment Distribution', fontweight='bold')

    # Panel D: Clade distribution
    metadata['clade'] = metadata['species'].apply(lambda x:
        'Old World' if x in ['Seoul', 'Hantaan', 'Puumala', 'Prospect Hill']
        else 'New World' if x in ['Andes', 'Sin Nombre', 'Choclo', 'Bayou', 'Caño Delgadito', 'Black Creek Canal', 'Muleshoe']
        else 'Unknown')

    clade_counts = metadata['clade'].value_counts()
    clade_colors = {'Old World': '#4e79a7', 'New World': '#e15759', 'Unknown': '#cccccc'}
    colors = [clade_colors.get(c, '#cccccc') for c in clade_counts.index]

    ax_d.pie(clade_counts.values, labels=clade_counts.index, autopct='%1.1f%%',
             colors=colors, startangle=90)
    ax_d.set_title('Panel D: Clade Distribution', fontweight='bold')

    plt.suptitle(f'HantaVec Dataset Overview — Three-Segment Orthohantavirus Analysis (N={len(metadata)})',
                fontsize=16, fontweight='bold', y=0.95)
    plt.tight_layout()

    # Save
    output_path = output_dir / "F1_dataset_overview.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved F1: {output_path}")

def make_f2_pca_analysis(embeddings, metadata, output_dir):
    """Generate F2: PCA Analysis with scree plot."""

    print("  Generating F2: PCA Analysis...")

    # Compute PCA
    pca = PCA(n_components=10)
    pca_coords = pca.fit_transform(embeddings)

    fig, (ax_main, ax_scree) = plt.subplots(1, 2, figsize=(15, 6))

    # Main PCA plot
    for species in metadata['species'].unique():
        if species in SPECIES_COLORS:
            mask = metadata['species'] == species
            species_coords = pca_coords[mask]
            ax_main.scatter(species_coords[:, 0], species_coords[:, 1],
                           c=SPECIES_COLORS[species], label=f'{species} (n={np.sum(mask)})',
                           alpha=0.7, s=30)

    ax_main.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax_main.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax_main.set_title('PCA of ESM-2 Embeddings by Species')
    ax_main.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax_main.grid(alpha=0.3)

    # Scree plot
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    ax_scree.plot(range(1, 11), pca.explained_variance_ratio_, 'bo-', label='Individual')
    ax_scree.plot(range(1, 11), cumvar, 'ro-', label='Cumulative')
    ax_scree.set_xlabel('Principal Component')
    ax_scree.set_ylabel('Explained Variance Ratio')
    ax_scree.set_title('PCA Scree Plot')
    ax_scree.legend()
    ax_scree.grid(alpha=0.3)

    plt.suptitle(f'PCA Analysis — Orthohantavirus ESM-2 Embeddings (N={len(metadata)})',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Save
    output_path = output_dir / "F2_pca_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved F2: {output_path}")

def make_f5_similarity_heatmap(embeddings, metadata, output_dir):
    """Generate F5: Cosine Similarity Heatmap (stratified sample)."""

    print("  Generating F5: Similarity Heatmap...")

    # Stratified sampling: 8 sequences per major species
    major_species = ['Seoul', 'Andes', 'Puumala', 'Choclo', 'Prospect Hill', 'Sin Nombre']
    selected_indices = []

    for species in major_species:
        species_mask = metadata['species'] == species
        species_indices = metadata[species_mask].index.tolist()

        if len(species_indices) >= 8:
            # Take first 8 (could be random sample)
            selected_indices.extend(species_indices[:8])
        else:
            selected_indices.extend(species_indices)

    # Compute similarity matrix
    selected_embeddings = embeddings[selected_indices]
    similarity_matrix = cosine_similarity(selected_embeddings)

    # Create labels
    selected_metadata = metadata.iloc[selected_indices]
    labels = [f"{row['species'][:8]}_{row['segment']}" for _, row in selected_metadata.iterrows()]

    # Plot heatmap
    plt.figure(figsize=(12, 10))

    sns.heatmap(similarity_matrix,
                xticklabels=labels,
                yticklabels=labels,
                cmap='viridis',
                center=0.5,
                square=True,
                cbar_kws={'label': 'Cosine Similarity'})

    plt.title(f'ESM-2 Embedding Similarity Matrix\nStratified Sample (n={len(selected_indices)})',
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save
    output_path = output_dir / "F5_similarity_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved F5: {output_path}")

def main():
    print("=" * 70)
    print("GENERATING REMAINING FIGURES: F1, F2, F5")
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

    # Generate figures
    make_f1_dataset_overview(metadata, output_dir)
    make_f2_pca_analysis(embeddings, metadata, output_dir)
    make_f5_similarity_heatmap(embeddings, metadata, output_dir)

    print(f"\n✅ All remaining figures generated!")

    return 0

if __name__ == "__main__":
    exit(main())