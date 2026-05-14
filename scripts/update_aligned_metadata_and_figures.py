#!/usr/bin/env python3
"""
Update aligned metadata with reclassified species and regenerate F1 and F4 figures.
"""

import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def update_aligned_metadata():
    """Update the aligned metadata with reclassified species."""

    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load updated main metadata
    main_metadata = pd.read_csv(processed_dir / "metadata_level1.tsv", sep="\t")

    # Load aligned metadata
    aligned_metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"
    aligned_metadata = pd.read_csv(aligned_metadata_path, sep="\t")

    print(f"Updating aligned metadata with reclassified species...")
    print(f"Main metadata: {len(main_metadata)} sequences")
    print(f"Aligned metadata: {len(aligned_metadata)} sequences")

    # Update species and clade in aligned metadata based on genbank_id
    updated_count = 0

    for idx, row in aligned_metadata.iterrows():
        genbank_id = row['genbank_id']

        # Find corresponding entry in main metadata
        main_row = main_metadata[main_metadata['genbank_id'] == genbank_id]

        if not main_row.empty:
            new_species = main_row.iloc[0]['species']
            new_clade = main_row.iloc[0]['clade']

            if row['species'] != new_species:
                aligned_metadata.at[idx, 'species'] = new_species
                aligned_metadata.at[idx, 'clade'] = new_clade
                updated_count += 1

    print(f"Updated {updated_count} entries in aligned metadata")

    # Save updated aligned metadata
    aligned_metadata.to_csv(aligned_metadata_path, sep="\t", index=False)
    print(f"✓ Saved updated aligned metadata")

    # Report new distribution
    print(f"\nUpdated aligned species distribution:")
    species_counts = aligned_metadata['species'].value_counts()
    for species, count in species_counts.items():
        print(f"  {species}: {count}")

    print(f"\nUpdated aligned clade distribution:")
    clade_counts = aligned_metadata['clade'].value_counts()
    total = len(aligned_metadata)
    for clade, count in clade_counts.items():
        pct = count/total*100
        print(f"  {clade}: {count} ({pct:.1f}%)")

    return aligned_metadata

def regenerate_f1_f4_figures():
    """Regenerate F1 and F4 figures with updated species classification."""

    print(f"\nRegenerating F1 and F4 figures...")

    # Import figure generation functions
    sys.path.insert(0, str(Path(__file__).parent))

    # Read the F1 and F4 functions from the main figure script
    import importlib.util
    import numpy as np

    # Load the figure script
    figure_script_path = Path(__file__).parent / "04_mvp_figures.py"
    spec = importlib.util.spec_from_file_location("mvp_figures", figure_script_path)
    mvp_figures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mvp_figures)

    # Load data
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    figures_dir = Path(__file__).parent.parent / "results" / "figures"

    # Load aligned data
    metadata = pd.read_csv(processed_dir / "metadata_level1_with_embeddings.tsv", sep="\t")
    embeddings = np.load(embeddings_dir / "embeddings_level1.npy")
    umap_coords = np.load(embeddings_dir / "umap_coords_level1.npy")

    print(f"Loaded {len(embeddings)} embeddings, {len(metadata)} metadata")

    # Generate F1: Dataset Overview
    print("Generating F1: Dataset Overview...")
    mvp_figures.make_f1_dataset_overview(metadata, figures_dir)
    print("✓ F1 generated")

    # Generate F4: UMAP Old/New World
    print("Generating F4: UMAP Old/New World...")
    mvp_figures.make_f4_umap_oldnewworld(umap_coords, metadata, figures_dir)
    print("✓ F4 generated")

def main():
    print("="*70)
    print("UPDATING ALIGNED METADATA AND REGENERATING F1/F4 FIGURES")
    print("="*70)

    # Update aligned metadata
    aligned_metadata = update_aligned_metadata()

    # Regenerate figures
    regenerate_f1_f4_figures()

    print(f"\n✅ Updated aligned metadata and regenerated F1 and F4 figures")
    print(f"Figures most affected by species counts have been updated!")

    return 0

if __name__ == "__main__":
    exit(main())