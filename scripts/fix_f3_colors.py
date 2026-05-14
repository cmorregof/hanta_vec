#!/usr/bin/env python3
"""
Fix F3 colors and regenerate with proper species color mapping.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from pathlib import Path
from scipy.spatial import ConvexHull

# Updated species color mapping for current dataset
UPDATED_SPECIES_COLORS = {
    "Seoul": "#e15759",           # Red
    "Andes": "#b07aa1",           # Purple
    "Puumala": "#4e79a7",         # Blue
    "Choclo": "#f28e2b",          # Orange
    "Sin Nombre": "#edc948",      # Yellow
    "Bayou": "#59a14f",           # Green
    "Caño Delgadito": "#76b7b2",  # Teal
    "Prospect Hill": "#af7aa1",   # Light purple
    "Black Creek Canal": "#9c755f", # Brown
    "Unknown": "#aaaaaa",         # Gray
    "other": "#cccccc",           # Light gray fallback
}

def make_f3_fixed(coords_umap: np.ndarray, metadata: pd.DataFrame, output_dir: Path):
    """Generate F3 with fixed colors, correct N count, and convex hulls."""

    print("  Generating F3: UMAP by Species with fixed colors...")

    # Get actual count for title
    n_sequences = len(metadata)

    print(f"  Species distribution:")
    species_counts = metadata['species'].value_counts()
    for species, count in species_counts.items():
        color = UPDATED_SPECIES_COLORS.get(species, UPDATED_SPECIES_COLORS["other"])
        print(f"    {species}: {count} (color: {color})")

    # Plotly (interactive)
    df_plot = metadata.copy()
    df_plot["umap_x"] = coords_umap[:, 0]
    df_plot["umap_y"] = coords_umap[:, 1]

    fig_plotly = px.scatter(
        df_plot,
        x="umap_x",
        y="umap_y",
        color="species",
        hover_data=["genbank_id", "clade", "segment_used", "isolate"],
        title=f"UMAP of ESM-2 Embeddings — Orthohantavirus Three Segments (N={n_sequences})",
        labels={"umap_x": "UMAP 1", "umap_y": "UMAP 2"},
        color_discrete_map=UPDATED_SPECIES_COLORS,
    )
    fig_plotly.update_traces(marker=dict(size=6))
    fig_plotly.update_layout(template="plotly_white", font=dict(size=11))

    html_path = output_dir / "large" / "F3_umap_species.html"
    html_path.parent.mkdir(parents=True, exist_ok=True)
    fig_plotly.write_html(html_path)
    print(f"  ✓ Saved interactive: {html_path}")

    # Static PNG with convex hulls
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot points by species with proper colors
    species_list = sorted(metadata["species"].unique())

    for species in species_list:
        mask = metadata["species"] == species
        species_coords = coords_umap[mask]

        if len(species_coords) == 0:
            continue

        color = UPDATED_SPECIES_COLORS.get(species, UPDATED_SPECIES_COLORS["other"])
        count = len(species_coords)

        # Plot points
        ax.scatter(
            species_coords[:, 0],
            species_coords[:, 1],
            c=color,
            label=f"{species} (n={count})",
            s=35,
            alpha=0.7,
            edgecolors="none",
            zorder=2
        )

        # Add convex hull for species with enough points
        if len(species_coords) >= 3 and species != "Unknown":
            try:
                hull = ConvexHull(species_coords)
                hull_points = species_coords[hull.vertices]
                hull_points = np.vstack([hull_points, hull_points[0]])  # Close the hull

                ax.plot(hull_points[:, 0], hull_points[:, 1],
                       color=color, linestyle='-', linewidth=1.5, alpha=0.8, zorder=1)
                ax.fill(hull_points[:, 0], hull_points[:, 1],
                       color=color, alpha=0.1, zorder=0)

            except Exception as e:
                print(f"    Could not create hull for {species}: {e}")

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(f"UMAP of ESM-2 Embeddings — Orthohantavirus Three Segments (N={n_sequences})",
                fontsize=13, fontweight="bold")

    # Legend outside plot area
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=9, framealpha=0.95)
    ax.grid(True, alpha=0.3)

    png_path = output_dir / "small" / "F3_umap_species.png"
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    print(f"  ✓ Saved static: {png_path}")
    plt.close(fig)

    print(f"  ✓ F3 colors fixed and regenerated with convex hulls")

def main():
    print("="*70)
    print("FIXING F3 COLORS AND REGENERATING")
    print("="*70)

    # Load data
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    figures_dir = Path(__file__).parent.parent / "results" / "figures"
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load aligned data
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"
    umap_coords_path = embeddings_dir / "umap_coords_level1.npy"

    if not all([p.exists() for p in [metadata_path, umap_coords_path]]):
        print("❌ Required files not found!")
        return 1

    metadata = pd.read_csv(metadata_path, sep="\t")
    umap_coords = np.load(umap_coords_path)

    print(f"Loaded {len(umap_coords)} UMAP coordinates, {len(metadata)} metadata")

    # Generate fixed F3
    make_f3_fixed(umap_coords, metadata, figures_dir)

    print(f"\n✅ F3 colors fixed! All species now have proper color mapping.")
    print(f"Species with convex hulls: {len([s for s in metadata['species'].value_counts().items() if s[1] >= 3 and s[0] != 'Unknown'])}")

    return 0

if __name__ == "__main__":
    exit(main())