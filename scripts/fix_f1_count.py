#!/usr/bin/env python3
"""
Fix F1 dataset overview with correct N count and updated species colors.
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

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

UPDATED_CLADE_COLORS = {
    "Old World": "#2166ac",       # Blue
    "New World": "#d6604d",       # Red
    "Unknown": "#aaaaaa",         # Gray
}

def make_f1_fixed(metadata: pd.DataFrame, output_dir: Path):
    """Generate F1 with fixed N count and updated colors."""

    print("  Generating F1: Dataset Overview with correct N count...")

    # Get actual count for title
    n_sequences = len(metadata)

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # Panel A: Species barplot (horizontal)
    ax_a = fig.add_subplot(gs[0, 0])
    species_counts = metadata["species"].value_counts().sort_values()
    colors_list = [
        UPDATED_SPECIES_COLORS.get(s, UPDATED_SPECIES_COLORS["other"]) for s in species_counts.index
    ]
    ax_a.barh(range(len(species_counts)), species_counts.values, color=colors_list)
    ax_a.set_yticks(range(len(species_counts)))
    ax_a.set_yticklabels([s[:30] for s in species_counts.index], fontsize=9)
    ax_a.set_xlabel("Count")
    ax_a.set_title("A. Species Distribution", fontweight="bold", loc="left")

    print(f"    Panel A: {len(species_counts)} species")

    # Panel B: Sequence length histogram
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.hist(metadata["length"], bins=20, color="steelblue", edgecolor="black")
    median_len = metadata["length"].median()
    ax_b.axvline(median_len, color="red", linestyle="--", linewidth=2, label=f"Median: {median_len:.0f}")
    ax_b.set_xlabel("Sequence length (aa)")
    ax_b.set_ylabel("Count")
    ax_b.set_title("B. Sequence Length Distribution", fontweight="bold", loc="left")
    ax_b.legend()

    print(f"    Panel B: Length range {metadata['length'].min()}-{metadata['length'].max()}, median {median_len:.0f}")

    # Panel C: Old World vs New World (use "clade" column)
    ax_c = fig.add_subplot(gs[1, 0])
    clade_counts = metadata["clade"].value_counts()
    clade_order = ["Old World", "New World", "Unknown"]
    clade_counts = clade_counts.reindex([x for x in clade_order if x in clade_counts.index])
    colors_clade = [UPDATED_CLADE_COLORS[x] for x in clade_counts.index]
    ax_c.bar(range(len(clade_counts)), clade_counts.values, color=colors_clade, edgecolor="black")
    ax_c.set_xticks(range(len(clade_counts)))
    ax_c.set_xticklabels(clade_counts.index, rotation=15, ha="right")
    ax_c.set_ylabel("Count")
    ax_c.set_title("C. Old World vs New World", fontweight="bold", loc="left")

    print(f"    Panel C: Old World {clade_counts.get('Old World', 0)}, New World {clade_counts.get('New World', 0)}")

    # Panel D: Extraction methods
    ax_d = fig.add_subplot(gs[1, 1])
    method_counts = metadata["method"].value_counts()
    ax_d.bar(range(len(method_counts)), method_counts.values, color="skyblue", edgecolor="black")
    ax_d.set_xticks(range(len(method_counts)))
    ax_d.set_xticklabels([m[:20] for m in method_counts.index], rotation=45, ha="right", fontsize=9)
    ax_d.set_ylabel("Count")
    ax_d.set_title("D. Extraction Methods", fontweight="bold", loc="left")

    print(f"    Panel D: {len(method_counts)} extraction methods")

    # Fixed title with correct N count
    fig.suptitle(f"HantaVec Dataset Overview — Three Segments (N={n_sequences})",
                fontsize=14, fontweight="bold", y=0.995)

    output_path = output_dir / "small" / "F1_dataset_overview.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"  ✓ Saved {output_path}")
    plt.close(fig)

    print(f"  ✓ F1 regenerated with correct N={n_sequences}")

def main():
    print("="*70)
    print("FIXING F1 DATASET OVERVIEW WITH CORRECT N COUNT")
    print("="*70)

    # Load data
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    figures_dir = Path(__file__).parent.parent / "results" / "figures"

    # Load aligned metadata
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    if not metadata_path.exists():
        print("❌ Required metadata file not found!")
        return 1

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded {len(metadata)} sequences")

    # Generate fixed F1
    make_f1_fixed(metadata, figures_dir)

    print(f"\n✅ F1 fixed! Title now shows correct N={len(metadata)} instead of stale N=398")
    print(f"All species colors updated to match current dataset")

    return 0

if __name__ == "__main__":
    exit(main())