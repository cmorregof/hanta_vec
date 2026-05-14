#!/usr/bin/env python3
"""
Re-classify Unknown sequences using updated organism patterns.
Add patterns for Choclo, Caño Delgadito, Bayou, Black Creek Canal.
"""

import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from data.fetch import classify_species_from_description

def main():
    print("="*70)
    print("RE-CLASSIFYING UNKNOWN SEQUENCES WITH UPDATED PATTERNS")
    print("="*70)

    # Load current metadata
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1.tsv"

    if not metadata_path.exists():
        print("❌ Metadata file not found!")
        return 1

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded metadata: {len(metadata)} sequences")

    # Count current Unknown sequences
    unknown_before = len(metadata[metadata['species'] == 'Unknown'])
    print(f"Unknown sequences before re-classification: {unknown_before}")

    if unknown_before == 0:
        print("✓ No Unknown sequences to reclassify")
        return 0

    # Show current species distribution
    print(f"\nCurrent species distribution:")
    species_counts = metadata['species'].value_counts()
    for species, count in species_counts.items():
        print(f"  {species}: {count}")

    # Re-classify all sequences based on description
    print(f"\nRe-classifying sequences...")
    reclassified_count = 0

    for idx, row in metadata.iterrows():
        description = row.get('description', '')
        current_species = row.get('species', 'Unknown')

        if description:
            new_species, new_clade = classify_species_from_description(description)

            # Update if classification changed
            if new_species != current_species:
                metadata.at[idx, 'species'] = new_species
                metadata.at[idx, 'clade'] = new_clade
                reclassified_count += 1

                if current_species == 'Unknown':
                    print(f"  {row['genbank_id']}: Unknown → {new_species}")

    print(f"\nReclassified {reclassified_count} sequences")

    # Count remaining Unknown sequences
    unknown_after = len(metadata[metadata['species'] == 'Unknown'])
    print(f"Unknown sequences after re-classification: {unknown_after}")
    print(f"Improvement: {unknown_before - unknown_after} sequences classified")

    # Show updated species distribution
    print(f"\nUpdated species distribution:")
    updated_species_counts = metadata['species'].value_counts()
    for species, count in updated_species_counts.items():
        print(f"  {species}: {count}")

    # Show clade distribution
    print(f"\nUpdated clade distribution:")
    clade_counts = metadata['clade'].value_counts()
    total = len(metadata)
    for clade, count in clade_counts.items():
        pct = count/total*100
        print(f"  {clade}: {count} ({pct:.1f}%)")

    # Save updated metadata
    if reclassified_count > 0:
        metadata.to_csv(metadata_path, sep="\t", index=False)
        print(f"\n✓ Updated metadata saved to {metadata_path}")
    else:
        print(f"\n✓ No changes needed")

    # Check improvement in New World representation
    new_world_count = len(metadata[metadata['clade'] == 'New World'])
    new_world_pct = new_world_count / total * 100
    print(f"\nNew World representation: {new_world_count}/{total} ({new_world_pct:.1f}%)")

    return 0

if __name__ == "__main__":
    exit(main())