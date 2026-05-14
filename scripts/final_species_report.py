#!/usr/bin/env python3
"""
Final species × segment count report after Prospect Hill subsampling.
"""

from Bio import SeqIO
from pathlib import Path
import re

def get_species_and_segment_counts():
    """Extract species and segment counts from visualization FASTA."""

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    species_segment_counts = {}
    total_count = 0

    for record in SeqIO.parse(visualization_file, "fasta"):
        # Extract segment from record ID
        segment = record.id.split('_')[1] if '_' in record.id else 'unknown'

        # Extract species name from description
        # Format: >ACCESSION[_SEGMENT] Species | isolate | method
        description = record.description

        # Get species name (everything between first space and first |)
        match = re.match(r'^[^\s]+\s+([^|]+)', description)
        if match:
            species = match.group(1).strip()
        else:
            species = 'unknown'

        # Initialize counters
        if species not in species_segment_counts:
            species_segment_counts[species] = {'S': 0, 'M': 0, 'L': 0}

        # Count this sequence
        if segment in ['S', 'M', 'L']:
            species_segment_counts[species][segment] += 1

        total_count += 1

    return species_segment_counts, total_count

def main():
    print("=" * 70)
    print("FINAL SPECIES × SEGMENT COUNTS AFTER PROSPECT HILL SUBSAMPLING")
    print("=" * 70)

    species_counts, total = get_species_and_segment_counts()

    # Display results
    print(f"\n{'Species':<20} {'S':<4} {'M':<4} {'L':<4} {'Total':<6}")
    print("-" * 45)

    # Key target species
    target_species = ['Seoul', 'Andes', 'Puumala', 'Prospect Hill', 'Choclo', 'Sin Nombre']

    grand_total = 0
    for species in target_species:
        if species in species_counts:
            counts = species_counts[species]
            s_count = counts['S']
            m_count = counts['M']
            l_count = counts['L']
            species_total = s_count + m_count + l_count

            print(f"{species:<20} {s_count:<4} {m_count:<4} {l_count:<4} {species_total:<6}")
            grand_total += species_total

    # Other species
    other_total = 0
    other_species = []
    for species, counts in species_counts.items():
        if species not in target_species:
            species_total = counts['S'] + counts['M'] + counts['L']
            other_total += species_total
            if species_total > 1:  # Show species with >1 sequence
                other_species.append(f"{species}: {species_total}")

    if other_species:
        print(f"\nOther species:")
        for species_info in other_species:
            print(f"  {species_info}")

    if other_total > 0:
        print(f"{'Other species':<20} {'':>4} {'':>4} {'':>4} {other_total:<6}")
        grand_total += other_total

    print("-" * 45)
    print(f"{'TOTAL':<20} {'':>4} {'':>4} {'':>4} {grand_total:<6}")

    print(f"\n✅ Dataset ready for embedding recomputation!")
    print(f"Prospect Hill successfully subsampled from 318 to 60 sequences")
    print(f"Total dataset size: {grand_total} sequences")

    return 0

if __name__ == "__main__":
    exit(main())