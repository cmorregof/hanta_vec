#!/usr/bin/env python3
"""
Subsample Prospect Hill to 60 sequences (20 per segment) using random sampling with seed=42.
"""

import random
import pandas as pd
from Bio import SeqIO
from pathlib import Path

def subsample_prospect_hill():
    """Subsample Prospect Hill sequences to 20 per segment."""

    print("Subsampling Prospect Hill sequences...")

    # Set random seed for reproducibility
    random.seed(42)

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")
    temp_file = Path("data/processed/sequences_level1_for_visualization_temp.fasta")

    # Collect sequences by category
    prospect_hill_sequences = {'S': [], 'M': [], 'L': []}
    other_sequences = []

    print("  Categorizing sequences...")

    for record in SeqIO.parse(visualization_file, "fasta"):
        if "prospect hill" in record.description.lower():
            # Extract segment from record ID
            segment = record.id.split('_')[1] if '_' in record.id else 'unknown'

            if segment in ['S', 'M', 'L']:
                prospect_hill_sequences[segment].append(record)
            else:
                print(f"    Warning: Prospect Hill sequence with unknown segment: {record.id}")
                other_sequences.append(record)
        else:
            other_sequences.append(record)

    # Report current Prospect Hill counts
    print(f"  Current Prospect Hill sequences:")
    total_ph = 0
    for segment, seqs in prospect_hill_sequences.items():
        count = len(seqs)
        total_ph += count
        print(f"    {segment}: {count} sequences")
    print(f"    Total: {total_ph} sequences")

    # Subsample to 20 per segment
    print(f"  Subsampling to 20 per segment...")
    subsampled_ph = {}
    final_ph_count = 0

    for segment, seqs in prospect_hill_sequences.items():
        if len(seqs) <= 20:
            # Keep all if 20 or fewer
            subsampled_ph[segment] = seqs
            print(f"    {segment}: kept all {len(seqs)} sequences")
        else:
            # Random sample 20
            subsampled = random.sample(seqs, 20)
            subsampled_ph[segment] = subsampled
            print(f"    {segment}: sampled {len(subsampled)} from {len(seqs)} sequences")

        final_ph_count += len(subsampled_ph[segment])

    print(f"  Final Prospect Hill count: {final_ph_count} sequences")

    # Write subsampled dataset
    print(f"  Writing subsampled dataset...")

    with open(temp_file, 'w') as f:
        # Write non-Prospect Hill sequences
        for record in other_sequences:
            SeqIO.write(record, f, "fasta")

        # Write subsampled Prospect Hill sequences
        for segment, seqs in subsampled_ph.items():
            for record in seqs:
                SeqIO.write(record, f, "fasta")

    # Replace original file
    temp_file.replace(visualization_file)

    # Verify final counts
    final_total = sum(1 for _ in SeqIO.parse(visualization_file, "fasta"))
    print(f"  Final dataset size: {final_total} sequences")

    return final_total

def report_final_species_counts():
    """Report exact final counts per species per segment."""

    print("\nReporting final species × segment counts...")

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    # Count by species and segment
    species_segment_counts = {}

    for record in SeqIO.parse(visualization_file, "fasta"):
        # Extract species and segment
        description = record.description
        segment = record.id.split('_')[1] if '_' in record.id else 'unknown'

        # Extract species name (between accession and first |)
        parts = description.split('|')
        if len(parts) >= 2:
            species_part = parts[0].strip()
            # Remove accession to get species
            species = ' '.join(species_part.split()[1:])  # Skip accession, take rest
        else:
            species = 'unknown'

        # Clean up species names
        species = species.strip()

        if species not in species_segment_counts:
            species_segment_counts[species] = {}

        if segment not in species_segment_counts[species]:
            species_segment_counts[species][segment] = 0

        species_segment_counts[species][segment] += 1

    # Report in organized format
    print(f"\nFinal species × segment table:")
    print(f"{'Species':<20} {'S':<4} {'M':<4} {'L':<4} {'Total':<6}")
    print("-" * 45)

    total_sequences = 0
    target_species = ['Seoul', 'Andes', 'Puumala', 'Prospect Hill', 'Choclo', 'Sin Nombre']

    for species in target_species:
        if species in species_segment_counts:
            counts = species_segment_counts[species]
            s_count = counts.get('S', 0)
            m_count = counts.get('M', 0)
            l_count = counts.get('L', 0)
            species_total = s_count + m_count + l_count

            print(f"{species:<20} {s_count:<4} {m_count:<4} {l_count:<4} {species_total:<6}")
            total_sequences += species_total

    # Report other species
    other_total = 0
    for species, counts in species_segment_counts.items():
        if species not in target_species and species != 'unknown':
            species_total = sum(counts.values())
            other_total += species_total
            if species_total > 5:  # Only show species with >5 sequences
                s_count = counts.get('S', 0)
                m_count = counts.get('M', 0)
                l_count = counts.get('L', 0)
                print(f"{species:<20} {s_count:<4} {m_count:<4} {l_count:<4} {species_total:<6}")

    if other_total > 0:
        print(f"{'Other species':<20} {'':+<4} {'':+<4} {'':+<4} {other_total:<6}")
        total_sequences += other_total

    print("-" * 45)
    print(f"{'TOTAL':<20} {'':+<4} {'':+<4} {'':+<4} {total_sequences:<6}")

    return total_sequences

def main():
    print("=" * 70)
    print("SUBSAMPLING PROSPECT HILL TO PRESERVE SCIENTIFIC VALUE")
    print("=" * 70)

    # Subsample Prospect Hill
    final_size = subsample_prospect_hill()

    # Report final counts
    total_sequences = report_final_species_counts()

    print(f"\n✅ Prospect Hill subsampling complete!")
    print(f"Dataset balanced for embedding recomputation: {total_sequences} sequences")
    print(f"Ready for fixed embedding pipeline (no zero padding)")

    return 0

if __name__ == "__main__":
    exit(main())