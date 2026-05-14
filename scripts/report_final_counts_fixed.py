#!/usr/bin/env python3
"""
Fixed species counting for final subsampling report.
"""

from Bio import SeqIO
from pathlib import Path

def report_final_species_counts_fixed():
    """Report exact final counts per species per segment with fixed parsing."""

    print("Reporting final species × segment counts (FIXED)...")

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    # Count by species and segment
    species_segment_counts = {}

    for record in SeqIO.parse(visualization_file, "fasta"):
        # Extract segment from record ID
        segment = record.id.split('_')[1] if '_' in record.id else 'unknown'

        # Extract species name - it's between accession and first |
        description = record.description

        # Format: >ACCESSION Species | isolate | method
        if '|' in description:
            # Split by | and get the part before first |
            before_pipe = description.split('|')[0].strip()
            # Remove accession to get species (everything after first space)
            if ' ' in before_pipe:
                species = ' '.join(before_pipe.split()[1:]).strip()
            else:
                species = 'unknown'
        else:
            species = 'unknown'

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
    other_details = []
    for species, counts in species_segment_counts.items():
        if species not in target_species and species != 'unknown':
            species_total = sum(counts.values())
            other_total += species_total
            if species_total > 3:  # Show species with >3 sequences
                s_count = counts.get('S', 0)
                m_count = counts.get('M', 0)
                l_count = counts.get('L', 0)
                other_details.append(f"{species}: {species_total} ({s_count}S+{m_count}M+{l_count}L)")

    # Show other species breakdown
    if other_details:
        print("\nOther significant species:")
        for detail in other_details:
            print(f"  {detail}")

    if other_total > 0:
        print(f"{'Other species':<20} {'':><4} {'':><4} {'':><4} {other_total:<6}")
        total_sequences += other_total

    print("-" * 45)
    print(f"{'TOTAL':<20} {'':><4} {'':><4} {'':><4} {total_sequences:<6}")

    return total_sequences, species_segment_counts

def main():
    total_sequences, species_counts = report_final_species_counts_fixed()

    print(f"\n✅ Final count verification complete!")
    print(f"Total sequences after Prospect Hill subsampling: {total_sequences}")

    return 0

if __name__ == "__main__":
    exit(main())