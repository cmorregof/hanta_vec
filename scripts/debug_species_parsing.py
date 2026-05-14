#!/usr/bin/env python3
"""
Debug species name parsing to see what's being extracted.
"""

from Bio import SeqIO
from pathlib import Path

def debug_species_parsing():
    """Debug what species names are being extracted."""

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    species_examples = {}
    count = 0

    for record in SeqIO.parse(visualization_file, "fasta"):
        description = record.description

        # Current parsing logic
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

        # Collect examples
        if species not in species_examples:
            species_examples[species] = []
        if len(species_examples[species]) < 3:
            species_examples[species].append(description[:80])

        count += 1
        if count > 100:  # Only process first 100 for debugging
            break

    print("Species names extracted and sample headers:")
    for species, examples in species_examples.items():
        print(f"\n'{species}':")
        for example in examples:
            print(f"  {example}")

if __name__ == "__main__":
    debug_species_parsing()