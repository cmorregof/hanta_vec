#!/usr/bin/env python3
"""
Split sequences_level1_for_visualization.fasta by segment tags into segment files.
"""

from Bio import SeqIO
from pathlib import Path
import re

def parse_segment_from_header(header):
    """Extract segment from header, handling both formats."""

    # Format 1: _S , _M , _L at end of header
    match = re.search(r'_([SML])\s', header)
    if match:
        return match.group(1)

    # Format 2: |S|, |M|, |L| in description
    match = re.search(r'\|([SML])\|', header)
    if match:
        return match.group(1)

    return None

def split_visualization_fasta():
    """Split visualization FASTA into segment-specific files."""

    print("Splitting sequences_level1_for_visualization.fasta by segment...")

    # Setup paths
    input_file = Path("data/processed/sequences_level1_for_visualization.fasta")
    output_dir = Path("data/processed")

    # Output files
    output_files = {
        'S': output_dir / "S_sequences_level1.fasta",
        'M': output_dir / "M_sequences_level1.fasta",
        'L': output_dir / "L_sequences_level1.fasta"
    }

    # Segment counters
    segment_counts = {'S': 0, 'M': 0, 'L': 0, 'unknown': 0}

    # Open output files
    output_handles = {}
    for segment in ['S', 'M', 'L']:
        output_handles[segment] = open(output_files[segment], 'w')

    try:
        # Process sequences
        for record in SeqIO.parse(input_file, "fasta"):
            segment = parse_segment_from_header(record.description)

            if segment in ['S', 'M', 'L']:
                # Write to appropriate segment file
                SeqIO.write(record, output_handles[segment], "fasta")
                segment_counts[segment] += 1
            else:
                segment_counts['unknown'] += 1
                print(f"  Unknown segment for: {record.id}")

    finally:
        # Close output files
        for handle in output_handles.values():
            handle.close()

    # Report results
    print(f"\nSegment file splitting complete:")
    for segment in ['S', 'M', 'L']:
        count = segment_counts[segment]
        output_file = output_files[segment]
        print(f"  {segment}: {count} sequences → {output_file}")

    print(f"  Unknown: {segment_counts['unknown']} sequences")
    print(f"  Total processed: {sum(segment_counts.values())}")

    return segment_counts

def verify_target_species():
    """Verify Seoul and Puumala appear in segment files."""

    print(f"\nVerifying Seoul and Puumala in segment files...")

    segment_files = {
        'S': Path("data/processed/S_sequences_level1.fasta"),
        'M': Path("data/processed/M_sequences_level1.fasta"),
        'L': Path("data/processed/L_sequences_level1.fasta")
    }

    target_species = ['Seoul', 'Puumala']

    for segment, file_path in segment_files.items():
        if not file_path.exists():
            print(f"  ❌ {segment}: File not found")
            continue

        species_found = {species: 0 for species in target_species}

        for record in SeqIO.parse(file_path, "fasta"):
            for species in target_species:
                if species.lower() in record.description.lower():
                    species_found[species] += 1

        print(f"  {segment} segment:")
        for species, count in species_found.items():
            status = "✓" if count > 0 else "❌"
            print(f"    {status} {species}: {count} sequences")

def main():
    print("=" * 60)
    print("SPLITTING VISUALIZATION FASTA BY SEGMENT")
    print("=" * 60)

    # Split the FASTA file
    segment_counts = split_visualization_fasta()

    # Verify target species
    verify_target_species()

    print(f"\n✅ FASTA splitting complete!")
    print(f"Ready to recompute embeddings for missing sequences.")

    return 0

if __name__ == "__main__":
    exit(main())