#!/usr/bin/env python3
"""
Restore the original 893-sequence visualization FASTA file.
"""

from Bio import SeqIO
from pathlib import Path

def restore_original():
    """Extract first 893 sequences from the inflated backup file."""

    backup_file = Path("data/processed/sequences_level1_for_visualization.fasta.backup")
    output_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    print(f"Restoring original visualization FASTA...")
    print(f"Reading from: {backup_file}")

    count = 0
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(backup_file, "fasta"):
            if count >= 893:
                break

            SeqIO.write(record, out, "fasta")
            count += 1

            if count % 100 == 0:
                print(f"  Restored {count} sequences...")

    print(f"✓ Restored {count} sequences to {output_file}")
    return count

def main():
    count = restore_original()
    print(f"Original visualization FASTA restored with {count} sequences")
    return 0

if __name__ == "__main__":
    exit(main())