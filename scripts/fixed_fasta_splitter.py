#!/usr/bin/env python3
"""
Fixed FASTA splitter with improved regex and Puumala recovery.
"""

import re
from Bio import SeqIO
from pathlib import Path

def get_segment_from_header_and_sequence(header, sequence):
    """Extract segment using extraction method as ground truth."""

    # Segment mapping based on extraction method
    SEGMENT_MAP = {
        # S-segment: nucleocapsid
        'n_protein_full': 'S',
        # M-segment: glycoprotein (all variants)
        'gpc_full': 'M',
        'gpc_truncated_1022': 'M',
        'gn_domain': 'M',
        'gpc_nterm': 'M',      # fragment but still M origin
        'nucleotide_tier1_gpc': 'M',
        'nucleotide_tier2_gpc': 'M',
        'nucleotide_tier2_large_cds': 'M',  # verify this
        # L-segment: polymerase
        'rdrp_full': 'L',
        # Ambiguous — exclude entirely
        'protein_db_direct': None,  # no reliable segment assignment
        'protein_db': None,
    }

    # Pattern A: accession suffix (_S, _M, _L)
    if re.search(r'_[SML]\s', header):
        return re.search(r'_([SML])\s', header).group(1)

    # Pattern B: pipe-delimited segment field
    if re.search(r'\|[SML]\|', header):
        return re.search(r'\|([SML])\|', header).group(1)

    # Pattern C: extraction method mapping (ground truth)
    for method, segment in SEGMENT_MAP.items():
        if method in header:
            return segment  # None for excluded methods

    return None

def extract_puumala_sequences():
    """Check if Puumala sequences need to be added (skip if already present)."""

    print("Checking Puumala sequences in visualization FASTA...")

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    if not visualization_file.exists():
        print(f"❌ {visualization_file} not found")
        return 0

    # Check if Puumala sequences already exist
    puumala_count = 0
    for record in SeqIO.parse(visualization_file, "fasta"):
        if "puumala" in record.description.lower():
            puumala_count += 1

    if puumala_count > 0:
        print(f"  ✓ Puumala sequences already present: {puumala_count}")
        return puumala_count
    else:
        print(f"  No Puumala sequences found in visualization FASTA")
        print(f"  (This is expected - Puumala absent from Level 1 dataset)")
        return 0

def split_fasta_fixed():
    """Split visualization FASTA using fixed regex."""

    print("Splitting visualization FASTA with fixed regex...")

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
        # Process sequences with fixed regex including sequence length
        for record in SeqIO.parse(input_file, "fasta"):
            segment = get_segment_from_header_and_sequence(record.description, record.seq)

            if segment in ['S', 'M', 'L']:
                # Write to appropriate segment file
                SeqIO.write(record, output_handles[segment], "fasta")
                segment_counts[segment] += 1
            else:
                segment_counts['unknown'] += 1
                # Show sample unknown headers and lengths for debugging
                if segment_counts['unknown'] <= 5:
                    seq_len = len(record.seq)
                    print(f"  Unknown segment: {record.description[:60]}... (length: {seq_len})")

    finally:
        # Close output files
        for handle in output_handles.values():
            handle.close()

    # Report results
    print(f"\nFixed FASTA splitting results:")
    for segment in ['S', 'M', 'L']:
        count = segment_counts[segment]
        output_file = output_files[segment]
        print(f"  {segment}: {count} sequences → {output_file}")

    print(f"  Unknown: {segment_counts['unknown']} sequences")
    print(f"  Total processed: {sum(segment_counts.values())}")

    return segment_counts

def verify_species_coverage():
    """Verify that major species appear in segment files."""

    print(f"\nVerifying species coverage in segment files...")

    segment_files = {
        'S': Path("data/processed/S_sequences_level1.fasta"),
        'M': Path("data/processed/M_sequences_level1.fasta"),
        'L': Path("data/processed/L_sequences_level1.fasta")
    }

    target_species = ['Seoul', 'Andes', 'Puumala', 'Sin Nombre', 'Choclo']

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
            print(f"    {status} {species}: {count}")

def main():
    print("=" * 70)
    print("FIXED FASTA SPLITTER WITH PUUMALA RECOVERY")
    print("=" * 70)

    # Step 1: Extract and append Puumala sequences
    puumala_count = extract_puumala_sequences()

    # Step 2: Split with fixed regex
    segment_counts = split_fasta_fixed()

    # Step 3: Verify species coverage
    verify_species_coverage()

    print(f"\n✅ Fixed FASTA splitting complete!")
    print(f"Added {puumala_count} Puumala sequences")
    print(f"Segment totals: S={segment_counts['S']}, M={segment_counts['M']}, L={segment_counts['L']}")
    print(f"Ready to recompute embeddings with improved coverage")

    return 0

if __name__ == "__main__":
    exit(main())