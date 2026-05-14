#!/usr/bin/env python3
"""
Recover missing Puumala sequences by extracting them from all_candidates.fasta
using accession IDs from Level 1 metadata.
"""

import pandas as pd
from Bio import SeqIO
from pathlib import Path

def extract_puumala_accessions():
    """Extract Puumala accession IDs from Level 1 metadata."""

    print("Extracting Puumala accession IDs from Level 1 metadata...")

    metadata_path = Path("data/processed/metadata_level1.tsv")
    metadata = pd.read_csv(metadata_path, sep="\t")

    # Find Puumala sequences
    puumala_metadata = metadata[metadata['species'].str.contains('Puumala', case=False, na=False)]
    print(f"  Found {len(puumala_metadata)} Puumala entries in metadata")

    # Extract accession IDs and segments
    accessions = []
    for _, row in puumala_metadata.iterrows():
        record_id = row['record_id']  # e.g., 'PV276153.1_S'
        base_accession = record_id.split('_')[0]  # Remove _S/_M/_L suffix
        segment = row['segment']

        accessions.append({
            'base_accession': base_accession,
            'segment': segment,
            'record_id': record_id,
            'isolate': row.get('isolate', ''),
            'extraction_method': row.get('method', '')
        })

    print(f"  Extracted {len(accessions)} unique accession IDs")

    # Show breakdown by segment
    segment_counts = {}
    for acc in accessions:
        seg = acc['segment']
        segment_counts[seg] = segment_counts.get(seg, 0) + 1

    for segment, count in segment_counts.items():
        print(f"    {segment}: {count} accessions")

    return accessions

def search_candidates_for_puumala(puumala_accessions):
    """Search all_candidates.fasta for Puumala accessions."""

    print("\nSearching all_candidates.fasta for Puumala sequences...")

    candidates_file = Path("data/raw/proteins/all_candidates.fasta")
    if not candidates_file.exists():
        print(f"❌ {candidates_file} not found")
        return []

    # Create lookup of base accessions
    accession_lookup = {acc['base_accession']: acc for acc in puumala_accessions}

    found_sequences = []
    missing_accessions = set(accession_lookup.keys())

    for record in SeqIO.parse(candidates_file, "fasta"):
        # Extract base accession from record ID
        if '.' in record.id:
            base_accession = record.id  # Keep full ID with version
        else:
            base_accession = record.id

        if base_accession in accession_lookup:
            acc_info = accession_lookup[base_accession]
            found_sequences.append({
                'record': record,
                'accession_info': acc_info
            })
            missing_accessions.discard(base_accession)

    print(f"  Found {len(found_sequences)} sequences in all_candidates.fasta")
    print(f"  Missing {len(missing_accessions)} sequences")

    if missing_accessions:
        print(f"  Missing accessions: {sorted(list(missing_accessions))[:5]}...")

    return found_sequences

def format_puumala_sequences(found_sequences):
    """Format Puumala sequences with correct headers."""

    print(f"\nFormatting {len(found_sequences)} Puumala sequences...")

    formatted_sequences = []

    for seq_info in found_sequences:
        record = seq_info['record']
        acc_info = seq_info['accession_info']

        # Create header in visualization format
        # >ACCESSION_SEGMENT Species | isolate | extraction_method
        header = f">{acc_info['record_id']} Puumala | {acc_info['isolate']} | {acc_info['extraction_method']}"

        # Update record
        record.id = acc_info['record_id']
        record.description = header[1:]  # Remove '>' prefix

        formatted_sequences.append(record)

    return formatted_sequences

def append_to_visualization_fasta(formatted_sequences):
    """Append formatted Puumala sequences to visualization FASTA."""

    print(f"\nAppending {len(formatted_sequences)} Puumala sequences to visualization FASTA...")

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    # Check current sequence count
    current_count = 0
    for record in SeqIO.parse(visualization_file, "fasta"):
        current_count += 1

    print(f"  Current visualization FASTA: {current_count} sequences")

    # Append Puumala sequences
    with open(visualization_file, "a") as f:
        for record in formatted_sequences:
            SeqIO.write(record, f, "fasta")

    # Verify final count
    final_count = 0
    for record in SeqIO.parse(visualization_file, "fasta"):
        final_count += 1

    print(f"  Final visualization FASTA: {final_count} sequences")
    print(f"  Added: {final_count - current_count} sequences")

def main():
    print("=" * 70)
    print("RECOVERING MISSING PUUMALA SEQUENCES")
    print("=" * 70)

    # Step 1: Extract Puumala accessions from metadata
    puumala_accessions = extract_puumala_accessions()

    if not puumala_accessions:
        print("No Puumala accessions found in metadata")
        return 0

    # Step 2: Search all_candidates.fasta
    found_sequences = search_candidates_for_puumala(puumala_accessions)

    if not found_sequences:
        print("No Puumala sequences found in all_candidates.fasta")
        return 1

    # Step 3: Format with correct headers
    formatted_sequences = format_puumala_sequences(found_sequences)

    # Step 4: Append to visualization FASTA
    append_to_visualization_fasta(formatted_sequences)

    print(f"\n✅ Puumala sequence recovery complete!")
    print(f"Recovered {len(found_sequences)} of {len(puumala_accessions)} expected sequences")

    return 0

if __name__ == "__main__":
    exit(main())