#!/usr/bin/env python3
"""
Targeted NCBI fetch for 136 missing Puumala accessions by exact ID.
"""

import pandas as pd
import time
from Bio import Entrez, SeqIO
from pathlib import Path
from io import StringIO

# Set your email for NCBI
Entrez.email = "claude.research@anthropic.com"

def extract_missing_puumala_accessions():
    """Extract the 136 missing Puumala accession IDs from Level 1 metadata."""

    print("Extracting missing Puumala accession IDs...")

    # Load metadata
    metadata_path = Path("data/processed/metadata_level1.tsv")
    metadata = pd.read_csv(metadata_path, sep="\t")

    # Get Puumala entries
    puumala_metadata = metadata[metadata['species'].str.contains('Puumala', case=False, na=False)]
    print(f"  Total Puumala in metadata: {len(puumala_metadata)}")

    # Load current visualization FASTA to see what we already have
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")
    existing_accessions = set()

    for record in SeqIO.parse(visualization_file, "fasta"):
        if "puumala" in record.description.lower():
            # Extract base accession
            base_acc = record.id.split('_')[0]
            existing_accessions.add(base_acc)

    print(f"  Existing Puumala sequences: {len(existing_accessions)}")

    # Find missing accessions
    missing_accessions = []

    for _, row in puumala_metadata.iterrows():
        record_id = row['record_id']  # e.g., 'PV276153.1_S'
        base_accession = record_id.split('_')[0]  # Remove _S/_M/_L suffix

        if base_accession not in existing_accessions:
            missing_accessions.append({
                'base_accession': base_accession,
                'segment': row['segment'],
                'record_id': record_id,
                'isolate': row.get('isolate', ''),
                'method': row.get('method', ''),
                'description': row.get('description', '')
            })

    print(f"  Missing accessions to fetch: {len(missing_accessions)}")

    # Show segment breakdown
    segment_counts = {}
    for acc in missing_accessions:
        seg = acc['segment']
        segment_counts[seg] = segment_counts.get(seg, 0) + 1

    print("  Missing by segment:")
    for segment, count in segment_counts.items():
        print(f"    {segment}: {count} accessions")

    return missing_accessions

def fetch_sequences_from_ncbi(accession_list, batch_size=10):
    """Fetch sequences from NCBI in batches, trying protein then nucleotide."""

    print(f"\nFetching {len(accession_list)} sequences from NCBI...")

    fetched_sequences = {}
    failed_accessions = []
    batch_count = 0

    # Process in batches
    for i in range(0, len(accession_list), batch_size):
        batch = accession_list[i:i+batch_size]
        batch_accessions = [acc['base_accession'] for acc in batch]
        batch_count += 1

        print(f"  Batch {batch_count}: {batch_accessions[:3]}...")

        # Try protein database first
        success = False
        for db_name in ["protein", "nucleotide"]:
            try:
                print(f"    Trying {db_name} database...")

                # Fetch from NCBI
                handle = Entrez.efetch(
                    db=db_name,
                    id=batch_accessions,
                    rettype="fasta",
                    retmode="text"
                )

                fasta_data = handle.read()
                handle.close()

                if fasta_data.strip():
                    # Parse FASTA data with different parser for comments
                    fasta_io = StringIO(fasta_data)
                    try:
                        records = list(SeqIO.parse(fasta_io, "fasta"))
                    except:
                        # Try alternative parser for files with comments
                        fasta_io = StringIO(fasta_data)
                        records = list(SeqIO.parse(fasta_io, "fasta-2line"))

                    print(f"      Retrieved {len(records)} sequences from {db_name}")

                    # Map records back to accession info
                    for record in records:
                        # Extract base accession more carefully
                        record_acc = record.id.split('|')[0] if '|' in record.id else record.id

                        # Find matching accession info
                        for acc_info in batch:
                            if acc_info['base_accession'] == record_acc:
                                fetched_sequences[record_acc] = {
                                    'record': record,
                                    'info': acc_info
                                }
                                break

                    success = True
                    break

                else:
                    print(f"      No sequences returned from {db_name}")

            except Exception as e:
                print(f"      Error with {db_name}: {e}")
                continue

        if not success:
            print(f"    Failed to fetch batch from both databases")
            failed_accessions.extend(batch_accessions)

        # Be nice to NCBI servers
        time.sleep(2)

    print(f"\nFetch summary:")
    print(f"  Successfully fetched: {len(fetched_sequences)}")
    print(f"  Failed/not found: {len(failed_accessions)}")

    if failed_accessions:
        print(f"  Failed accessions: {failed_accessions[:5]}...")

    return fetched_sequences, failed_accessions

def format_and_append_sequences(fetched_sequences):
    """Format sequences with correct headers and append to visualization FASTA."""

    print(f"\nFormatting and appending {len(fetched_sequences)} sequences...")

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    # Get current count
    current_count = sum(1 for _ in SeqIO.parse(visualization_file, "fasta"))
    print(f"  Current visualization FASTA: {current_count} sequences")

    formatted_count = 0
    with open(visualization_file, "a") as f:
        for base_acc, seq_data in fetched_sequences.items():
            record = seq_data['record']
            info = seq_data['info']

            # Create header in visualization format
            # >ACCESSION_SEGMENT Puumala | isolate | extraction_method
            header = f"{info['record_id']} Puumala | {info['isolate']} | {info['method']}"

            # Update record
            record.id = info['record_id']
            record.description = header

            # Write to file
            SeqIO.write(record, f, "fasta")
            formatted_count += 1

    # Verify final count
    final_count = sum(1 for _ in SeqIO.parse(visualization_file, "fasta"))
    print(f"  Final visualization FASTA: {final_count} sequences")
    print(f"  Added: {final_count - current_count} sequences")

    return formatted_count

def main():
    print("=" * 70)
    print("TARGETED NCBI FETCH FOR MISSING PUUMALA SEQUENCES")
    print("=" * 70)

    # Step 1: Extract missing accession IDs
    missing_accessions = extract_missing_puumala_accessions()

    if not missing_accessions:
        print("No missing Puumala accessions to fetch")
        return 0

    # Step 2: Fetch from NCBI
    fetched_sequences, failed_accessions = fetch_sequences_from_ncbi(missing_accessions)

    if not fetched_sequences:
        print("No sequences successfully fetched from NCBI")
        return 1

    # Step 3: Format and append to visualization FASTA
    added_count = format_and_append_sequences(fetched_sequences)

    # Final report
    print(f"\n{'='*50}")
    print("TARGETED FETCH COMPLETE")
    print(f"{'='*50}")
    print(f"Requested: {len(missing_accessions)} Puumala sequences")
    print(f"Fetched successfully: {len(fetched_sequences)}")
    print(f"Failed/not found: {len(failed_accessions)}")
    print(f"Added to dataset: {added_count}")

    if failed_accessions:
        print(f"\nFailed accessions (first 10):")
        for acc in failed_accessions[:10]:
            print(f"  {acc}")

    success_rate = len(fetched_sequences) / len(missing_accessions) * 100
    print(f"\nSuccess rate: {success_rate:.1f}%")

    return 0

if __name__ == "__main__":
    exit(main())