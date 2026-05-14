#!/usr/bin/env python3
"""
Integrate validated New World sequences into the main dataset.
Focus on Sin Nombre sequences with correct taxid 3052499.
"""

import json
import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def main():
    print("="*70)
    print("INTEGRATING NEW WORLD SEQUENCES INTO MAIN DATASET")
    print("="*70)

    # Load new world results
    new_world_file = Path(__file__).parent.parent / "data" / "raw" / "new_world_targeted" / "new_world_results.json"
    if not new_world_file.exists():
        print("❌ New World results file not found!")
        return 1

    with open(new_world_file) as f:
        new_world_results = json.load(f)

    print(f"New World sequences available:")
    for segment, sequences in new_world_results.items():
        print(f"  {segment}-segment: {len(sequences)} sequences")

    # Filter for validated species (focus on Sin Nombre for now)
    validated_species = ["Sin Nombre"]

    filtered_results = {"S": {}, "M": {}, "L": {}}
    species_counts = {}

    for segment, sequences in new_world_results.items():
        for acc_id, seq_data in sequences.items():
            species = seq_data.get('species', 'Unknown')

            if species in validated_species:
                filtered_results[segment][acc_id] = seq_data

                if species not in species_counts:
                    species_counts[species] = {"S": 0, "M": 0, "L": 0}
                species_counts[species][segment] += 1

    print(f"\nFiltered for validated species:")
    for species, counts in species_counts.items():
        total = sum(counts.values())
        print(f"  {species}: {total} total (S={counts['S']}, M={counts['M']}, L={counts['L']})")

    # Load existing metadata to check for duplicates
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    existing_metadata = pd.read_csv(processed_dir / "metadata_level1.tsv", sep="\t")
    existing_accessions = set(existing_metadata['genbank_id'].values)

    print(f"\nExisting dataset: {len(existing_accessions)} sequences")

    # Check for duplicates and add new sequences
    new_sequences = {"S": {}, "M": {}, "L": {}}
    duplicate_count = 0
    new_count = 0

    for segment, sequences in filtered_results.items():
        for acc_id, seq_data in sequences.items():
            if acc_id in existing_accessions:
                duplicate_count += 1
            else:
                new_sequences[segment][acc_id] = seq_data
                new_count += 1

    print(f"Duplicates with existing dataset: {duplicate_count}")
    print(f"New sequences to add: {new_count}")

    if new_count == 0:
        print("✓ No new sequences to add - all sequences already in dataset")
        return 0

    # Add new sequences to FASTA files
    for segment in ['S', 'M', 'L']:
        if not new_sequences[segment]:
            print(f"No new {segment}-segment sequences to add")
            continue

        # Load existing FASTA
        fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
        existing_records = list(SeqIO.parse(fasta_path, "fasta"))
        print(f"\n{segment}-segment: Adding {len(new_sequences[segment])} new sequences")

        # Create new records
        new_records = []
        for acc_id, seq_data in new_sequences[segment].items():
            record = SeqRecord(
                Seq(seq_data['sequence']),
                id=acc_id,
                description=f"{seq_data['species']} | {seq_data['isolate']} | {seq_data['method']}"
            )
            new_records.append(record)

        # Combine and save
        all_records = existing_records + new_records
        SeqIO.write(all_records, fasta_path, "fasta")
        print(f"  ✓ Updated {fasta_path} ({len(existing_records)} → {len(all_records)})")

    # Update metadata
    print(f"\nUpdating metadata...")
    new_metadata_rows = []

    for segment in ['S', 'M', 'L']:
        for acc_id, seq_data in new_sequences[segment].items():
            row = {
                'record_id': f"{acc_id}_{segment}",
                'genbank_id': acc_id,
                'segment': segment.upper(),
                'species': seq_data['species'],
                'taxon_id': seq_data.get('taxon_id', '3052499'),  # Sin Nombre taxid
                'clade': seq_data['clade'],
                'isolate': seq_data['isolate'],
                'strain': seq_data['strain'],
                'length': len(seq_data['sequence']),
                'method': seq_data['method'],
                'description': seq_data['description'][:100] + "..." if len(seq_data['description']) > 100 else seq_data['description']
            }
            new_metadata_rows.append(row)

    if new_metadata_rows:
        new_metadata_df = pd.DataFrame(new_metadata_rows)

        # Combine with existing
        combined_metadata = pd.concat([existing_metadata, new_metadata_df], ignore_index=True)
        combined_metadata.to_csv(processed_dir / "metadata_level1.tsv", sep="\t", index=False)

        print(f"  ✓ Updated metadata ({len(existing_metadata)} → {len(combined_metadata)} rows)")

        # Species distribution update
        species_dist = combined_metadata.groupby(['species', 'segment']).size().unstack(fill_value=0)
        print(f"\nUpdated species distribution:")
        for species in species_dist.index:
            total = species_dist.loc[species].sum()
            print(f"  {species}: {total} total (S={species_dist.loc[species, 'S']}, M={species_dist.loc[species, 'M']}, L={species_dist.loc[species, 'L']})")

    print(f"\n✅ Successfully integrated {new_count} new sequences")
    print(f"Dataset now contains: {len(combined_metadata)} total sequences")

    return 0

if __name__ == "__main__":
    sys.exit(main())