#!/usr/bin/env python3
"""
Integrate the 101 new New World sequences into the main dataset.
"""

import json
import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():
    print("="*70)
    print("INTEGRATING 101 NEW WORLD SEQUENCES INTO MAIN DATASET")
    print("="*70)

    # Load new world results
    enhanced_file = Path(__file__).parent.parent / "data" / "raw" / "enhanced_new_world" / "enhanced_new_world_results.json"
    if not enhanced_file.exists():
        print("❌ Enhanced New World results file not found!")
        return 1

    with open(enhanced_file) as f:
        enhanced_results = json.load(f)

    print(f"Loaded enhanced results for {len(enhanced_results)} species")

    # Count total sequences
    total_new = 0
    for species, segments in enhanced_results.items():
        species_total = sum(len(seqs) for seqs in segments.values())
        total_new += species_total
        print(f"  {species}: {species_total} sequences")

    print(f"Total new sequences to integrate: {total_new}")

    # Load existing metadata
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1.tsv"
    existing_metadata = pd.read_csv(metadata_path, sep="\t")

    print(f"Existing dataset: {len(existing_metadata)} sequences")

    # Prepare new metadata rows
    new_metadata_rows = []

    for species, segments in enhanced_results.items():
        for segment, sequences in segments.items():
            for acc_id, seq_data in sequences.items():
                row = {
                    'record_id': f"{acc_id}_{segment}",
                    'genbank_id': acc_id,
                    'segment': segment.upper(),
                    'species': seq_data.get('species_classified', species),
                    'taxon_id': 'Unknown',  # Will need to map later
                    'clade': seq_data.get('clade', 'New World'),
                    'isolate': seq_data.get('isolate', ''),
                    'strain': seq_data.get('strain', ''),
                    'length': len(seq_data['sequence']),
                    'method': seq_data.get('method', 'unknown'),
                    'description': seq_data.get('description', '')[:100] + "..." if len(seq_data.get('description', '')) > 100 else seq_data.get('description', '')
                }
                new_metadata_rows.append(row)

    if new_metadata_rows:
        new_metadata_df = pd.DataFrame(new_metadata_rows)

        # Combine with existing
        combined_metadata = pd.concat([existing_metadata, new_metadata_df], ignore_index=True)

        print(f"Updated metadata: {len(existing_metadata)} → {len(combined_metadata)} (+{len(new_metadata_rows)})")

        # Save updated metadata
        combined_metadata.to_csv(metadata_path, sep="\t", index=False)
        print(f"✓ Updated metadata saved")

        # Update FASTA files
        for segment in ['S', 'M', 'L']:
            fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"

            # Load existing sequences
            existing_records = list(SeqIO.parse(fasta_path, "fasta")) if fasta_path.exists() else []

            # Add new sequences
            new_records_added = 0
            for species, segments in enhanced_results.items():
                if segment in segments:
                    for acc_id, seq_data in segments[segment].items():
                        record = SeqRecord(
                            Seq(seq_data['sequence']),
                            id=acc_id,
                            description=f"{species} | {seq_data.get('isolate', '')} | {seq_data.get('method', '')}"
                        )
                        existing_records.append(record)
                        new_records_added += 1

            if new_records_added > 0:
                # Save updated FASTA
                SeqIO.write(existing_records, fasta_path, "fasta")
                print(f"✓ Updated {segment}-segment FASTA (+{new_records_added} sequences)")

        # Generate final species report
        print(f"\n📊 FINAL SPECIES DISTRIBUTION:")
        species_counts = combined_metadata.groupby(['species', 'clade']).size().unstack(fill_value=0)
        total_counts = combined_metadata['species'].value_counts()

        for species in total_counts.index:
            count = total_counts[species]
            clade = combined_metadata[combined_metadata['species'] == species]['clade'].iloc[0]
            print(f"  {species} ({clade}): {count}")

        # Clade summary
        clade_summary = combined_metadata['clade'].value_counts()
        print(f"\n📊 CLADE SUMMARY:")
        for clade, count in clade_summary.items():
            pct = count/len(combined_metadata)*100
            print(f"  {clade}: {count} ({pct:.1f}%)")

        print(f"\n✅ Successfully integrated {total_new} new sequences!")
        print(f"Dataset expanded: {len(existing_metadata)} → {len(combined_metadata)} sequences")
        print(f"Ready for embedding computation and figure regeneration")

        return 0

    else:
        print("❌ No new metadata to integrate")
        return 1

if __name__ == "__main__":
    exit(main())