#!/usr/bin/env python3
"""
Re-fetch Andes virus using the correct taxid 1980456.
Replace the contaminated Andes sequences with genuine ones.
"""

import sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import (
    setup_ncbi, search_ncbi, fetch_batch, load_config,
    extract_s_segment_protein, extract_m_segment_protein, extract_l_segment_protein,
    extract_isolate_strain, classify_species_from_description
)

def simple_sequence_qc(sequence: str, min_length: int = 50) -> bool:
    """Simple QC check for sequence quality."""
    if len(sequence) < min_length:
        return False
    ambiguous_count = sum(1 for aa in sequence if aa in 'XBZ*')
    ambiguous_pct = ambiguous_count / len(sequence) * 100
    return ambiguous_pct <= 10

def main():
    print("="*70)
    print("RE-FETCH ANDES VIRUS WITH CORRECT TAXID 1980456")
    print("="*70)

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # Use correct Andes taxid
    correct_andes_taxid = "1980456"

    print(f"Fetching Andes virus with correct taxid {correct_andes_taxid}...")

    try:
        query = f"txid{correct_andes_taxid}[Organism:exp]"
        nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=1000)
        print(f"  Found {len(nucleotide_ids)} nucleotide records")

        if not nucleotide_ids:
            print("❌ No records found!")
            return 1

        # Download records
        print(f"  Downloading records...")
        andes_records = {}

        for records, failed in fetch_batch(
            nucleotide_ids,
            db="nucleotide",
            batch_size=config["ncbi"]["batch_size"],
            sleep_time=config["ncbi"]["sleep_between_batches"],
        ):
            for record in records:
                andes_records[record.id] = record

        print(f"  ✓ Downloaded {len(andes_records)} records")

        # Check organisms
        organism_counts = {}
        for record in andes_records.values():
            organism = record.annotations.get('organism', 'Unknown')
            organism_counts[organism] = organism_counts.get(organism, 0) + 1

        print(f"  Organisms found:")
        for org, count in sorted(organism_counts.items()):
            print(f"    {org}: {count}")

        # Extract segments
        extractors = {
            'S': extract_s_segment_protein,
            'M': extract_m_segment_protein,
            'L': extract_l_segment_protein
        }

        andes_results = {"S": [], "M": [], "L": []}
        total_extracted = 0

        for record_id, record in andes_records.items():
            species, clade = classify_species_from_description(record.description)

            for segment, extractor in extractors.items():
                protein_data = extractor(record)

                if protein_data:
                    min_length = {"S": 200, "M": 400, "L": 1000}.get(segment, 100)
                    qc_passed = simple_sequence_qc(protein_data['sequence'], min_length)

                    if qc_passed:
                        isolate, strain = extract_isolate_strain(record)

                        andes_data = {
                            'record_id': f"{record.id}_{segment}",
                            'genbank_id': record.id,
                            'segment': segment,
                            'species': species,
                            'taxon_id': correct_andes_taxid,
                            'clade': clade,
                            'isolate': isolate,
                            'strain': strain,
                            'length': len(protein_data['sequence']),
                            'method': protein_data['method'],
                            'description': record.description,
                            'sequence': protein_data['sequence']
                        }

                        andes_results[segment].append(andes_data)
                        total_extracted += 1

        print(f"\nExtraction results:")
        for segment in ['S', 'M', 'L']:
            count = len(andes_results[segment])
            print(f"  {segment}-segment: {count} sequences")

        print(f"  Total extracted: {total_extracted}")

        if total_extracted == 0:
            print("❌ No sequences extracted!")
            return 1

        # Remove old contaminated Andes sequences from dataset
        print(f"\nRemoving contaminated Andes sequences...")
        processed_dir = Path(__file__).parent.parent / "data" / "processed"
        metadata_path = processed_dir / "metadata_level1.tsv"

        metadata = pd.read_csv(metadata_path, sep="\t")
        print(f"  Before: {len(metadata)} total sequences")

        # Remove old Andes sequences (they have wrong organisms)
        metadata_clean = metadata[metadata['species'] != 'Andes'].copy()
        print(f"  Removed: {len(metadata) - len(metadata_clean)} contaminated Andes sequences")

        # Add new genuine Andes sequences
        new_andes_rows = []
        for segment in ['S', 'M', 'L']:
            for andes_data in andes_results[segment]:
                new_andes_rows.append({k: v for k, v in andes_data.items() if k != 'sequence'})

        if new_andes_rows:
            new_andes_df = pd.DataFrame(new_andes_rows)
            metadata_updated = pd.concat([metadata_clean, new_andes_df], ignore_index=True)

            print(f"  Added: {len(new_andes_rows)} genuine Andes sequences")
            print(f"  After: {len(metadata_updated)} total sequences")

            # Save updated metadata
            metadata_updated.to_csv(metadata_path, sep="\t", index=False)
            print(f"  ✓ Updated metadata saved")

            # Update FASTA files
            from Bio import SeqIO
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq

            for segment in ['S', 'M', 'L']:
                if andes_results[segment]:
                    fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"

                    # Load existing sequences (excluding old Andes)
                    existing_records = []
                    if fasta_path.exists():
                        for record in SeqIO.parse(fasta_path, "fasta"):
                            # Keep non-Andes sequences
                            if not any(old_andes_id in record.id for old_andes_id in [
                                r['genbank_id'] for r in metadata[metadata['species'] == 'Andes'].to_dict('records')
                            ]):
                                existing_records.append(record)

                    # Add new Andes sequences
                    for andes_data in andes_results[segment]:
                        record = SeqRecord(
                            Seq(andes_data['sequence']),
                            id=andes_data['genbank_id'],
                            description=f"Andes | {andes_data['isolate']} | {andes_data['method']}"
                        )
                        existing_records.append(record)

                    # Save updated FASTA
                    SeqIO.write(existing_records, fasta_path, "fasta")
                    print(f"  ✓ Updated {segment}-segment FASTA")

        print(f"\n✅ Successfully replaced contaminated Andes sequences with {total_extracted} genuine ones")
        print(f"Need to re-run embedding computation and figure generation")

        return 0

    except Exception as e:
        print(f"ERROR: {e}")
        return 1

if __name__ == "__main__":
    exit(main())