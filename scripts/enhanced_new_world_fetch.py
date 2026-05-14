#!/usr/bin/env python3
"""
Enhanced New World species fetch using direct NCBI searches.
Since ViPR API is unavailable, use comprehensive NCBI organism searches.
"""

import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from typing import Dict, List
import json

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import (
    setup_ncbi, search_ncbi, fetch_batch, load_config,
    extract_s_segment_protein, extract_m_segment_protein, extract_l_segment_protein,
    extract_isolate_strain, classify_species_from_description
)

def sequence_identity(seq1: str, seq2: str) -> float:
    """Compute simple sequence identity."""
    if not seq1 or not seq2:
        return 0.0

    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0

    matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
    return matches / min_len

def deduplicate_against_existing(new_sequences: Dict, existing_sequences: Dict, threshold: float = 0.99) -> Dict:
    """Deduplicate new sequences against existing dataset."""
    print(f"  Deduplicating against {len(existing_sequences)} existing sequences (threshold={threshold})...")

    deduplicated = {}
    duplicate_count = 0

    for new_acc, new_data in new_sequences.items():
        new_seq = new_data.get('sequence', '')
        is_duplicate = False

        for existing_acc, existing_seq in existing_sequences.items():
            identity = sequence_identity(new_seq, existing_seq)

            if identity >= threshold:
                is_duplicate = True
                duplicate_count += 1
                print(f"    Duplicate: {new_acc} vs {existing_acc} ({identity:.3f})")
                break

        if not is_duplicate:
            deduplicated[new_acc] = new_data

    print(f"  Removed {duplicate_count} duplicates, kept {len(deduplicated)} unique")
    return deduplicated

def simple_qc(sequence: str, segment: str) -> bool:
    """Simple QC for extracted sequences."""
    min_lengths = {"S": 200, "M": 400, "L": 1000}
    min_len = min_lengths.get(segment, 100)

    if len(sequence) < min_len:
        return False

    ambiguous_count = sum(1 for aa in sequence if aa in 'XBZ*')
    ambiguous_pct = ambiguous_count / len(sequence) * 100

    return ambiguous_pct <= 10

def main():
    print("="*70)
    print("ENHANCED NEW WORLD HANTAVIRUS FETCH (NCBI DIRECT)")
    print("="*70)

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # Target New World species with various search terms
    target_species = {
        "Black Creek Canal": [
            "Black Creek Canal virus",
            "Black Creek Canal hantavirus",
            "Orthohantavirus nigrorivense"
        ],
        "Caño Delgadito": [
            "Caño Delgadito virus",
            "Cano Delgadito virus",
            "Caño Delgadito hantavirus"
        ],
        "Choclo": [
            "Choclo virus",
            "Choclo hantavirus",
            "Orthohantavirus chocloense"
        ],
        "Bayou": [
            "Bayou virus",
            "Bayou hantavirus",
            "Orthohantavirus bayoui"
        ],
        "Muleshoe": [
            "Muleshoe virus",
            "Muleshoe hantavirus"
        ]
    }

    # Load existing sequences for deduplication
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    existing_sequences = {}

    for segment in ['S', 'M', 'L']:
        fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
        if fasta_path.exists():
            for record in SeqIO.parse(fasta_path, "fasta"):
                existing_sequences[record.id] = str(record.seq)

    print(f"Loaded {len(existing_sequences)} existing sequences for deduplication")

    # Results collection
    all_results = {}
    total_new_sequences = 0

    extractors = {
        'S': extract_s_segment_protein,
        'M': extract_m_segment_protein,
        'L': extract_l_segment_protein
    }

    for species, search_terms in target_species.items():
        print(f"\n{'-'*50}")
        print(f"SPECIES: {species}")
        print(f"{'-'*50}")

        species_sequences = {"S": {}, "M": {}, "L": {}}
        all_records = {}

        for search_term in search_terms:
            print(f"  Searching: '{search_term}'")

            try:
                # Try organism-based search
                query = f'"{search_term}"'
                nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=500)
                print(f"    Found {len(nucleotide_ids)} nucleotide records")

                if nucleotide_ids:
                    # Download records
                    for records, failed in fetch_batch(
                        nucleotide_ids,
                        db="nucleotide",
                        batch_size=min(100, len(nucleotide_ids)),
                        sleep_time=1
                    ):
                        for record in records:
                            if record.id not in all_records:
                                all_records[record.id] = record

            except Exception as e:
                print(f"    Error: {e}")

        print(f"  Downloaded {len(all_records)} unique records")

        # Extract segments from all records
        extracted_count = 0

        for record_id, record in all_records.items():
            # Classify species
            classified_species, clade = classify_species_from_description(record.description)

            for segment, extractor in extractors.items():
                protein_data = extractor(record)

                if protein_data:
                    sequence = protein_data['sequence']

                    if simple_qc(sequence, segment):
                        isolate, strain = extract_isolate_strain(record)

                        sequence_data = {
                            'sequence': sequence,
                            'description': record.description,
                            'organism': record.annotations.get('organism', 'Unknown'),
                            'species_classified': classified_species,
                            'clade': clade,
                            'isolate': isolate,
                            'strain': strain,
                            'method': protein_data['method'],
                            'original_length': protein_data['original_length'],
                            'search_term': search_term
                        }

                        species_sequences[segment][record.id] = sequence_data
                        extracted_count += 1

        print(f"  Extracted {extracted_count} total sequences:")
        for segment in ['S', 'M', 'L']:
            count = len(species_sequences[segment])
            print(f"    {segment}-segment: {count}")

        # Deduplicate each segment
        deduplicated_sequences = {"S": {}, "M": {}, "L": {}}

        for segment in ['S', 'M', 'L']:
            if species_sequences[segment]:
                # Create map of existing sequences for this segment
                segment_existing = {}
                for acc, seq in existing_sequences.items():
                    if f"_{segment}" in acc or segment in acc:
                        segment_existing[acc] = seq

                # Deduplicate
                deduplicated = deduplicate_against_existing(
                    species_sequences[segment],
                    segment_existing,
                    threshold=0.99
                )

                deduplicated_sequences[segment] = deduplicated

        # Count new unique sequences
        species_new_count = sum(len(seqs) for seqs in deduplicated_sequences.values())

        if species_new_count > 0:
            all_results[species] = deduplicated_sequences
            total_new_sequences += species_new_count
            print(f"  ✅ {species_new_count} NEW unique sequences after deduplication")
        else:
            print(f"  ❌ No new unique sequences after deduplication")

    print(f"\n{'='*70}")
    print("ENHANCED NEW WORLD FETCH SUMMARY")
    print(f"{'='*70}")
    print(f"Total new sequences: {total_new_sequences}")

    if all_results:
        print(f"\nPer-species breakdown:")
        for species, segments in all_results.items():
            s_count = len(segments['S'])
            m_count = len(segments['M'])
            l_count = len(segments['L'])
            total = s_count + m_count + l_count
            print(f"  {species}: {total} (S:{s_count}, M:{m_count}, L:{l_count})")

        # Save results for integration
        output_dir = Path(__file__).parent.parent / "data" / "raw" / "enhanced_new_world"
        output_dir.mkdir(parents=True, exist_ok=True)

        with open(output_dir / "enhanced_new_world_results.json", "w") as f:
            json.dump(all_results, f, indent=2)

        print(f"\n✓ Results saved to {output_dir}")

        if total_new_sequences >= 5:
            print(f"✅ SUCCESS: Found {total_new_sequences} new sequences (≥5)")
            return 0
        else:
            print(f"⚠ PARTIAL: Found {total_new_sequences} new sequences (<5)")
            return 1
    else:
        print("❌ No new sequences found")
        return 1

if __name__ == "__main__":
    exit(main())