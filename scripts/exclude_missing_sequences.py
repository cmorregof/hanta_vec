#!/usr/bin/env python3
"""
Exclude the 20 missing NC_ reference sequences to maintain dataset consistency.
"""

import json
import re
from Bio import SeqIO
from pathlib import Path

def get_segment_from_extraction_method(header):
    """Extract segment using extraction method as ground truth."""
    SEGMENT_MAP = {
        'n_protein_full': 'S',
        'gpc_full': 'M', 'gpc_truncated_1022': 'M', 'gn_domain': 'M', 'gpc_nterm': 'M',
        'nucleotide_tier1_gpc': 'M', 'nucleotide_tier2_gpc': 'M', 'nucleotide_tier2_large_cds': 'M',
        'rdrp_full': 'L',
    }

    if re.search(r'_[SML]\s', header):
        return re.search(r'_([SML])\s', header).group(1)

    for method, segment in SEGMENT_MAP.items():
        if method in header:
            return segment
    return None

def extract_base_accession(record_id):
    """Fixed accession extraction that handles periods correctly."""
    if '_' in record_id and record_id.split('_')[-1] in ['S', 'M', 'L']:
        base_with_version = record_id.rsplit('_', 1)[0]
    else:
        base_with_version = record_id
    return base_with_version.split()[0]

def extract_species_from_header(header):
    """Extract species name from header."""
    match = re.match(r'^[^\s]+\s+([^|]+)', header)
    if match:
        return match.group(1).strip()
    return 'unknown'

def main():
    print("=" * 70)
    print("EXCLUDING MISSING NC_ SEQUENCES FOR DATASET CONSISTENCY")
    print("=" * 70)

    # Load embedding indices to know which accessions we have
    embedded_accessions = {'S': set(), 'M': set(), 'L': set()}

    for segment in ['S', 'M', 'L']:
        index_file = Path(f"results/embeddings/embeddings_{segment}_index.json")
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)
            if 'accession_ids' in index:
                embedded_accessions[segment] = set(index['accession_ids'])

    print(f"📋 Loaded embedding indices:")
    for segment in ['S', 'M', 'L']:
        print(f"  {segment}: {len(embedded_accessions[segment])} embedded accessions")

    # Filter FASTA to only keep sequences with embeddings
    print(f"\n📋 Filtering FASTA sequences...")

    input_file = Path("data/processed/sequences_level1_for_visualization.fasta")
    output_file = Path("data/processed/sequences_level1_for_visualization_filtered.fasta")

    kept_sequences = []
    excluded_sequences = []

    for record in SeqIO.parse(input_file, "fasta"):
        segment = get_segment_from_extraction_method(record.description)
        if segment is None:
            excluded_sequences.append((record.id, "no_segment"))
            continue

        base_accession = extract_base_accession(record.id)

        if base_accession in embedded_accessions[segment]:
            kept_sequences.append(record)
        else:
            excluded_sequences.append((record.id, f"no_embedding_{segment}"))

    # Write filtered FASTA
    with open(output_file, 'w') as f:
        SeqIO.write(kept_sequences, f, "fasta")

    print(f"  Original sequences: {len(kept_sequences) + len(excluded_sequences)}")
    print(f"  Kept sequences: {len(kept_sequences)}")
    print(f"  Excluded sequences: {len(excluded_sequences)}")

    # Show excluded sequences
    if excluded_sequences:
        print(f"\n❌ Excluded sequences:")
        for seq_id, reason in excluded_sequences:
            print(f"  {seq_id} ({reason})")

    # Replace original file
    output_file.replace(input_file)
    print(f"\n✅ Updated {input_file}")

    # Generate final species counts
    print(f"\n📊 FINAL SPECIES COUNTS:")

    species_counts = {}

    for record in kept_sequences:
        segment = get_segment_from_extraction_method(record.description)
        species = extract_species_from_header(record.description)

        if species not in species_counts:
            species_counts[species] = {'S': 0, 'M': 0, 'L': 0}
        species_counts[species][segment] += 1

    target_species = ['Seoul', 'Andes', 'Puumala', 'Prospect Hill', 'Choclo', 'Sin Nombre']

    print(f"{'Species':<15} | {'S':<3} | {'M':<3} | {'L':<3} | {'Total':<5}")
    print("-" * 50)

    total_sequences = 0
    for species in target_species:
        if species in species_counts:
            counts = species_counts[species]
            s_count = counts['S']
            m_count = counts['M']
            l_count = counts['L']
            species_total = s_count + m_count + l_count

            print(f"{species:<15} | {s_count:<3} | {m_count:<3} | {l_count:<3} | {species_total:<5}")
            total_sequences += species_total

    # Other species
    other_total = 0
    for species, counts in species_counts.items():
        if species not in target_species:
            species_total = sum(counts.values())
            other_total += species_total

    if other_total > 0:
        print(f"{'Other species':<15} | {'':>3} | {'':>3} | {'':>3} | {other_total:<5}")
        total_sequences += other_total

    print("-" * 50)
    print(f"{'TOTAL':<15} | {'':>3} | {'':>3} | {'':>3} | {total_sequences:<5}")

    print(f"\n✅ Dataset consistency achieved!")
    print(f"Final dataset: {total_sequences} sequences")
    print(f"All sequences now have corresponding embeddings")

    return 0

if __name__ == "__main__":
    exit(main())