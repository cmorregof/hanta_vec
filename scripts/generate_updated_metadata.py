#!/usr/bin/env python3
"""
Generate updated metadata for filtered dataset (767 sequences, 100% coverage).
"""

import numpy as np
import pandas as pd
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

def extract_isolate_from_header(header):
    """Extract isolate from header."""
    # Format: >ACCESSION Species | isolate | method
    parts = header.split('|')
    if len(parts) >= 2:
        isolate = parts[1].strip()
        # Clean up isolate names
        if 'Old World' in isolate or 'New World' in isolate:
            return ''  # These are not real isolate names
        return isolate
    return ''

def classify_clade(species):
    """Classify species into Old World vs New World clades."""
    old_world = ['Seoul', 'Hantaan', 'Puumala', 'Dobrava']
    new_world = ['Andes', 'Sin Nombre', 'Choclo', 'Bayou', 'Caño Delgadito', 'Black Creek Canal', 'Muleshoe']

    if species in old_world:
        return 'Old World'
    elif species in new_world:
        return 'New World'
    else:
        return 'Unknown'

def main():
    print("=" * 70)
    print("GENERATING UPDATED METADATA FOR FILTERED DATASET")
    print("=" * 70)

    # Load filtered FASTA
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    # Load embedding indices to get the correct order
    embedding_indices = {}
    for segment in ['S', 'M', 'L']:
        index_file = Path(f"results/embeddings/embeddings_{segment}_index.json")
        with open(index_file, 'r') as f:
            index = json.load(f)
        embedding_indices[segment] = index['accession_ids']

    # Create metadata for each segment
    all_metadata = []
    embedding_idx = 0

    for segment in ['S', 'M', 'L']:
        print(f"\n📋 Processing {segment} segment...")

        segment_metadata = []
        segment_sequences = []

        # Get sequences for this segment
        for record in SeqIO.parse(visualization_file, "fasta"):
            record_segment = get_segment_from_extraction_method(record.description)
            if record_segment != segment:
                continue

            base_accession = extract_base_accession(record.id)
            species = extract_species_from_header(record.description)
            isolate = extract_isolate_from_header(record.description)
            clade = classify_clade(species)

            segment_sequences.append({
                'accession': base_accession,
                'record_id': record.id,
                'species': species,
                'isolate': isolate,
                'clade': clade,
                'sequence': str(record.seq)
            })

        # Sort by embedding index order
        accession_to_idx = {acc: i for i, acc in enumerate(embedding_indices[segment])}
        segment_sequences.sort(key=lambda x: accession_to_idx.get(x['accession'], 999))

        print(f"  Found {len(segment_sequences)} sequences")

        # Create metadata
        for i, seq_info in enumerate(segment_sequences):
            metadata_row = {
                'genbank_id': seq_info['accession'],
                'record_id': seq_info['record_id'],
                'species': seq_info['species'],
                'isolate': seq_info['isolate'],
                'clade': seq_info['clade'],
                'segment': segment,
                'segment_used': segment,
                'embedding_index': embedding_idx,
                'segment_index': i,
                'sequence_length': len(seq_info['sequence'])
            }

            segment_metadata.append(metadata_row)
            all_metadata.append(metadata_row)
            embedding_idx += 1

        # Save segment-specific metadata
        segment_df = pd.DataFrame(segment_metadata)
        segment_file = Path(f"results/embeddings/metadata_{segment}.csv")
        segment_df.to_csv(segment_file, index=False)
        print(f"  Saved: {segment_file}")

    # Save combined metadata
    all_df = pd.DataFrame(all_metadata)
    combined_file = Path("results/embeddings/metadata_combined.csv")
    all_df.to_csv(combined_file, index=False)
    print(f"\n✅ Saved combined metadata: {combined_file}")

    # Summary
    print(f"\n📊 METADATA SUMMARY:")
    print(f"Total sequences: {len(all_metadata)}")

    species_counts = all_df['species'].value_counts()
    for species, count in species_counts.head(8).items():
        print(f"  {species}: {count}")

    segment_counts = all_df['segment'].value_counts()
    print(f"\nSegment distribution:")
    for segment, count in segment_counts.items():
        print(f"  {segment}: {count}")

    print(f"\n✅ Metadata generation complete!")
    return 0

if __name__ == "__main__":
    exit(main())