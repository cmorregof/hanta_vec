#!/usr/bin/env python3
"""
Debug the accession mismatch between FASTA and embeddings.
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

def debug_accession_matching():
    """Debug accession extraction and matching logic."""

    print("=" * 70)
    print("DEBUG: ACCESSION MISMATCH ANALYSIS")
    print("=" * 70)

    # Load FASTA sequences
    print("\n📋 Analyzing FASTA accession extraction...")
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    fasta_by_segment = {'S': set(), 'M': set(), 'L': set()}
    extraction_examples = []

    count = 0
    for record in SeqIO.parse(visualization_file, "fasta"):
        segment = get_segment_from_extraction_method(record.description)
        if segment is None:
            continue

        base_accession = record.id.split('_')[0] if '_' in record.id else record.id
        fasta_by_segment[segment].add(base_accession)

        # Show examples of extraction logic
        if count < 10:
            extraction_examples.append({
                'record_id': record.id,
                'base_accession': base_accession,
                'segment': segment,
                'description': record.description[:60] + '...'
            })
        count += 1

    print(f"  Total FASTA sequences processed: {count}")
    for segment in ['S', 'M', 'L']:
        print(f"  {segment}: {len(fasta_by_segment[segment])} unique accessions")

    print(f"\n📝 Accession extraction examples:")
    for ex in extraction_examples:
        print(f"  {ex['record_id']} -> {ex['base_accession']} ({ex['segment']})")

    # Load embedding indices
    print(f"\n📋 Analyzing embedding indices...")
    embeddings_dir = Path("results/embeddings")

    embedded_by_segment = {'S': set(), 'M': set(), 'L': set()}

    for segment in ['S', 'M', 'L']:
        index_file = embeddings_dir / f"embeddings_{segment}_index.json"
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)

            if 'accession_ids' in index:
                embedded_by_segment[segment] = set(index['accession_ids'])
                print(f"  {segment}: {len(embedded_by_segment[segment])} embedded accessions")

                # Show examples
                sample_embedded = list(embedded_by_segment[segment])[:5]
                print(f"    Sample: {sample_embedded}")

    # Find overlaps and mismatches
    print(f"\n🔍 Overlap analysis:")
    total_overlap = 0
    total_fasta = 0
    total_embedded = 0

    for segment in ['S', 'M', 'L']:
        fasta_set = fasta_by_segment[segment]
        embedded_set = embedded_by_segment[segment]

        overlap = fasta_set.intersection(embedded_set)
        fasta_only = fasta_set - embedded_set
        embedded_only = embedded_set - fasta_set

        print(f"\n  {segment} segment:")
        print(f"    FASTA accessions: {len(fasta_set)}")
        print(f"    Embedded accessions: {len(embedded_set)}")
        print(f"    Overlap (matched): {len(overlap)}")
        print(f"    FASTA only (missing embeddings): {len(fasta_only)}")
        print(f"    Embedded only (orphaned embeddings): {len(embedded_only)}")

        if fasta_only:
            print(f"    Sample missing: {list(fasta_only)[:3]}")
        if embedded_only:
            print(f"    Sample orphaned: {list(embedded_only)[:3]}")

        total_overlap += len(overlap)
        total_fasta += len(fasta_set)
        total_embedded += len(embedded_set)

    overall_coverage = (total_overlap / total_fasta) * 100 if total_fasta > 0 else 0
    print(f"\n📊 Summary:")
    print(f"Total FASTA accessions: {total_fasta}")
    print(f"Total embedded accessions: {total_embedded}")
    print(f"Total matched: {total_overlap}")
    print(f"Overall coverage: {overall_coverage:.1f}%")

    return 0

if __name__ == "__main__":
    exit(debug_accession_matching())