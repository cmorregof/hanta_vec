#!/usr/bin/env python3
"""
Debug accession matching between FASTA and embedding indices.
"""

import json
from Bio import SeqIO
from pathlib import Path

def debug_accession_matching():
    """Debug why accessions aren't matching between FASTA and embeddings."""

    print("=" * 70)
    print("DEBUG: ACCESSION MATCHING")
    print("=" * 70)

    # Check visualization FASTA accessions
    print("\n📋 VISUALIZATION FASTA SAMPLE ACCESSIONS:")
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")
    viz_accessions = []

    count = 0
    for record in SeqIO.parse(visualization_file, "fasta"):
        # Extract base accession
        base_acc = record.id.split('_')[0] if '_' in record.id else record.id
        viz_accessions.append(base_acc)

        if count < 10:
            print(f"  {record.id} -> {base_acc}")
        count += 1
        if count >= 20:  # Limit sample
            break

    print(f"Total visualization FASTA accessions (sample): {len(viz_accessions)}")

    # Check S segment FASTA accessions
    print(f"\n📋 S SEGMENT FASTA SAMPLE ACCESSIONS:")
    s_file = Path("data/processed/S_sequences_level1.fasta")
    s_accessions = []

    count = 0
    for record in SeqIO.parse(s_file, "fasta"):
        base_acc = record.id.split('_')[0] if '_' in record.id else record.id
        s_accessions.append(base_acc)

        if count < 10:
            print(f"  {record.id} -> {base_acc}")
        count += 1
        if count >= 20:
            break

    print(f"Total S segment FASTA accessions (sample): {len(s_accessions)}")

    # Check S embedding index
    print(f"\n📋 S EMBEDDING INDEX ACCESSIONS:")
    index_file = Path("results/embeddings/embeddings_S_final_index.json")

    with open(index_file, 'r') as f:
        index = json.load(f)

    # Get accession list from index
    if 'accession_ids' in index:
        embedded_accs = index['accession_ids'][:20]  # Sample
        print(f"  Sample embedded accessions: {embedded_accs[:10]}")
        print(f"Total embedded accessions: {len(index['accession_ids'])}")
    else:
        # Check if accessions are stored as keys
        accession_keys = [k for k in index.keys() if not k.startswith(('accession_ids', 'embedding_shape', 'clean_count'))]
        print(f"  Sample accession keys: {accession_keys[:10]}")
        print(f"Total accession keys: {len(accession_keys)}")
        embedded_accs = accession_keys[:20]

    # Check overlaps
    print(f"\n🔍 OVERLAP ANALYSIS:")
    viz_set = set(viz_accessions[:20])  # Use sample for comparison
    s_set = set(s_accessions[:20])
    embedded_set = set(embedded_accs)

    viz_embedded_overlap = viz_set.intersection(embedded_set)
    s_embedded_overlap = s_set.intersection(embedded_set)

    print(f"Visualization FASTA ∩ S embeddings: {len(viz_embedded_overlap)} / {len(viz_set)}")
    print(f"S segment FASTA ∩ S embeddings: {len(s_embedded_overlap)} / {len(s_set)}")

    if viz_embedded_overlap:
        print(f"  Sample matches: {list(viz_embedded_overlap)[:5]}")
    else:
        print(f"  NO MATCHES found between visualization FASTA and embeddings")

    if s_embedded_overlap:
        print(f"  S segment matches: {list(s_embedded_overlap)[:5]}")
    else:
        print(f"  NO MATCHES found between S segment FASTA and embeddings")

    return 0

if __name__ == "__main__":
    exit(debug_accession_matching())