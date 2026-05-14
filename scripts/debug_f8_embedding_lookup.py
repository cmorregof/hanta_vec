#!/usr/bin/env python3
"""
Debug F8 embedding lookup for Seoul, Andes, Puumala isolates.
Check if embedding indices match across segments.
"""

import numpy as np
import pandas as pd
from pathlib import Path

def find_matched_isolates_debug(metadata):
    """Find isolates with sequences in multiple segments - DEBUG VERSION."""

    print("  Finding matched isolates for debugging...")

    # Group by isolate and species to find cross-segment matches
    isolate_segment_map = {}

    for _, row in metadata.iterrows():
        isolate = row.get('isolate', '')
        species = row['species']
        segment = row['segment_used']
        genbank_id = row['genbank_id']

        # Handle NaN isolates by using 'empty' as placeholder
        if pd.isna(isolate) or isolate == '':
            isolate = 'empty'

        # Create a key for matching
        key = f"{species}_{isolate}"

        if key not in isolate_segment_map:
            isolate_segment_map[key] = {}

        if segment not in isolate_segment_map[key]:
            isolate_segment_map[key][segment] = []

        isolate_segment_map[key][segment].append({
            'genbank_id': genbank_id,
            'embedding_index': row.name,  # pandas index
            'record_id': row.get('record_id', ''),
            'row_data': row
        })

    # Find isolates with multiple segments
    matched_isolates = {}

    for key, segments in isolate_segment_map.items():
        if len(segments) >= 2:  # Has at least 2 segments
            # Extract species from key
            species = key.split('_')[0]

            if species == 'Unknown':
                continue  # Skip Unknown species for cleaner analysis

            if species not in matched_isolates:
                matched_isolates[species] = []

            matched_isolates[species].append({
                'isolate_key': key,
                'segments': segments
            })

    return matched_isolates

def debug_embedding_lookup(metadata, embeddings_dir):
    """Debug embedding lookup for major species."""

    print("  Debugging embedding lookup for Seoul, Andes, Puumala...")

    # Check if we have separate segment embeddings
    embeddings_files = {
        'S': embeddings_dir / "embeddings_S.npy",
        'M': embeddings_dir / "embeddings_M.npy",
        'L': embeddings_dir / "embeddings_L.npy"
    }

    # Also check for single embedding file
    single_embedding_file = embeddings_dir / "embeddings_level1.npy"

    print(f"  Embedding files status:")
    for seg, path in embeddings_files.items():
        exists = "✓" if path.exists() else "❌"
        print(f"    {seg}: {path} {exists}")

    single_exists = "✓" if single_embedding_file.exists() else "❌"
    print(f"    Single: {single_embedding_file} {single_exists}")

    # Load embeddings
    embeddings = {}
    if single_embedding_file.exists():
        print(f"  Loading single embedding file...")
        single_embeddings = np.load(single_embedding_file)
        print(f"    Shape: {single_embeddings.shape}")
        # Use single embeddings for all segments
        embeddings['single'] = single_embeddings

    for seg, path in embeddings_files.items():
        if path.exists():
            embeddings[seg] = np.load(path)
            print(f"    Loaded {seg}: shape {embeddings[seg].shape}")

    # Find matched isolates
    matched_isolates = find_matched_isolates_debug(metadata)

    # Debug Seoul isolates specifically
    target_species = ['Seoul', 'Andes', 'Puumala']

    for species in target_species:
        if species not in matched_isolates:
            print(f"\n  ❌ No matched isolates found for {species}")
            continue

        species_isolates = matched_isolates[species]
        print(f"\n  🔍 DEBUGGING {species} ({len(species_isolates)} matched isolates)")
        print(f"  {'='*60}")

        # Show first 5 isolates
        for i, isolate_info in enumerate(species_isolates[:5]):
            print(f"\n  Isolate {i+1}: {isolate_info['isolate_key']}")
            segments = isolate_info['segments']

            for seg_name, seq_list in segments.items():
                seq_info = seq_list[0]  # Take first sequence
                genbank_id = seq_info['genbank_id']
                embedding_idx = seq_info['embedding_index']
                record_id = seq_info.get('record_id', '')

                print(f"    {seg_name}: {genbank_id} (idx={embedding_idx}, record={record_id})")

                # Check if embedding exists
                if 'single' in embeddings:
                    # Using single embedding file - check bounds
                    if embedding_idx < len(embeddings['single']):
                        emb = embeddings['single'][embedding_idx]
                        print(f"      ✓ Single embedding found: shape {emb.shape}, norm {np.linalg.norm(emb):.3f}")
                    else:
                        print(f"      ❌ Index {embedding_idx} out of bounds for single embeddings (max {len(embeddings['single'])-1})")

                if seg_name in embeddings:
                    # Check segment-specific embedding
                    seg_embeddings = embeddings[seg_name]

                    # Try to find matching embedding by index
                    if embedding_idx < len(seg_embeddings):
                        emb = seg_embeddings[embedding_idx]
                        print(f"      ✓ {seg_name} embedding found: shape {emb.shape}, norm {np.linalg.norm(emb):.3f}")
                    else:
                        print(f"      ❌ Index {embedding_idx} out of bounds for {seg_name} embeddings (max {len(seg_embeddings)-1})")

            # Try to compute similarities for this isolate
            print(f"    Similarity computation test:")
            segment_pairs = [('S', 'M'), ('S', 'L'), ('M', 'L')]

            for seg1, seg2 in segment_pairs:
                if seg1 in segments and seg2 in segments:
                    try:
                        # Get embedding indices
                        emb1_idx = segments[seg1][0]['embedding_index']
                        emb2_idx = segments[seg2][0]['embedding_index']

                        # Use single embeddings if available, otherwise segment-specific
                        if 'single' in embeddings:
                            if emb1_idx < len(embeddings['single']) and emb2_idx < len(embeddings['single']):
                                emb1 = embeddings['single'][emb1_idx]
                                emb2 = embeddings['single'][emb2_idx]

                                # Compute cosine similarity
                                cos_sim = np.dot(emb1, emb2) / (np.linalg.norm(emb1) * np.linalg.norm(emb2))

                                if np.isfinite(cos_sim):
                                    print(f"      {seg1}-{seg2}: similarity = {cos_sim:.6f} ✓")
                                else:
                                    print(f"      {seg1}-{seg2}: similarity = {cos_sim} ❌ (not finite)")
                            else:
                                print(f"      {seg1}-{seg2}: index out of bounds ❌")
                        else:
                            print(f"      {seg1}-{seg2}: no single embeddings available ❌")

                    except Exception as e:
                        print(f"      {seg1}-{seg2}: ERROR - {e}")

def main():
    print("="*70)
    print("DEBUGGING F8 EMBEDDING LOOKUP")
    print("="*70)

    # Load data
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load aligned metadata
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"
    metadata = pd.read_csv(metadata_path, sep="\t")

    print(f"Loaded metadata: {len(metadata)} sequences")

    # Debug embedding lookup
    debug_embedding_lookup(metadata, embeddings_dir)

    return 0

if __name__ == "__main__":
    exit(main())