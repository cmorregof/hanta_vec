#!/usr/bin/env python3
"""
Debug F8 cross-segment isolate matching.
Check if isolate and strain fields are populated for different species.
"""

import pandas as pd
from pathlib import Path

def debug_isolate_matching():
    """Debug the isolate matching logic for F8."""

    print("="*70)
    print("DEBUGGING F8 CROSS-SEGMENT ISOLATE MATCHING")
    print("="*70)

    # Load aligned metadata
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded metadata: {len(metadata)} sequences")

    # Check isolate and strain field population by species
    print(f"\nISOLATE FIELD POPULATION BY SPECIES:")
    print(f"{'-'*60}")

    for species in metadata['species'].unique():
        if species == 'Unknown':
            continue

        species_meta = metadata[metadata['species'] == species]

        # Count non-empty isolate fields
        isolate_populated = species_meta['isolate'].notna() & (species_meta['isolate'] != '')
        strain_populated = species_meta['strain'].notna() & (species_meta['strain'] != '')

        isolate_count = isolate_populated.sum()
        strain_count = strain_populated.sum()
        total = len(species_meta)

        print(f"{species}: {total} sequences")
        print(f"  Isolate populated: {isolate_count}/{total} ({isolate_count/total*100:.1f}%)")
        print(f"  Strain populated: {strain_count}/{total} ({strain_count/total*100:.1f}%)")

        # Show samples
        if isolate_count > 0:
            sample_isolates = species_meta[isolate_populated]['isolate'].head(3).tolist()
            print(f"  Sample isolates: {sample_isolates}")
        if strain_count > 0:
            sample_strains = species_meta[strain_populated]['strain'].head(3).tolist()
            print(f"  Sample strains: {sample_strains}")
        print()

    # Debug the current matching logic
    print(f"TESTING CURRENT MATCHING LOGIC:")
    print(f"{'-'*60}")

    # Group by isolate and species to find cross-segment matches
    isolate_segment_map = {}

    for _, row in metadata.iterrows():
        isolate = row.get('isolate', '')
        species = row['species']
        segment = row['segment_used']
        genbank_id = row['genbank_id']

        # Current logic: Create key from species + isolate
        key = f"{species}_{isolate}".strip('_')

        if key not in isolate_segment_map:
            isolate_segment_map[key] = {}

        if segment not in isolate_segment_map[key]:
            isolate_segment_map[key][segment] = []

        isolate_segment_map[key][segment].append({
            'genbank_id': genbank_id,
            'embedding_index': row.name
        })

    # Find isolates with multiple segments
    matched_isolates = {}

    for key, segments in isolate_segment_map.items():
        if len(segments) >= 2:  # Has at least 2 segments
            species = key.split('_')[0] if '_' in key else 'Unknown'

            if species not in matched_isolates:
                matched_isolates[species] = []

            matched_isolates[species].append({
                'isolate_key': key,
                'segments': list(segments.keys()),
                'segment_counts': {seg: len(seqs) for seg, seqs in segments.items()}
            })

    print(f"CURRENT MATCHING RESULTS:")
    total_matched = 0
    for species, isolates in matched_isolates.items():
        count = len(isolates)
        total_matched += count
        print(f"  {species}: {count} matched isolates")

    print(f"Total matched isolates: {total_matched}")

    # Show sample matched isolates
    print(f"\nSAMPLE OF FIRST 10 MATCHED ISOLATES:")
    print(f"{'-'*60}")

    sample_count = 0
    for species, isolates in matched_isolates.items():
        for isolate_info in isolates[:5]:  # Max 5 per species
            if sample_count >= 10:
                break

            isolate_key = isolate_info['isolate_key']
            segments = isolate_info['segments']
            counts = isolate_info['segment_counts']

            segment_str = '+'.join(segments)
            count_str = ', '.join([f"{seg}:{counts[seg]}" for seg in segments])

            print(f"{sample_count+1:2d}. {species}: {isolate_key}")
            print(f"    Segments: {segment_str} ({count_str})")

            sample_count += 1

        if sample_count >= 10:
            break

    # Alternative matching strategies
    print(f"\nALTERNATIVE MATCHING STRATEGIES:")
    print(f"{'-'*60}")

    # Strategy 1: Match by strain instead of isolate
    print("Strategy 1: Match by strain field")
    strain_matches = 0
    for species in ['Seoul', 'Andes', 'Puumala']:
        species_meta = metadata[metadata['species'] == species]
        strains_with_segments = {}

        for _, row in species_meta.iterrows():
            strain = row.get('strain', '')
            if strain and strain.strip():
                segment = row['segment_used']
                if strain not in strains_with_segments:
                    strains_with_segments[strain] = set()
                strains_with_segments[strain].add(segment)

        multi_segment_strains = {strain: segments for strain, segments in strains_with_segments.items() if len(segments) >= 2}
        strain_matches += len(multi_segment_strains)

        print(f"  {species}: {len(multi_segment_strains)} strains with multiple segments")
        if multi_segment_strains:
            for strain, segments in list(multi_segment_strains.items())[:3]:
                print(f"    {strain}: {'+'.join(segments)}")

    print(f"Total strain-based matches: {strain_matches}")

    # Strategy 2: Match by GenBank accession prefix
    print(f"\nStrategy 2: Match by GenBank accession prefix")
    prefix_matches = 0
    for species in ['Seoul', 'Andes', 'Puumala']:
        species_meta = metadata[metadata['species'] == species]
        prefixes_with_segments = {}

        for _, row in species_meta.iterrows():
            genbank_id = row['genbank_id']
            # Extract prefix (everything before the dot)
            prefix = genbank_id.split('.')[0] if '.' in genbank_id else genbank_id
            segment = row['segment_used']

            if prefix not in prefixes_with_segments:
                prefixes_with_segments[prefix] = set()
            prefixes_with_segments[prefix].add(segment)

        multi_segment_prefixes = {prefix: segments for prefix, segments in prefixes_with_segments.items() if len(segments) >= 2}
        prefix_matches += len(multi_segment_prefixes)

        print(f"  {species}: {len(multi_segment_prefixes)} accession prefixes with multiple segments")
        if multi_segment_prefixes:
            for prefix, segments in list(multi_segment_prefixes.items())[:3]:
                print(f"    {prefix}: {'+'.join(segments)}")

    print(f"Total prefix-based matches: {prefix_matches}")

def main():
    debug_isolate_matching()
    return 0

if __name__ == "__main__":
    exit(main())