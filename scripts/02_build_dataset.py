#!/usr/bin/env python3
"""
Phase 1: Build three-segment curated dataset from NCBI sequences.

Three-segment strategy:
- S-segment: Nucleocapsid protein (N, ~430 aa)
- M-segment: Glycoprotein precursor (GPC, ~900-1150 aa)
- L-segment: RNA polymerase (RdRp, ~2100-2200 aa)

Steps:
1. Fetch sequences from all three segments for all taxa
2. Apply segment-specific QC filters
3. Remove duplicates within each segment (not across segments)
4. Create Level 1 dataset with cross-segment matching
5. Generate comprehensive reports

Outputs:
- data/processed/{S,M,L}_sequences_level1.fasta
- data/processed/metadata_level1.tsv (combined all segments)
- results/manifests/dataset_manifest_level1.tsv
- results/manifests/qc_report.json
"""

import json
import sys
import logging
from pathlib import Path
from typing import Dict, Tuple

import yaml
import pandas as pd
from Bio import Entrez, SeqIO

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import (
    load_config,
    fetch_three_segment_sequences,
)


def setup_logging(log_dir: Path) -> logging.Logger:
    """Setup logging."""
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "phase1.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def check_sequence_quality(sequence: str, segment: str) -> bool:
    """Apply segment-specific QC filters."""
    # Check for ambiguous amino acids (15% threshold)
    ambiguous_count = sum(sequence.count(aa) for aa in ['X', 'B', 'Z', 'J'])
    if (ambiguous_count / len(sequence)) > 0.15:
        return False

    # Check for stop codons
    if '*' in sequence.replace('*', '', len(sequence) - 1):  # Allow terminal stop
        return False

    # Segment-specific length checks
    length_ranges = {
        'S': (380, 500),    # Nucleocapsid
        'M': (400, 1300),   # Glycoprotein
        'L': (1800, 2300),  # Polymerase
    }

    min_len, max_len = length_ranges[segment]
    return min_len <= len(sequence) <= max_len


def remove_near_duplicates_within_segment(segment_data: Dict, threshold: float = 0.99) -> Tuple[Dict, Dict]:
    """Remove near-duplicates within a segment using sequence identity."""
    import difflib

    sequences = [(record_id, data['sequence']) for record_id, data in segment_data.items()]
    to_remove = set()
    removed_data = {}

    for i, (id1, seq1) in enumerate(sequences):
        if id1 in to_remove:
            continue

        for j, (id2, seq2) in enumerate(sequences[i+1:], i+1):
            if id2 in to_remove:
                continue

            # Calculate sequence similarity
            similarity = difflib.SequenceMatcher(None, seq1, seq2).ratio()

            if similarity >= threshold:
                # Keep the first one, remove the second
                to_remove.add(id2)
                removed_data[id2] = segment_data[id2]

    # Return filtered data
    kept_data = {k: v for k, v in segment_data.items() if k not in to_remove}
    return kept_data, removed_data


def find_cross_segment_matches(segment_data_dict: Dict) -> Dict:
    """Find isolates that appear in multiple segments."""
    # Build mapping of isolate/strain -> segments -> record_ids
    isolate_map = {}

    for segment, data in segment_data_dict.items():
        for record_id, record_data in data.items():
            isolate = record_data.get('isolate')
            strain = record_data.get('strain')
            species = record_data['species']

            # Use isolate or strain as key
            key = isolate or strain
            if key:
                key = f"{species}_{key}"
                if key not in isolate_map:
                    isolate_map[key] = {}
                if segment not in isolate_map[key]:
                    isolate_map[key][segment] = []
                isolate_map[key][segment].append(record_id)

    return isolate_map


def main():
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    log_dir = Path(__file__).parent.parent / "logs"
    logger = setup_logging(log_dir)

    logger.info("="*80)
    logger.info("PHASE 1: BUILD THREE-SEGMENT DATASET")
    logger.info("="*80)

    config = load_config(config_path)

    # Step 1: Fetch sequences from all three segments
    logger.info("\n1. FETCH SEQUENCES FROM ALL THREE SEGMENTS")
    output_dir = Path(__file__).parent.parent / "data" / "raw" / "proteins"
    raw_data = fetch_three_segment_sequences(config_path, output_dir, use_broad_search=True)

    total_fetched = sum(len(segment_data) for segment_data in raw_data.values())
    logger.info(f"\nTotal fetched across all segments: {total_fetched}")

    # Step 2: Apply segment-specific QC
    logger.info("\n2. APPLY SEGMENT-SPECIFIC QC")
    qc_passed = {}
    qc_failed = {}

    for segment, segment_data in raw_data.items():
        passed = {}
        failed = {}

        for record_id, record_data in segment_data.items():
            sequence = record_data['sequence']

            if check_sequence_quality(sequence, segment):
                passed[record_id] = record_data
            else:
                failed[record_id] = record_data

        qc_passed[segment] = passed
        qc_failed[segment] = failed

        logger.info(f"  {segment}-segment: {len(passed)} passed, {len(failed)} failed QC")

    # Step 3: Remove near-duplicates within each segment
    logger.info("\n3. REMOVE NEAR-DUPLICATES (≥99% identity within segment)")
    final_data = {}
    removed_dups = {}

    for segment, segment_data in qc_passed.items():
        kept, removed = remove_near_duplicates_within_segment(segment_data, threshold=0.99)
        final_data[segment] = kept
        removed_dups[segment] = removed

        logger.info(f"  {segment}-segment: {len(kept)} kept, {len(removed)} duplicates removed")

    # Step 4: Cross-segment analysis
    logger.info("\n4. CROSS-SEGMENT MATCHING")
    cross_matches = find_cross_segment_matches(final_data)

    # Count isolates with all three segments
    all_three_segments = 0
    for isolate, segments in cross_matches.items():
        if len(segments) == 3:
            all_three_segments += 1

    logger.info(f"  Matched isolates: {len(cross_matches)}")
    logger.info(f"  Isolates with all three segments: {all_three_segments}")

    # Step 5: Build combined metadata and save files
    logger.info("\n5. SAVE PROCESSED DATA")
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    # Save FASTA files per segment
    all_metadata = []
    for segment, segment_data in final_data.items():
        fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"

        with open(fasta_path, 'w') as f:
            for record_id, record_data in segment_data.items():
                # Write FASTA record
                header = f">{record_id}_{segment} {record_data['species']} | {record_data['clade']} | {record_data['method']}"
                f.write(f"{header}\n{record_data['sequence']}\n")

                # Collect metadata
                metadata_row = {
                    'record_id': f"{record_id}_{segment}",
                    'genbank_id': record_id,
                    'segment': segment,
                    'species': record_data['species'],
                    'taxon_id': record_data['taxon_id'],
                    'clade': record_data['clade'],
                    'isolate': record_data.get('isolate', ''),
                    'strain': record_data.get('strain', ''),
                    'length': record_data['length'],
                    'method': record_data['method'],
                    'description': record_data['description'],
                }
                all_metadata.append(metadata_row)

        logger.info(f"  Saved {fasta_path.name}: {len(segment_data)} sequences")

    # Save combined metadata
    import pandas as pd
    metadata_df = pd.DataFrame(all_metadata)
    metadata_path = processed_dir / "metadata_level1.tsv"
    metadata_df.to_csv(metadata_path, sep='\t', index=False)
    logger.info(f"  Saved {metadata_path.name}: {len(metadata_df)} records")

    # Step 6: Generate comprehensive QC report
    logger.info("\n6. GENERATE QC REPORT")

    # Species × Segment cross-table
    species_segment_table = metadata_df.pivot_table(
        index='species', columns='segment', values='record_id', aggfunc='count', fill_value=0
    )

    # Method breakdown per segment
    method_breakdown = {}
    for segment in ['S', 'M', 'L']:
        segment_df = metadata_df[metadata_df['segment'] == segment]
        method_breakdown[segment] = segment_df['method'].value_counts().to_dict()

    # Cross-segment matches
    cross_segment_stats = {
        'total_isolates': len(cross_matches),
        'with_all_three': all_three_segments,
        'with_two_segments': sum(1 for segments in cross_matches.values() if len(segments) == 2),
        'single_segment_only': sum(1 for segments in cross_matches.values() if len(segments) == 1),
    }

    qc_report = {
        'total_fetched': total_fetched,
        'total_after_qc': sum(len(segment_data) for segment_data in final_data.values()),
        'per_segment_counts': {segment: len(data) for segment, data in final_data.items()},
        'species_segment_table': species_segment_table.to_dict(),
        'method_breakdown': method_breakdown,
        'cross_segment_matches': cross_segment_stats,
        'qc_passed_per_segment': {segment: len(data) for segment, data in qc_passed.items()},
        'qc_failed_per_segment': {segment: len(data) for segment, data in qc_failed.items()},
        'duplicates_removed_per_segment': {segment: len(data) for segment, data in removed_dups.items()},
    }

    # Save QC report
    manifests_dir = Path(__file__).parent.parent / "results" / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    with open(manifests_dir / "qc_report.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    logger.info("\n" + "="*50)
    logger.info("QC REPORT SUMMARY")
    logger.info("="*50)
    logger.info(f"Total sequences after QC: {qc_report['total_after_qc']}")
    for segment, count in qc_report['per_segment_counts'].items():
        logger.info(f"  {segment}-segment: {count} sequences")
    logger.info(f"\nCross-segment matches:")
    logger.info(f"  All three segments: {cross_segment_stats['with_all_three']} isolates")
    logger.info(f"  Two segments: {cross_segment_stats['with_two_segments']} isolates")

    logger.info("\n" + "="*80)
    logger.info("✓ PHASE 1 COMPLETE")
    logger.info("="*80)

    # Success criteria for three segments
    total_sequences = qc_report['total_after_qc']
    if total_sequences < 1500:
        logger.warning(f"⚠ Total sequences ({total_sequences}) < 1500 target")
        return 1

    min_per_segment = min(qc_report['per_segment_counts'].values())
    if min_per_segment < 500:
        logger.warning(f"⚠ Minimum per segment ({min_per_segment}) < 500 target")
        return 1

    if cross_segment_stats['with_all_three'] < 200:
        logger.warning(f"⚠ Cross-segment matches ({cross_segment_stats['with_all_three']}) < 200 target")
        return 1

    logger.info("✓ Three-segment dataset ready for Phase 2 (embeddings)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
