#!/usr/bin/env python3
"""
Step 2: Identify what needs recomputation.
Compare clean embedding index against full metadata.
"""

import json
import pandas as pd
from pathlib import Path

def load_clean_index(embeddings_dir, segment):
    """Load clean embedding index for a segment."""

    index_file = embeddings_dir / f"embeddings_{segment}_index.json"
    if not index_file.exists():
        print(f"❌ Index file not found: {index_file}")
        return []

    with open(index_file, 'r') as f:
        index_data = json.load(f)

    accession_ids = index_data.get('accession_ids', [])
    print(f"  {segment}: Loaded {len(accession_ids)} valid accession IDs")
    return accession_ids

def identify_missing_sequences(metadata, valid_embeddings_per_segment):
    """Identify which sequences need recomputation per segment."""

    print(f"\n🔍 IDENTIFYING MISSING SEQUENCES")
    print("=" * 50)

    missing_per_segment = {}

    for segment in ['S', 'M', 'L']:
        print(f"\n📊 {segment} Segment Analysis:")

        # Get all sequences for this segment from metadata
        segment_metadata = metadata[metadata['segment_used'] == segment]
        all_accessions = set(segment_metadata['genbank_id'].tolist())

        # Get valid embeddings for this segment
        valid_accessions = set(valid_embeddings_per_segment.get(segment, []))

        # Find missing
        missing_accessions = all_accessions - valid_accessions

        # Report
        print(f"  Total sequences in metadata: {len(all_accessions)}")
        print(f"  Valid embeddings found: {len(valid_accessions)}")
        print(f"  Missing embeddings: {len(missing_accessions)}")
        print(f"  Recomputation needed: {len(missing_accessions)}/{len(all_accessions)} ({len(missing_accessions)/len(all_accessions)*100:.1f}%)")

        # Sample missing sequences
        if len(missing_accessions) > 0:
            sample_missing = sorted(list(missing_accessions))[:5]
            print(f"  Sample missing: {sample_missing}")

        # Sample valid sequences
        if len(valid_accessions) > 0:
            sample_valid = sorted(list(valid_accessions))[:3]
            print(f"  Sample valid: {sample_valid}")

        missing_per_segment[segment] = {
            'missing_accessions': sorted(list(missing_accessions)),
            'valid_accessions': sorted(list(valid_accessions)),
            'missing_count': len(missing_accessions),
            'valid_count': len(valid_accessions),
            'total_count': len(all_accessions)
        }

    return missing_per_segment

def save_recomputation_plan(missing_per_segment, embeddings_dir):
    """Save the recomputation plan to file."""

    print(f"\n💾 SAVING RECOMPUTATION PLAN")
    print("=" * 50)

    plan_file = embeddings_dir / "recomputation_plan.json"

    # Create comprehensive plan
    plan = {
        'summary': {
            'total_sequences_needed': sum(data['missing_count'] for data in missing_per_segment.values()),
            'total_sequences_valid': sum(data['valid_count'] for data in missing_per_segment.values()),
            'segments': {}
        },
        'missing_by_segment': missing_per_segment
    }

    for segment, data in missing_per_segment.items():
        plan['summary']['segments'][segment] = {
            'missing': data['missing_count'],
            'valid': data['valid_count'],
            'total': data['total_count'],
            'percent_missing': round(data['missing_count'] / data['total_count'] * 100, 1)
        }

    with open(plan_file, 'w') as f:
        json.dump(plan, f, indent=2)

    print(f"✓ Saved recomputation plan: {plan_file}")
    return plan

def validate_no_overlaps(missing_per_segment):
    """Ensure no sequence appears in multiple segments (sanity check)."""

    print(f"\n✅ VALIDATION: CHECKING FOR OVERLAPS")
    print("=" * 50)

    all_missing = []
    all_valid = []

    for segment, data in missing_per_segment.items():
        all_missing.extend(data['missing_accessions'])
        all_valid.extend(data['valid_accessions'])

    # Check for duplicates within missing
    missing_set = set(all_missing)
    if len(all_missing) != len(missing_set):
        print(f"⚠️  Warning: {len(all_missing) - len(missing_set)} duplicate accessions in missing list")
    else:
        print(f"✓ No duplicate missing accessions")

    # Check for duplicates within valid
    valid_set = set(all_valid)
    if len(all_valid) != len(valid_set):
        print(f"⚠️  Warning: {len(all_valid) - len(valid_set)} duplicate accessions in valid list")
    else:
        print(f"✓ No duplicate valid accessions")

    # Check for overlap between missing and valid (should be zero)
    overlap = missing_set & valid_set
    if len(overlap) > 0:
        print(f"❌ Critical error: {len(overlap)} accessions appear in both missing and valid lists!")
        print(f"   Sample overlaps: {sorted(list(overlap))[:5]}")
    else:
        print(f"✓ No overlap between missing and valid (good)")

def main():
    print("=" * 70)
    print("STEP 2: IDENTIFY SEQUENCES NEEDING RECOMPUTATION")
    print("=" * 70)

    # Setup paths
    project_root = Path(__file__).parent.parent
    embeddings_dir = project_root / "results" / "embeddings"
    processed_dir = project_root / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    # Load metadata
    if not metadata_path.exists():
        print(f"❌ Metadata not found: {metadata_path}")
        return 1

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded metadata: {len(metadata)} sequences")

    # Load clean embedding indices
    print(f"\n📂 LOADING CLEAN EMBEDDING INDICES")
    print("=" * 50)

    valid_embeddings_per_segment = {}
    for segment in ['S', 'M', 'L']:
        accession_ids = load_clean_index(embeddings_dir, segment)
        valid_embeddings_per_segment[segment] = accession_ids

    total_valid = sum(len(ids) for ids in valid_embeddings_per_segment.values())
    print(f"\nTotal valid embeddings loaded: {total_valid}")

    # Identify missing sequences
    missing_per_segment = identify_missing_sequences(metadata, valid_embeddings_per_segment)

    # Create and save recomputation plan
    plan = save_recomputation_plan(missing_per_segment, embeddings_dir)

    # Validate for overlaps
    validate_no_overlaps(missing_per_segment)

    # Final summary
    print(f"\n📋 RECOMPUTATION SUMMARY")
    print("=" * 50)
    total_missing = plan['summary']['total_sequences_needed']
    total_sequences = len(metadata)

    print(f"Sequences with valid embeddings: {total_valid}/{total_sequences} ({total_valid/total_sequences*100:.1f}%)")
    print(f"Sequences needing recomputation: {total_missing}/{total_sequences} ({total_missing/total_sequences*100:.1f}%)")

    print(f"\nPer segment breakdown:")
    for segment in ['S', 'M', 'L']:
        seg_data = plan['summary']['segments'][segment]
        print(f"  {segment}: {seg_data['missing']}/{seg_data['total']} missing ({seg_data['percent_missing']}%)")

    print(f"\n✅ Step 2 complete! Recomputation plan created.")
    print(f"Expected computation time: 15-30 minutes for ~{total_missing} sequences")
    print(f"Next: Run Step 3 to fix the embedding script before recomputing")

    return 0

if __name__ == "__main__":
    exit(main())