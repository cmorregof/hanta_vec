#!/usr/bin/env python3
"""
Comprehensive diagnostic of the embedding pipeline to identify
why 88.8% of embeddings are zero vectors.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

def analyze_embedding_file(file_path, segment_name, metadata=None):
    """Analyze a single embedding file for zero vectors."""

    print(f"\n🔍 ANALYZING {segment_name} EMBEDDINGS: {file_path}")
    print("=" * 60)

    if not file_path.exists():
        print(f"❌ File does not exist: {file_path}")
        return

    # Load embeddings
    embeddings = np.load(file_path)
    print(f"Shape: {embeddings.shape}")

    # Calculate norms
    norms = np.linalg.norm(embeddings, axis=1)
    zero_threshold = 1e-6

    zero_rows = np.sum(norms < zero_threshold)
    non_zero_rows = len(embeddings) - zero_rows

    print(f"Zero rows (norm < {zero_threshold}): {zero_rows}/{len(embeddings)} ({zero_rows/len(embeddings)*100:.1f}%)")
    print(f"Non-zero rows: {non_zero_rows}/{len(embeddings)} ({non_zero_rows/len(embeddings)*100:.1f}%)")

    if non_zero_rows > 0:
        non_zero_norms = norms[norms >= zero_threshold]
        print(f"Mean norm of non-zero rows: {non_zero_norms.mean():.6f}")
        print(f"Min/Max non-zero norms: {non_zero_norms.min():.6f} / {non_zero_norms.max():.6f}")

        # Sample some non-zero embeddings
        non_zero_indices = np.where(norms >= zero_threshold)[0]
        sample_indices = non_zero_indices[:5] if len(non_zero_indices) >= 5 else non_zero_indices

        print(f"Sample non-zero embeddings (first {len(sample_indices)}):")
        for i, idx in enumerate(sample_indices):
            norm_val = norms[idx]
            print(f"  Index {idx}: norm = {norm_val:.6f}")

            # Try to match with metadata if provided
            if metadata is not None:
                # Filter metadata for this segment
                seg_metadata = metadata[metadata['segment_used'] == segment_name]
                if len(seg_metadata) > idx:
                    row = seg_metadata.iloc[idx]
                    genbank_id = row.get('genbank_id', 'Unknown')
                    print(f"    → {genbank_id}")

    # Check for patterns in zero embeddings
    zero_indices = np.where(norms < zero_threshold)[0]
    if len(zero_indices) > 0:
        print(f"Zero embedding indices (first 10): {zero_indices[:10].tolist()}")
        print(f"Zero embedding indices (last 10): {zero_indices[-10:].tolist()}")

        # Check if zeros are clustered
        if len(zero_indices) > 1:
            consecutive_gaps = np.diff(zero_indices)
            has_consecutive = np.any(consecutive_gaps == 1)
            print(f"Has consecutive zero embeddings: {has_consecutive}")

            # Check if zeros are at the end (padding pattern)
            last_non_zero = np.where(norms >= zero_threshold)[0]
            if len(last_non_zero) > 0:
                last_valid_idx = last_non_zero[-1]
                trailing_zeros = len(embeddings) - 1 - last_valid_idx
                print(f"Trailing zeros after last valid embedding: {trailing_zeros}")

def check_embedding_reports(results_dir):
    """Check embedding generation reports for failures."""

    print(f"\n📊 CHECKING EMBEDDING REPORTS")
    print("=" * 60)

    # Check main embedding report
    report_path = results_dir / "embeddings" / "embedding_report.json"
    if report_path.exists():
        with open(report_path, 'r') as f:
            report = json.load(f)

        print(f"Found embedding report: {report_path}")
        print(f"Report keys: {list(report.keys())}")

        # Check for failure indicators
        for key, value in report.items():
            if isinstance(value, dict):
                if 'failed' in value or 'errors' in value or 'failures' in value:
                    print(f"⚠️  {key}: {value}")
                elif 'successful' in value or 'processed' in value:
                    print(f"✓ {key}: {value}")

        # Look for segment-specific reports
        for segment in ['S', 'M', 'L']:
            if segment in report:
                seg_data = report[segment]
                print(f"{segment} segment report: {seg_data}")
    else:
        print(f"❌ No embedding report found at {report_path}")

    # Check for segment-specific reports
    for segment in ['S', 'M', 'L']:
        seg_report_path = results_dir / "embeddings" / f"embedding_report_{segment}.json"
        if seg_report_path.exists():
            with open(seg_report_path, 'r') as f:
                seg_report = json.load(f)
            print(f"Found {segment} segment report: {seg_report}")

def check_embedding_cache(results_dir):
    """Check embedding cache for issues."""

    print(f"\n💾 CHECKING EMBEDDING CACHE")
    print("=" * 60)

    cache_dir = results_dir / "embeddings" / "cache"
    if not cache_dir.exists():
        print(f"❌ No cache directory found at {cache_dir}")
        return

    print(f"Cache directory exists: {cache_dir}")

    # List cache files
    cache_files = list(cache_dir.glob("*.npy"))
    print(f"Cache files found: {len(cache_files)}")

    if len(cache_files) > 0:
        print(f"Sample cache files:")
        for i, cache_file in enumerate(cache_files[:5]):
            print(f"  {cache_file.name}")

            # Check if cached embedding is zero
            try:
                cached_emb = np.load(cache_file)
                norm = np.linalg.norm(cached_emb)
                print(f"    Shape: {cached_emb.shape}, Norm: {norm:.6f}")

                if norm < 1e-6:
                    print(f"    ⚠️  Zero cached embedding!")

            except Exception as e:
                print(f"    ❌ Error loading cache: {e}")

    # Check cache statistics
    zero_cache_count = 0
    valid_cache_count = 0

    for cache_file in cache_files:
        try:
            cached_emb = np.load(cache_file)
            norm = np.linalg.norm(cached_emb)

            if norm < 1e-6:
                zero_cache_count += 1
            else:
                valid_cache_count += 1

        except:
            continue

    total_cache = zero_cache_count + valid_cache_count
    if total_cache > 0:
        print(f"Cache statistics:")
        print(f"  Valid cached embeddings: {valid_cache_count}/{total_cache} ({valid_cache_count/total_cache*100:.1f}%)")
        print(f"  Zero cached embeddings: {zero_cache_count}/{total_cache} ({zero_cache_count/total_cache*100:.1f}%)")

def check_metadata_alignment(metadata_path, embeddings_dir):
    """Check if metadata aligns with embeddings."""

    print(f"\n🔗 CHECKING METADATA-EMBEDDING ALIGNMENT")
    print("=" * 60)

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Total metadata entries: {len(metadata)}")

    # Check segment distribution in metadata
    segment_counts = metadata['segment_used'].value_counts()
    print(f"Metadata segment distribution:")
    for segment, count in segment_counts.items():
        print(f"  {segment}: {count}")

    # Check if embedding files match metadata counts
    for segment in ['S', 'M', 'L']:
        emb_file = embeddings_dir / f"embeddings_{segment}.npy"
        if emb_file.exists():
            embeddings = np.load(emb_file)
            metadata_count = segment_counts.get(segment, 0)
            embedding_count = len(embeddings)

            print(f"{segment}: metadata={metadata_count}, embeddings={embedding_count}", end="")
            if metadata_count == embedding_count:
                print(" ✓")
            else:
                print(f" ❌ MISMATCH!")
                print(f"  Difference: {embedding_count - metadata_count}")

def identify_likely_cause():
    """Analyze results to identify the most likely cause."""

    print(f"\n🎯 LIKELY CAUSE ANALYSIS")
    print("=" * 60)
    print("Based on the diagnostic results, the most likely causes are:")
    print("1. ESM-2 forward pass returning zeros for long sequences")
    print("2. Append operation creating index misalignment with zero padding")
    print("3. Consolidation script using np.zeros padding for unmatched sequences")
    print("4. Cache corruption with zero embeddings being saved")
    print("\nReview the diagnostic output above to determine which scenario matches the evidence.")

def main():
    print("=" * 70)
    print("COMPREHENSIVE EMBEDDING PIPELINE DIAGNOSTIC")
    print("=" * 70)
    print("Diagnosing why 88.8% of embeddings are zero vectors...")

    # Setup paths
    project_root = Path(__file__).parent.parent
    results_dir = project_root / "results"
    embeddings_dir = results_dir / "embeddings"
    processed_dir = project_root / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    # Load metadata for reference
    metadata = None
    if metadata_path.exists():
        metadata = pd.read_csv(metadata_path, sep="\t")

    # 1. Analyze individual embedding files
    embedding_files = {
        'S': embeddings_dir / "embeddings_S.npy",
        'M': embeddings_dir / "embeddings_M.npy",
        'L': embeddings_dir / "embeddings_L.npy"
    }

    for segment, file_path in embedding_files.items():
        analyze_embedding_file(file_path, segment, metadata)

    # Also check the consolidated file
    consolidated_file = embeddings_dir / "embeddings_level1.npy"
    if consolidated_file.exists():
        analyze_embedding_file(consolidated_file, "CONSOLIDATED", metadata)

    # 2. Check embedding reports
    check_embedding_reports(results_dir)

    # 3. Check embedding cache
    check_embedding_cache(results_dir)

    # 4. Check metadata alignment
    if metadata is not None:
        check_metadata_alignment(metadata_path, embeddings_dir)

    # 5. Provide analysis of likely causes
    identify_likely_cause()

    print(f"\n" + "=" * 70)
    print("DIAGNOSTIC COMPLETE")
    print("=" * 70)
    print("🚨 DO NOT RECOMPUTE EMBEDDINGS until root cause is identified and fixed!")

    return 0

if __name__ == "__main__":
    exit(main())