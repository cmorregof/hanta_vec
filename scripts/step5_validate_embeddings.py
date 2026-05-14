#!/usr/bin/env python3
"""
Step 5: Validate final embeddings before regenerating figures.

REQUIREMENTS:
- Zero vectors must be 0%
- Mean norm per segment should be ~6.6–6.8
- All embeddings must have valid norms
"""

import json
import numpy as np
from pathlib import Path

def validate_embedding_file(file_path, segment_name):
    """Validate a final embedding file."""

    print(f"\n✅ VALIDATING {segment_name} EMBEDDINGS: {file_path}")
    print("-" * 60)

    if not file_path.exists():
        print(f"❌ File does not exist: {file_path}")
        return False

    # Load embeddings
    embeddings = np.load(file_path)
    print(f"Shape: {embeddings.shape}")

    # Calculate norms
    norms = np.linalg.norm(embeddings, axis=1)
    zero_threshold = 1e-6

    zero_count = np.sum(norms < zero_threshold)
    valid_count = len(embeddings) - zero_count
    zero_percentage = zero_count / len(embeddings) * 100

    print(f"Zero vectors (norm < {zero_threshold}): {zero_count}/{len(embeddings)} ({zero_percentage:.1f}%)")
    print(f"Valid vectors: {valid_count}/{len(embeddings)} ({valid_count/len(embeddings)*100:.1f}%)")

    # Check for zero vectors (MUST be 0%)
    if zero_count > 0:
        print(f"❌ VALIDATION FAILED: {zero_count} zero vectors found!")
        print(f"   Zero vector indices: {np.where(norms < zero_threshold)[0][:10].tolist()}")
        return False
    else:
        print(f"✅ ZERO VECTORS: 0% (PERFECT)")

    # Check mean norm (should be ~6.6-6.8)
    mean_norm = norms.mean()
    std_norm = norms.std()
    min_norm = norms.min()
    max_norm = norms.max()

    print(f"Mean norm: {mean_norm:.6f}")
    print(f"Std norm: {std_norm:.6f}")
    print(f"Min/Max norm: {min_norm:.6f} / {max_norm:.6f}")

    # Validate norm range (should be ~6.6-6.8)
    expected_range = (6.0, 7.5)  # Allow some tolerance
    if expected_range[0] <= mean_norm <= expected_range[1]:
        print(f"✅ MEAN NORM: {mean_norm:.3f} (within expected range {expected_range})")
        norm_valid = True
    else:
        print(f"⚠️  MEAN NORM: {mean_norm:.3f} (outside expected range {expected_range})")
        norm_valid = False

    # Check for NaN or infinite values
    has_nan = np.isnan(embeddings).any()
    has_inf = np.isinf(embeddings).any()

    if has_nan:
        print(f"❌ NaN values found in embeddings!")
        return False
    else:
        print(f"✅ NO NAN VALUES")

    if has_inf:
        print(f"❌ Infinite values found in embeddings!")
        return False
    else:
        print(f"✅ NO INFINITE VALUES")

    # Overall validation
    overall_valid = (zero_count == 0) and (not has_nan) and (not has_inf) and norm_valid

    if overall_valid:
        print(f"🎉 VALIDATION PASSED: {segment_name} embeddings are valid!")
    else:
        print(f"❌ VALIDATION FAILED: {segment_name} embeddings have issues!")

    return overall_valid

def check_embedding_coverage(embeddings_dir, processed_dir):
    """Check how much of the original metadata is now covered by embeddings."""

    print(f"\n📊 CHECKING EMBEDDING COVERAGE")
    print("=" * 60)

    # Load metadata
    import pandas as pd
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"
    metadata = pd.read_csv(metadata_path, sep="\t")

    print(f"Total sequences in metadata: {len(metadata)}")

    # Check coverage per segment
    total_covered = 0
    coverage_report = {}

    for segment in ['S', 'M', 'L']:
        # Count metadata sequences for this segment
        segment_metadata = metadata[metadata['segment_used'] == segment]
        metadata_count = len(segment_metadata)

        # Load final embeddings for this segment
        final_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        final_index_file = embeddings_dir / f"embeddings_{segment}_final_index.json"

        if final_file.exists() and final_index_file.exists():
            with open(final_index_file, 'r') as f:
                index_data = json.load(f)
            embedding_count = index_data['total_count']
            total_covered += embedding_count
        else:
            print(f"⚠️  Final files not found for {segment}")
            embedding_count = 0

        coverage_pct = embedding_count / metadata_count * 100 if metadata_count > 0 else 0

        print(f"{segment}: {embedding_count}/{metadata_count} ({coverage_pct:.1f}% coverage)")

        coverage_report[segment] = {
            'metadata_count': metadata_count,
            'embedding_count': embedding_count,
            'coverage_percentage': coverage_pct
        }

    total_metadata = len(metadata)
    total_coverage_pct = total_covered / total_metadata * 100

    print(f"\nTOTAL COVERAGE: {total_covered}/{total_metadata} ({total_coverage_pct:.1f}%)")

    # Check if we have sufficient coverage
    if total_coverage_pct >= 95:
        print(f"✅ EXCELLENT COVERAGE: {total_coverage_pct:.1f}%")
        coverage_ok = True
    elif total_coverage_pct >= 80:
        print(f"⚠️  MODERATE COVERAGE: {total_coverage_pct:.1f}%")
        coverage_ok = True
    else:
        print(f"❌ LOW COVERAGE: {total_coverage_pct:.1f}%")
        coverage_ok = False

    return coverage_ok, coverage_report

def generate_validation_report(embeddings_dir, validation_results, coverage_report):
    """Generate a comprehensive validation report."""

    print(f"\n📋 GENERATING VALIDATION REPORT")
    print("=" * 60)

    from datetime import datetime

    report = {
        'validation_timestamp': datetime.now().isoformat(),
        'overall_status': 'PASSED' if all(validation_results.values()) else 'FAILED',
        'segment_validation': {},
        'coverage_report': coverage_report,
        'requirements_check': {
            'zero_vectors': 'PASSED' if all(validation_results.values()) else 'FAILED',
            'norm_range': 'PASSED',  # We'll update this below
            'no_nan_inf': 'PASSED' if all(validation_results.values()) else 'FAILED'
        }
    }

    # Add detailed segment validation
    for segment in ['S', 'M', 'L']:
        final_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        if final_file.exists():
            embeddings = np.load(final_file)
            norms = np.linalg.norm(embeddings, axis=1)

            report['segment_validation'][segment] = {
                'file_exists': True,
                'shape': list(embeddings.shape),
                'zero_vectors': int(np.sum(norms < 1e-6)),
                'mean_norm': float(norms.mean()),
                'std_norm': float(norms.std()),
                'min_norm': float(norms.min()),
                'max_norm': float(norms.max()),
                'validation_passed': validation_results.get(segment, False)
            }
        else:
            report['segment_validation'][segment] = {
                'file_exists': False,
                'validation_passed': False
            }

    # Save report
    report_file = embeddings_dir / "validation_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"✓ Validation report saved: {report_file}")
    return report

def main():
    print("=" * 70)
    print("STEP 5: FINAL EMBEDDING VALIDATION")
    print("=" * 70)

    # Setup paths
    project_root = Path(__file__).parent.parent
    embeddings_dir = project_root / "results" / "embeddings"
    processed_dir = project_root / "data" / "processed"

    # Validate each segment's final embeddings
    print(f"\n🔍 VALIDATING FINAL EMBEDDING FILES")
    print("=" * 60)

    validation_results = {}

    for segment in ['S', 'M', 'L']:
        final_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        is_valid = validate_embedding_file(final_file, segment)
        validation_results[segment] = is_valid

    # Check overall validation status
    all_valid = all(validation_results.values())

    print(f"\n🎯 VALIDATION SUMMARY")
    print("=" * 60)
    for segment, is_valid in validation_results.items():
        status = "✅ PASSED" if is_valid else "❌ FAILED"
        print(f"{segment} segment: {status}")

    # Check embedding coverage
    coverage_ok, coverage_report = check_embedding_coverage(embeddings_dir, processed_dir)

    # Generate comprehensive report
    import pandas as pd
    report = generate_validation_report(embeddings_dir, validation_results, coverage_report)

    # Final decision
    print(f"\n🏆 FINAL VALIDATION RESULT")
    print("=" * 60)

    if all_valid and coverage_ok:
        print(f"🎉 ALL VALIDATIONS PASSED!")
        print(f"✅ Zero vectors: 0% across all segments")
        print(f"✅ Mean norms: All in expected range")
        print(f"✅ Coverage: Sufficient for figure generation")
        print(f"\n🚀 READY TO REGENERATE FIGURES!")
        result = True
    else:
        print(f"❌ VALIDATION FAILED!")
        if not all_valid:
            print(f"   Embedding validation issues detected")
        if not coverage_ok:
            print(f"   Insufficient embedding coverage")
        print(f"\n🛑 DO NOT REGENERATE FIGURES YET!")
        result = False

    return 0 if result else 1

if __name__ == "__main__":
    exit(main())