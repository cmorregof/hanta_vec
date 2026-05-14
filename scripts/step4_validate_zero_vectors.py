#!/usr/bin/env python3
"""
Step 4: Validate that zero vectors = 0% after fixed embedding computation.
"""

import numpy as np
import json
from pathlib import Path

def validate_zero_vectors():
    """Check for zero vectors in final embedding files."""

    print("=" * 70)
    print("STEP 4: ZERO VECTOR VALIDATION")
    print("=" * 70)

    segments = ['S', 'M', 'L']
    embeddings_dir = Path("results/embeddings")

    total_embeddings = 0
    total_zero_vectors = 0

    for segment in segments:
        print(f"\n🧬 Validating {segment} segment embeddings...")

        # Load final embeddings
        embedding_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        index_file = embeddings_dir / f"embeddings_{segment}_final_index.json"

        if not embedding_file.exists():
            print(f"  ❌ {segment} embedding file not found: {embedding_file}")
            continue

        if not index_file.exists():
            print(f"  ❌ {segment} index file not found: {index_file}")
            continue

        # Load data
        embeddings = np.load(embedding_file)
        with open(index_file, 'r') as f:
            index = json.load(f)

        # Check dimensions
        n_sequences, embed_dim = embeddings.shape
        print(f"  📊 Loaded: {n_sequences} embeddings × {embed_dim} dimensions")

        # Validate zero vectors
        zero_vectors = 0
        low_norm_vectors = 0
        mean_norms = []

        for i, embedding in enumerate(embeddings):
            norm = np.linalg.norm(embedding)
            mean_norms.append(norm)

            if norm == 0:
                zero_vectors += 1
            elif norm < 1e-6:
                low_norm_vectors += 1

        # Calculate statistics
        zero_percentage = (zero_vectors / n_sequences) * 100
        low_norm_percentage = (low_norm_vectors / n_sequences) * 100
        mean_norm = np.mean(mean_norms)

        print(f"  ✓ Zero vectors: {zero_vectors}/{n_sequences} ({zero_percentage:.1f}%)")
        print(f"  ✓ Low-norm vectors (<1e-6): {low_norm_vectors}/{n_sequences} ({low_norm_percentage:.1f}%)")
        print(f"  ✓ Mean embedding norm: {mean_norm:.6f}")

        # Accumulate totals
        total_embeddings += n_sequences
        total_zero_vectors += zero_vectors

        # Show sample accessions for verification
        sample_accessions = list(index.keys())[:3]
        print(f"  📝 Sample accessions: {sample_accessions}")

    # Final summary
    print(f"\n{'='*50}")
    print("ZERO VECTOR VALIDATION SUMMARY")
    print(f"{'='*50}")

    total_zero_percentage = (total_zero_vectors / total_embeddings) * 100
    print(f"Total embeddings: {total_embeddings}")
    print(f"Total zero vectors: {total_zero_vectors}")
    print(f"Zero vector percentage: {total_zero_percentage:.1f}%")

    if total_zero_percentage == 0:
        print(f"✅ SUCCESS: 0% zero vectors achieved!")
        print(f"Fixed embedding pipeline working correctly")
        print(f"Ready to regenerate figures with high-quality embeddings")
    else:
        print(f"❌ ISSUE: {total_zero_percentage:.1f}% zero vectors detected")
        print(f"Embedding pipeline needs further investigation")

    return total_zero_percentage

def main():
    zero_percentage = validate_zero_vectors()
    return 0 if zero_percentage == 0 else 1

if __name__ == "__main__":
    exit(main())