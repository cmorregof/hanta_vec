"""
Baseline: Simple pairwise sequence identity (no learned representation).

Serves as control for evaluating ESM-2 embeddings.
Computes percent identity from global pairwise alignment.
"""

import numpy as np
from Bio import Align


def pairwise_identity(seq1: str, seq2: str, aligner=None) -> float:
    """
    Compute percent sequence identity via global pairwise alignment.

    Args:
        seq1: First protein sequence (str)
        seq2: Second protein sequence (str)
        aligner: Optional Bio.Align.PairwiseAligner instance

    Returns:
        Identity as fraction [0.0, 1.0]

    Example:
        >>> pairwise_identity("MKFTK", "MKFTK")
        1.0
        >>> pairwise_identity("MKFTK", "MKFSX")
        0.8
    """
    if aligner is None:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0

    try:
        alignments = aligner.align(seq1.upper(), seq2.upper())
        if not alignments:
            return 0.0

        best = alignments[0]
        matches = sum(1 for a, b in zip(best.target, best.query) if a == b and a != '-')
        total = max(len(seq1), len(seq2))

        return matches / total if total > 0 else 0.0

    except Exception:
        return 0.0


def pairwise_identity_matrix(sequences: dict, sample_size: int = None) -> tuple:
    """
    Compute all-pairs sequence identity matrix.

    Args:
        sequences: {accession: sequence_str} dict
        sample_size: If set, randomly sample this many sequences (for speed)

    Returns:
        (identities, accessions) where:
        - identities: (N, N) symmetric matrix of percent identities
        - accessions: List of accession IDs in matrix order
    """
    accessions = list(sequences.keys())

    if sample_size and len(accessions) > sample_size:
        np.random.seed(42)
        idx = np.random.choice(len(accessions), sample_size, replace=False)
        accessions = [accessions[i] for i in sorted(idx)]

    n = len(accessions)
    identities = np.zeros((n, n))

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'

    for i in range(n):
        for j in range(i + 1, n):
            seq_i = sequences[accessions[i]]
            seq_j = sequences[accessions[j]]

            identity = pairwise_identity(seq_i, seq_j, aligner)
            identities[i, j] = identity
            identities[j, i] = identity

    np.fill_diagonal(identities, 1.0)

    return identities, accessions
