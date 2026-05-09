"""
Quality control filters for sequences.

Filters applied in order:
1. Length (350-600 aa for Gn)
2. Ambiguous amino acids (max 5%)
3. Internal stop codons
4. Exact duplicates (SHA256)
5. Near-duplicates (≥99% identity)
"""

import hashlib
from typing import Dict, List, Set, Tuple

import pandas as pd


def seq_to_hash(seq_str: str) -> str:
    """Compute SHA256 hash of sequence."""
    return hashlib.sha256(seq_str.encode()).hexdigest()


def check_length(seq: str, min_len: int = 350, max_len: int = 600) -> Tuple[bool, str]:
    """Check if sequence length is within range."""
    seq_len = len(seq)
    if min_len <= seq_len <= max_len:
        return True, ""
    return False, f"length_fail ({seq_len} aa)"


def check_ambiguous(seq: str, max_frac: float = 0.05) -> Tuple[bool, str]:
    """Check if ambiguous AA fraction is acceptable."""
    ambiguous_chars = set("XBZJ")
    ambig_count = sum(1 for c in seq if c in ambiguous_chars)
    ambig_frac = ambig_count / len(seq) if seq else 0

    if ambig_frac <= max_frac:
        return True, ""
    return False, f"ambiguous_aa_fail ({ambig_frac:.2%})"


def check_internal_stops(seq: str) -> Tuple[bool, str]:
    """Check for internal stop codons (*)."""
    if "*" in seq[:-1]:  # Allow * at end only
        return False, "internal_stop_fail"
    return True, ""


def apply_qc_filters(
    records_with_seqs: Dict[str, Dict],
    config_qc: Dict,
) -> Tuple[Dict[str, Dict], Dict[str, List[str]]]:
    """
    Apply QC filters to sequences.

    Returns:
        (passed_records, failed_records_by_reason)
    """
    passed = {}
    failed = {}

    for accession, data in records_with_seqs.items():
        seq = data["seq"]
        qc_flags = []

        # Filter 1: Length
        ok, msg = check_length(seq, config_qc["gn_length_min"], config_qc["gn_length_max"])
        if not ok:
            qc_flags.append(msg)

        # Filter 2: Ambiguous AA
        ok, msg = check_ambiguous(seq, config_qc["max_ambiguous_fraction"])
        if not ok:
            qc_flags.append(msg)

        # Filter 3: Internal stops
        ok, msg = check_internal_stops(seq)
        if not ok:
            qc_flags.append(msg)

        # If failed any filter, skip
        if qc_flags:
            failed[accession] = qc_flags
            continue

        # Passed QC
        data["metadata"]["qc_passed"] = True
        data["metadata"]["qc_flags"] = None
        data["seq_hash"] = seq_to_hash(seq)
        passed[accession] = data

    return passed, failed


def remove_exact_duplicates(
    records: Dict[str, Dict],
) -> Tuple[Dict[str, Dict], Dict[str, str]]:
    """
    Remove exact duplicates (same sequence hash).
    Keep first occurrence (by accession lexicographically).

    Returns:
        (unique_records, duplicates_map {dup_accession: original_accession})
    """
    seq_hash_to_accession = {}
    unique = {}
    duplicates = {}

    # Sort by accession to ensure deterministic "first" selection
    for accession in sorted(records.keys()):
        data = records[accession]
        seq_hash = data.get("seq_hash")

        if seq_hash in seq_hash_to_accession:
            # Duplicate
            original = seq_hash_to_accession[seq_hash]
            duplicates[accession] = original
        else:
            # First occurrence
            seq_hash_to_accession[seq_hash] = accession
            unique[accession] = data

    return unique, duplicates


def compute_pairwise_identity(seq1: str, seq2: str) -> float:
    """
    Compute pairwise sequence identity.
    Simple character-by-character comparison.
    """
    if not seq1 or not seq2:
        return 0.0

    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    identity = matches / max(len(seq1), len(seq2))
    return identity


def remove_near_duplicates(
    records: Dict[str, Dict],
    identity_threshold: float = 0.99,
) -> Tuple[Dict[str, Dict], Dict[str, str]]:
    """
    Remove near-duplicates (≥99% identity).
    Keep first occurrence (by accession lexicographically).

    Args:
        records: dict of {accession: {seq, ...}}
        identity_threshold: min identity to consider duplicate

    Returns:
        (unique_records, duplicates_map {dup_accession: original_accession})
    """
    # Sort accessions for deterministic behavior
    sorted_accessions = sorted(records.keys())
    seq_dict = {acc: records[acc]["seq"] for acc in sorted_accessions}

    kept = set()
    duplicates = {}

    for acc_idx, acc in enumerate(sorted_accessions):
        if acc in duplicates:
            continue

        kept.add(acc)
        seq = seq_dict[acc]

        # Compare with all later sequences
        for other_acc in sorted_accessions[acc_idx + 1 :]:
            if other_acc in duplicates:
                continue

            other_seq = seq_dict[other_acc]
            identity = compute_pairwise_identity(seq, other_seq)

            if identity >= identity_threshold:
                duplicates[other_acc] = acc

    # Build result dict
    unique = {acc: records[acc] for acc in kept}
    return unique, duplicates


def generate_qc_report(
    total_fetched: int,
    qc_passed: int,
    failed_by_reason: Dict[str, List[str]],
    exact_dups: Dict,
    near_dups: Dict,
    metadata_df: pd.DataFrame,
) -> Dict:
    """Generate QC report JSON."""
    species_counts = metadata_df["species_clean"].value_counts().to_dict()
    old_new_counts = metadata_df["old_new_world"].value_counts().to_dict()

    old_count = old_new_counts.get("Old World", 0)
    new_count = old_new_counts.get("New World", 0)
    total_with_world = old_count + new_count

    old_frac = old_count / total_with_world if total_with_world > 0 else 0
    new_frac = new_count / total_with_world if total_with_world > 0 else 0

    # Count failures by reason
    reason_counts = {}
    for reason_list in failed_by_reason.values():
        for reason in reason_list:
            # Extract reason without details
            reason_base = reason.split("(")[0].strip()
            reason_counts[reason_base] = reason_counts.get(reason_base, 0) + 1

    return {
        "total_fetched": total_fetched,
        "passed_qc": qc_passed,
        "failed_qc": total_fetched - qc_passed - len(exact_dups) - len(near_dups),
        "exact_duplicates": len(exact_dups),
        "near_duplicates": len(near_dups),
        "failure_reasons": reason_counts,
        "species_distribution": species_counts,
        "old_world_count": old_count,
        "new_world_count": new_count,
        "old_world_fraction": round(old_frac, 3),
        "new_world_fraction": round(new_frac, 3),
    }
