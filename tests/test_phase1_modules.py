#!/usr/bin/env python3
"""
Smoke tests for Phase 1 modules (no NCBI calls).

Tests:
- Metadata extraction and normalization
- QC filters
- Deduplication
- Manifest generation
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.metadata import (
    parse_country,
    parse_year,
    categorize_host,
    clean_species_name,
)
from data.qc import (
    check_length,
    check_ambiguous,
    check_internal_stops,
    seq_to_hash,
    remove_exact_duplicates,
)


def test_parse_country():
    """Test country parsing."""
    tests = [
        ("Argentina: Patagonia", ("Argentina: Patagonia", "Argentina", "South America", "New World")),
        ("China", ("China", "China", "East Asia", "Old World")),
        ("USA", ("USA", "USA", "North America", "New World")),
    ]

    for input_val, expected in tests:
        result = parse_country(input_val)
        assert result == expected, f"Failed: {input_val} → {result}"
        print(f"✓ parse_country({input_val!r})")


def test_parse_year():
    """Test year parsing."""
    tests = [
        ("2018", 2018),
        ("2018-03", 2018),
        ("2018-03-15", 2018),
        ("199X", None),
        ("", None),
    ]

    for input_val, expected in tests:
        result = parse_year(input_val)
        assert result == expected, f"Failed: {input_val} → {result}"
        print(f"✓ parse_year({input_val!r})")


def test_categorize_host():
    """Test host categorization."""
    tests = [
        ("Homo sapiens", "human"),
        ("Peromyscus maniculatus", "rodent"),
        ("unknown", "unknown"),
        ("", "unknown"),
    ]

    for input_val, expected in tests:
        result = categorize_host(input_val)
        assert result == expected, f"Failed: {input_val} → {result}"
        print(f"✓ categorize_host({input_val!r})")


def test_clean_species_name():
    """Test species name cleaning."""
    tests = [
        ("Andes virus segment M GPC", "Andes_virus"),
        ("Hantaan virus", "Hantaan_virus"),
        ("Unknown species", "Unknown_species"),
    ]

    for input_val, expected in tests:
        result = clean_species_name(input_val)
        assert result == expected, f"Failed: {input_val} → {result}"
        print(f"✓ clean_species_name({input_val!r})")


def test_qc_filters():
    """Test QC filter functions."""
    # Length filter
    ok, msg = check_length("M" * 400)
    assert ok, f"400 aa should pass: {msg}"
    print("✓ check_length(400 aa) passes")

    ok, msg = check_length("M" * 300)
    assert not ok, f"300 aa should fail: {msg}"
    print("✓ check_length(300 aa) fails")

    # Ambiguous AA filter
    ok, msg = check_ambiguous("M" * 390 + "X" * 10)  # 10/400 = 2.5%
    assert ok, f"2.5% ambiguous should pass: {msg}"
    print("✓ check_ambiguous(2.5%) passes")

    ok, msg = check_ambiguous("M" * 380 + "X" * 20)  # 20/400 = 5%
    assert ok, f"5% ambiguous should pass: {msg}"
    print("✓ check_ambiguous(5%) passes")

    ok, msg = check_ambiguous("M" * 370 + "X" * 30)  # 30/400 = 7.5%
    assert not ok, f"7.5% ambiguous should fail: {msg}"
    print("✓ check_ambiguous(7.5%) fails")

    # Stop codon filter
    ok, msg = check_internal_stops("MLKDP*")
    assert ok, f"Terminal stop should pass: {msg}"
    print("✓ check_internal_stops (terminal) passes")

    ok, msg = check_internal_stops("MLK*DP")
    assert not ok, f"Internal stop should fail: {msg}"
    print("✓ check_internal_stops (internal) fails")


def test_hash():
    """Test sequence hashing."""
    seq1 = "MLKDPQR"
    seq2 = "MLKDPQR"
    seq3 = "MLKDPQ"

    hash1 = seq_to_hash(seq1)
    hash2 = seq_to_hash(seq2)
    hash3 = seq_to_hash(seq3)

    assert hash1 == hash2, "Same sequences should hash same"
    assert hash1 != hash3, "Different sequences should hash different"
    print(f"✓ seq_to_hash deterministic and unique")


def test_remove_exact_duplicates():
    """Test exact duplicate removal."""
    records = {
        "seq1": {"seq": "MLKDP", "seq_hash": seq_to_hash("MLKDP")},
        "seq2": {"seq": "MLKDP", "seq_hash": seq_to_hash("MLKDP")},  # duplicate
        "seq3": {"seq": "MQPRL", "seq_hash": seq_to_hash("MQPRL")},
    }

    unique, dups = remove_exact_duplicates(records)

    assert len(unique) == 2, f"Should have 2 unique, got {len(unique)}"
    assert len(dups) == 1, f"Should have 1 duplicate, got {len(dups)}"
    assert "seq1" in unique, "seq1 (first lexicographically) should be kept"
    assert "seq2" in dups, "seq2 should be marked as duplicate"
    print("✓ remove_exact_duplicates works correctly")


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("PHASE 1 MODULE TESTS")
    print("="*60 + "\n")

    try:
        test_parse_country()
        print()
        test_parse_year()
        print()
        test_categorize_host()
        print()
        test_clean_species_name()
        print()
        test_qc_filters()
        print()
        test_hash()
        print()
        test_remove_exact_duplicates()

        print("\n" + "="*60)
        print("✓ ALL TESTS PASSED")
        print("="*60 + "\n")
        return 0

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
