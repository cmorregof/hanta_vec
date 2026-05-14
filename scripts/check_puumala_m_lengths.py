#!/usr/bin/env python3
"""
Check Puumala M-segment sequence length distribution.
"""

from Bio import SeqIO
from pathlib import Path
import collections

def analyze_puumala_m_lengths():
    """Analyze length distribution of Puumala M-segment sequences."""

    print("Checking Puumala M-segment sequence lengths...")

    m_file = Path("data/processed/M_sequences_level1.fasta")
    if not m_file.exists():
        print(f"❌ M-segment file not found: {m_file}")
        return

    # Collect Puumala sequence data
    puumala_lengths = []
    puumala_methods = []
    other_lengths = []

    for record in SeqIO.parse(m_file, "fasta"):
        seq_len = len(record.seq)
        description = record.description

        if "puumala" in description.lower():
            puumala_lengths.append(seq_len)
            # Extract method
            if "|" in description:
                parts = description.split("|")
                method = parts[-1].strip() if len(parts) > 1 else "unknown"
            else:
                method = "unknown"
            puumala_methods.append(method)
        else:
            other_lengths.append(seq_len)

    print(f"\nPuumala M-segment sequences: {len(puumala_lengths)}")
    print(f"Other M-segment sequences: {len(other_lengths)}")

    if len(puumala_lengths) == 0:
        print("No Puumala sequences found in M-segment file")
        return

    # Length distribution
    print(f"\nPuumala M-segment length distribution:")
    print(f"Min length: {min(puumala_lengths)}")
    print(f"Max length: {max(puumala_lengths)}")
    print(f"Mean length: {sum(puumala_lengths)/len(puumala_lengths):.1f}")

    # Binned distribution
    ranges = [
        ("400-500", 400, 500),
        ("500-700", 500, 700),
        ("700-900", 700, 900),
        ("900-1100", 900, 1100),
        ("1100-1300", 1100, 1300),
        ("1300+", 1300, 9999)
    ]

    print(f"\nLength ranges:")
    for label, min_len, max_len in ranges:
        if max_len == 9999:
            count = sum(1 for length in puumala_lengths if length >= min_len)
        else:
            count = sum(1 for length in puumala_lengths if min_len <= length < max_len)
        pct = count / len(puumala_lengths) * 100 if len(puumala_lengths) > 0 else 0
        print(f"  {label:8} aa: {count:3d} sequences ({pct:4.1f}%)")

    # Method distribution
    print(f"\nExtraction methods:")
    method_counts = collections.Counter(puumala_methods)
    for method, count in method_counts.most_common():
        pct = count / len(puumala_methods) * 100
        print(f"  {method:30} {count:3d} ({pct:4.1f}%)")

    # Sample sequences at different length ranges
    print(f"\nSample sequences by length:")
    samples = [
        ("Short (<600aa)", [l for l in puumala_lengths if l < 600]),
        ("Medium (600-1000aa)", [l for l in puumala_lengths if 600 <= l < 1000]),
        ("Long (1000+aa)", [l for l in puumala_lengths if l >= 1000])
    ]

    for label, lengths in samples:
        if lengths:
            print(f"  {label}: {len(lengths)} sequences, lengths: {sorted(lengths)[:5]}")

def main():
    print("=" * 60)
    print("CHECKING PUUMALA M-SEGMENT LENGTH DISTRIBUTION")
    print("=" * 60)

    analyze_puumala_m_lengths()
    return 0

if __name__ == "__main__":
    exit(main())