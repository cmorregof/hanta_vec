#!/usr/bin/env python3
"""
Verify Andes virus taxid and check why sequences were dropped.
Test both current taxid (188539) and user-suggested correct taxid (290432).
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import setup_ncbi, search_ncbi, load_config

def test_taxid(taxid, description):
    """Test what organisms a taxid returns."""
    print(f"\n{'-'*50}")
    print(f"TESTING TAXID {taxid} ({description})")
    print(f"{'-'*50}")

    try:
        query = f"txid{taxid}[Organism:exp]"
        nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=20)
        print(f"Query: {query}")
        print(f"Results: {len(nucleotide_ids)} nucleotide records")

        if nucleotide_ids:
            # Fetch a few records to check organisms
            from data.fetch import fetch_batch
            for records, failed in fetch_batch(nucleotide_ids[:5], db="nucleotide", batch_size=5, sleep_time=0.5):
                print(f"\nSample organisms found:")
                for record in records:
                    organism = record.annotations.get('organism', 'Unknown')
                    desc = record.description[:80] + "..." if len(record.description) > 80 else record.description
                    print(f"  {organism}")
                    print(f"    {record.id}: {desc}")
                break
        else:
            print("  ❌ No records found")

        return len(nucleotide_ids) > 0

    except Exception as e:
        print(f"  ERROR: {e}")
        return False

def main():
    print("="*70)
    print("ANDES VIRUS TAXID VERIFICATION")
    print("="*70)

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # Test current taxid in our data
    current_works = test_taxid("188539", "Current in our metadata")

    # Test problematic taxid mentioned earlier
    problem_works = test_taxid("1980519", "Problematic (CCHF)")

    # Test user-suggested correct taxid
    correct_works = test_taxid("290432", "User-suggested correct")

    # Also try direct organism search
    print(f"\n{'-'*50}")
    print("DIRECT ANDES VIRUS ORGANISM SEARCH")
    print(f"{'-'*50}")

    try:
        query = '"Andes virus"'
        nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=10)
        print(f"Query: {query}")
        print(f"Results: {len(nucleotide_ids)} records")

        if nucleotide_ids:
            from data.fetch import fetch_batch
            for records, failed in fetch_batch(nucleotide_ids[:3], db="nucleotide", batch_size=3, sleep_time=0.5):
                print(f"\nDirect search organisms:")
                for record in records:
                    organism = record.annotations.get('organism', 'Unknown')
                    print(f"  {organism} - {record.id}")
                break

    except Exception as e:
        print(f"ERROR in direct search: {e}")

    print(f"\n{'='*70}")
    print("RECOMMENDATIONS:")
    print(f"{'='*70}")

    if current_works:
        print("✓ Current taxid 188539 appears to work")
    else:
        print("❌ Current taxid 188539 is broken")

    if correct_works:
        print("✓ User-suggested taxid 290432 appears to work")
        print("  → Should use this for Andes virus fetch")
    else:
        print("❌ User-suggested taxid 290432 is also problematic")

    if problem_works:
        print("⚠️ Problematic taxid 1980519 still returns results")
    else:
        print("✓ Problematic taxid 1980519 confirmed broken")

if __name__ == "__main__":
    main()