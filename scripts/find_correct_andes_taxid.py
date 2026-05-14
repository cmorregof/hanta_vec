#!/usr/bin/env python3
"""
Find the correct Andes virus taxid by searching NCBI taxonomy directly.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import setup_ncbi, search_ncbi, load_config
from Bio import Entrez

def search_taxonomy(term):
    """Search NCBI taxonomy database for organism."""
    print(f"Searching taxonomy for: {term}")
    try:
        handle = Entrez.esearch(db="taxonomy", term=term, retmax=10)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            for taxid in record["IdList"]:
                # Get details for each taxid
                handle = Entrez.efetch(db="taxonomy", id=taxid)
                tax_record = Entrez.read(handle)
                handle.close()

                if tax_record:
                    name = tax_record[0].get("ScientificName", "Unknown")
                    lineage = tax_record[0].get("Lineage", "Unknown")
                    print(f"  Taxid {taxid}: {name}")
                    if "hantavirus" in name.lower() or "orthohantavirus" in name.lower():
                        print(f"    ✓ MATCH: {lineage}")
                        return taxid, name
                    else:
                        print(f"    - Not hantavirus: {lineage}")
        else:
            print(f"  No taxonomy results for {term}")

    except Exception as e:
        print(f"  Error searching taxonomy: {e}")

    return None, None

def test_organism_search_variants():
    """Try various organism search terms to find Andes virus."""
    search_terms = [
        "Andes virus",
        "Andes hantavirus",
        "Orthohantavirus andesense",
        "Andes orthohantavirus",
        '"Andes virus"[Organism]',
    ]

    results = {}

    for term in search_terms:
        print(f"\n{'-'*40}")
        print(f"Testing: {term}")
        print(f"{'-'*40}")

        try:
            query = term
            if not ("[" in query and "]" in query):
                query = f'"{term}"'

            nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=10)
            print(f"  Nucleotide results: {len(nucleotide_ids)}")

            if nucleotide_ids:
                # Fetch a few to check organisms
                from data.fetch import fetch_batch
                organisms = set()
                for records, failed in fetch_batch(nucleotide_ids[:3], db="nucleotide", batch_size=3, sleep_time=0.5):
                    for record in records:
                        organism = record.annotations.get('organism', 'Unknown')
                        organisms.add(organism)
                        print(f"    {organism}")

                results[term] = list(organisms)

        except Exception as e:
            print(f"    Error: {e}")

    return results

def main():
    print("="*70)
    print("FINDING CORRECT ANDES VIRUS TAXID")
    print("="*70)

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # Try taxonomy search first
    print("\n1. TAXONOMY DATABASE SEARCH:")
    taxid, name = search_taxonomy("Andes virus")

    if not taxid:
        taxid, name = search_taxonomy("Andes hantavirus")

    if not taxid:
        taxid, name = search_taxonomy("Orthohantavirus andesense")

    # Try organism searches
    print("\n2. ORGANISM SEARCH VARIANTS:")
    organism_results = test_organism_search_variants()

    print("\n" + "="*70)
    print("SUMMARY:")
    print("="*70)

    if taxid:
        print(f"✅ Found potential correct taxid: {taxid}")
        print(f"   Organism: {name}")

        # Verify this taxid works
        print(f"\nVerifying taxid {taxid}...")
        query = f"txid{taxid}[Organism:exp]"
        test_ids = search_ncbi(query, db="nucleotide", retmax=5)
        print(f"   Test query results: {len(test_ids)}")

    else:
        print("❌ Could not find correct Andes virus taxid in taxonomy")

    print(f"\nOrganism search summary:")
    for term, orgs in organism_results.items():
        if orgs:
            print(f"  {term}: {', '.join(orgs)}")

if __name__ == "__main__":
    main()