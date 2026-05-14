#!/usr/bin/env python3
"""
Fetch additional New World sequences from ViPR (Virus Pathogen Resource).
Target species: Andes, Sin Nombre, Black Creek Canal, Caño Delgadito, Choclo
Use ViPR REST API and deduplicate at 0.99 threshold.
"""

import requests
import json
import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from typing import Dict, List
import time

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def search_vipr_sequences(species_name: str, segment: str = None) -> List[Dict]:
    """Search ViPR for hantavirus sequences."""

    # ViPR search API endpoint
    base_url = "https://www.viprbrc.org/brc/api"

    # Try different search approaches
    search_terms = [
        f"{species_name} virus",
        f"{species_name} hantavirus",
        f"Orthohantavirus {species_name.lower()}",
    ]

    all_results = []

    for term in search_terms:
        try:
            print(f"    Searching ViPR for: '{term}'")

            # Search parameters
            params = {
                'q': term,
                'family': 'Hantaviridae',
                'wt': 'json',
                'rows': 1000,  # Max results
                'fl': 'accession,organism,country,collection_date,segment,host,length,definition'
            }

            response = requests.get(f"{base_url}/search", params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()

                if 'response' in data and 'docs' in data['response']:
                    results = data['response']['docs']
                    print(f"      Found {len(results)} results")
                    all_results.extend(results)
                else:
                    print(f"      No results structure found")

            else:
                print(f"      HTTP error: {response.status_code}")

        except Exception as e:
            print(f"      Error: {e}")

        # Rate limiting
        time.sleep(1)

    # Deduplicate by accession
    seen_accessions = set()
    unique_results = []

    for result in all_results:
        acc = result.get('accession', '')
        if acc and acc not in seen_accessions:
            seen_accessions.add(acc)
            unique_results.append(result)

    print(f"    Unique results after deduplication: {len(unique_results)}")
    return unique_results

def fetch_sequences_from_vipr(accessions: List[str]) -> Dict[str, str]:
    """Fetch actual sequence data from ViPR/NCBI."""

    print(f"  Fetching sequences for {len(accessions)} accessions...")

    # For now, use NCBI since ViPR sequence API may be complex
    # In production, would implement ViPR sequence API

    sequences = {}

    # Use the existing NCBI fetch infrastructure
    sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
    from data.fetch import fetch_batch, setup_ncbi, load_config

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # Batch fetch from NCBI
    batch_size = 50
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]

        try:
            for records, failed in fetch_batch(
                batch,
                db="nucleotide",
                batch_size=len(batch),
                sleep_time=1
            ):
                for record in records:
                    sequences[record.id] = {
                        'sequence': str(record.seq),
                        'description': record.description,
                        'organism': record.annotations.get('organism', 'Unknown')
                    }

        except Exception as e:
            print(f"    Batch fetch error: {e}")

    print(f"  Successfully fetched {len(sequences)} sequences")
    return sequences

def deduplicate_sequences(new_sequences: Dict, existing_sequences: Dict, threshold: float = 0.99) -> Dict:
    """Deduplicate new sequences against existing at given threshold."""

    print(f"  Deduplicating {len(new_sequences)} new vs {len(existing_sequences)} existing (threshold={threshold})...")

    # Simple sequence identity check (could be improved with alignment)
    deduplicated = {}
    duplicate_count = 0

    for new_acc, new_data in new_sequences.items():
        new_seq = new_data['sequence']
        is_duplicate = False

        # Check against existing sequences
        for existing_seq in existing_sequences.values():
            if isinstance(existing_seq, str):
                existing_seq_str = existing_seq
            else:
                existing_seq_str = str(existing_seq)

            # Simple identity check
            if len(new_seq) > 0 and len(existing_seq_str) > 0:
                identity = sum(a == b for a, b in zip(new_seq, existing_seq_str)) / max(len(new_seq), len(existing_seq_str))

                if identity >= threshold:
                    is_duplicate = True
                    duplicate_count += 1
                    break

        if not is_duplicate:
            deduplicated[new_acc] = new_data

    print(f"  Removed {duplicate_count} duplicates, kept {len(deduplicated)} new sequences")
    return deduplicated

def main():
    print("="*70)
    print("ViPR FETCH FOR NEW WORLD HANTAVIRUS SEQUENCES")
    print("="*70)

    # Target New World species
    target_species = [
        "Andes",
        "Sin Nombre",
        "Black Creek Canal",
        "Caño Delgadito",
        "Choclo"
    ]

    print(f"Targeting {len(target_species)} New World species...")

    # Load existing sequences for deduplication
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    existing_sequences = {}

    for segment in ['S', 'M', 'L']:
        fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
        if fasta_path.exists():
            for record in SeqIO.parse(fasta_path, "fasta"):
                existing_sequences[record.id] = str(record.seq)

    print(f"Loaded {len(existing_sequences)} existing sequences for deduplication")

    # Collect results
    vipr_results = {}
    total_found = 0

    for species in target_species:
        print(f"\n{'-'*50}")
        print(f"SPECIES: {species}")
        print(f"{'-'*50}")

        # Search ViPR
        results = search_vipr_sequences(species)

        if results:
            # Extract accessions
            accessions = []
            for result in results:
                acc = result.get('accession', '')
                if acc and acc not in existing_sequences:
                    accessions.append(acc)

            print(f"  New accessions to fetch: {len(accessions)}")

            if accessions:
                # Fetch sequences
                sequences = fetch_sequences_from_vipr(accessions[:100])  # Limit to 100 per species

                # Deduplicate
                deduplicated = deduplicate_sequences(sequences, existing_sequences, threshold=0.99)

                if deduplicated:
                    vipr_results[species] = {
                        'sequences': deduplicated,
                        'metadata': results
                    }
                    total_found += len(deduplicated)
                    print(f"  ✓ Added {len(deduplicated)} unique sequences for {species}")
                else:
                    print(f"  No unique sequences after deduplication")
            else:
                print(f"  No new accessions to fetch")
        else:
            print(f"  ❌ No ViPR results for {species}")

    print(f"\n{'='*70}")
    print("ViPR FETCH SUMMARY")
    print(f"{'='*70}")
    print(f"Total new sequences found: {total_found}")

    if vipr_results:
        print(f"\nPer-species breakdown:")
        for species, data in vipr_results.items():
            count = len(data['sequences'])
            print(f"  {species}: {count} sequences")

        # Save results for integration
        output_dir = Path(__file__).parent.parent / "data" / "raw" / "vipr_new_world"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save detailed results
        with open(output_dir / "vipr_results.json", "w") as f:
            # Make serializable
            serializable_results = {}
            for species, data in vipr_results.items():
                serializable_results[species] = {
                    'sequences': data['sequences'],
                    'count': len(data['sequences'])
                }
            json.dump(serializable_results, f, indent=2)

        print(f"\n✓ Results saved to {output_dir}")

        if total_found >= 10:
            print(f"✅ SUCCESS: Retrieved {total_found} new sequences (≥10)")
            return 0
        else:
            print(f"⚠ PARTIAL: Retrieved {total_found} new sequences (<10)")
            return 1
    else:
        print("❌ No new sequences retrieved from ViPR")
        return 1

if __name__ == "__main__":
    exit(main())