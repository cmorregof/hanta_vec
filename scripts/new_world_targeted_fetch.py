#!/usr/bin/env python3
"""
Enhanced New World species targeted fetch with corrected taxids.
Focus on increasing New World representation in the dataset.

Uses verified taxids and broad search strategies to maximize yield.
"""

import json
import sys
from pathlib import Path
from typing import Dict
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import (
    setup_ncbi, search_ncbi, fetch_batch, load_config,
    extract_s_segment_protein, extract_m_segment_protein, extract_l_segment_protein,
    extract_isolate_strain, classify_species_from_description
)

def simple_sequence_qc(sequence: str, min_length: int = 50) -> bool:
    """Simple QC check for sequence quality."""
    if len(sequence) < min_length:
        return False

    # Check for excessive ambiguous amino acids
    ambiguous_count = sum(1 for aa in sequence if aa in 'XBZ*')
    ambiguous_pct = ambiguous_count / len(sequence) * 100

    return ambiguous_pct <= 10  # Allow up to 10% ambiguous

def main():
    print("="*70)
    print("NEW WORLD HANTAVIRUS TARGETED FETCH (ENHANCED)")
    print("="*70)

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # New World species with verified/corrected taxids
    new_world_species = {
        "Sin Nombre": "3052499",  # CORRECTED from 64008
        "Andes": "1980519",      # From config
        "Maporal": "64005",      # From config
        # Add more New World species with organism name searches
        "Black Creek Canal": None,
        "Bayou": None,
        "Monongahela": None,
        "New York virus": None,  # Found in Sin Nombre search
        "Blue River": None,      # Found in Sin Nombre search
    }

    print(f"\nTargeting {len(new_world_species)} New World species...")

    # Results organized by segment
    all_results = {"S": {}, "M": {}, "L": {}}
    species_stats = {}

    extractors = {
        'S': extract_s_segment_protein,
        'M': extract_m_segment_protein,
        'L': extract_l_segment_protein
    }

    for species_name, taxid in new_world_species.items():
        print(f"\n{'-'*50}")
        print(f"SPECIES: {species_name}")
        print(f"{'TAXID: ' + taxid if taxid else 'ORGANISM SEARCH'}")
        print(f"{'-'*50}")

        try:
            # Search strategy
            if taxid:
                # Use taxid search
                query = f"txid{taxid}[Organism:exp]"
                search_method = f"taxid {taxid}"
            else:
                # Use organism name search
                query = f'"{species_name}" AND hantavirus'
                search_method = f"organism name"

            print(f"Query: {query}")
            print(f"Search method: {search_method}")

            nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=1000)
            print(f"  Found: {len(nucleotide_ids)} nucleotide records")

            if not nucleotide_ids:
                print(f"  ❌ No records found for {species_name}")
                species_stats[species_name] = {"records": 0, "segments": {"S": 0, "M": 0, "L": 0}}
                continue

            # Download records
            print(f"  Downloading {len(nucleotide_ids)} records...")
            species_records = {}

            for records, failed in fetch_batch(
                nucleotide_ids,
                db="nucleotide",
                batch_size=config["ncbi"]["batch_size"],
                sleep_time=config["ncbi"]["sleep_between_batches"],
            ):
                for record in records:
                    species_records[record.id] = record

            print(f"  ✓ Downloaded {len(species_records)} records")

            # Analyze organisms in this search
            organism_counts = {}
            for record in species_records.values():
                organism = record.annotations.get('organism', 'Unknown')
                organism_counts[organism] = organism_counts.get(organism, 0) + 1

            print(f"  Organisms found:")
            for org, count in sorted(organism_counts.items()):
                print(f"    {org}: {count}")

            # Extract segments
            segment_counts = {"S": 0, "M": 0, "L": 0}

            for record_id, record in species_records.items():
                # Classify species from description
                species_classified, clade = classify_species_from_description(record.description)

                for segment, extractor in extractors.items():
                    protein_data = extractor(record)

                    if protein_data:
                        # Basic QC
                        min_length = {"S": 200, "M": 400, "L": 1000}.get(segment, 100)
                        qc_passed = simple_sequence_qc(protein_data['sequence'], min_length)

                        if qc_passed:
                            isolate, strain = extract_isolate_strain(record)

                            # Store with full metadata
                            result_data = {
                                'accession': record.id,
                                'organism': record.annotations.get('organism', 'Unknown'),
                                'description': record.description,
                                'species': species_classified,
                                'clade': clade,
                                'isolate': isolate,
                                'strain': strain,
                                'sequence': protein_data['sequence'],
                                'method': protein_data['method'],
                                'original_length': protein_data['original_length'],
                                'search_species': species_name,
                                'search_method': search_method
                            }

                            all_results[segment][record.id] = result_data
                            segment_counts[segment] += 1

            species_stats[species_name] = {
                "records": len(species_records),
                "segments": segment_counts,
                "total_extracted": sum(segment_counts.values())
            }

            print(f"  Extracted: S={segment_counts['S']}, M={segment_counts['M']}, L={segment_counts['L']}")

        except Exception as e:
            print(f"  ERROR for {species_name}: {e}")
            species_stats[species_name] = {"records": 0, "segments": {"S": 0, "M": 0, "L": 0}, "error": str(e)}

    # Summary report
    print(f"\n{'='*70}")
    print("NEW WORLD TARGETED FETCH SUMMARY")
    print(f"{'='*70}")

    total_by_segment = {seg: len(results) for seg, results in all_results.items()}
    total_sequences = sum(total_by_segment.values())

    print(f"Total sequences extracted: {total_sequences}")
    print(f"  S-segment: {total_by_segment['S']}")
    print(f"  M-segment: {total_by_segment['M']}")
    print(f"  L-segment: {total_by_segment['L']}")

    print(f"\nPer-species breakdown:")
    for species, stats in species_stats.items():
        if "error" in stats:
            print(f"  {species}: ERROR - {stats['error']}")
        else:
            total = stats.get('total_extracted', sum(stats['segments'].values()))
            print(f"  {species}: {total} total (S={stats['segments']['S']}, M={stats['segments']['M']}, L={stats['segments']['L']})")

    # New World species classification check
    new_world_count = 0
    species_breakdown = {}

    for segment_results in all_results.values():
        for result in segment_results.values():
            species = result['species']
            if result['clade'] == 'New World' or 'New World' in result['clade']:
                new_world_count += 1

            species_breakdown[species] = species_breakdown.get(species, 0) + 1

    print(f"\nNew World sequences: {new_world_count}/{total_sequences}")
    print(f"Species classification breakdown:")
    for species, count in sorted(species_breakdown.items()):
        print(f"  {species}: {count}")

    # Save results for potential integration
    output_dir = Path(__file__).parent.parent / "data" / "raw" / "new_world_targeted"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save detailed results
    with open(output_dir / "new_world_results.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Save summary statistics
    summary = {
        "total_sequences": total_sequences,
        "by_segment": total_by_segment,
        "new_world_count": new_world_count,
        "species_stats": species_stats,
        "species_breakdown": species_breakdown
    }

    with open(output_dir / "fetch_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n✓ Results saved to {output_dir}")

    if total_sequences >= 50:  # Success threshold
        print(f"✅ SUCCESS: Retrieved {total_sequences} sequences (≥50)")
        return 0
    else:
        print(f"⚠ PARTIAL: Retrieved {total_sequences} sequences (<50)")
        return 1

if __name__ == "__main__":
    sys.exit(main())