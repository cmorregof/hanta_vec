#!/usr/bin/env python3
"""
Targeted fetch for Sin Nombre virus using the correct taxid 3052499.
This script will test the corrected taxid and report raw hits before QC.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import (
    setup_ncbi, search_ncbi, fetch_batch, load_config,
    extract_s_segment_protein, extract_m_segment_protein, extract_l_segment_protein,
    extract_isolate_strain, classify_species_from_description
)

def main():
    print("="*60)
    print("TARGETED SIN NOMBRE VIRUS FETCH")
    print("Testing corrected taxid 3052499")
    print("="*60)

    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    config = load_config(config_path)

    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    # Use the correct Sin Nombre taxid
    taxid = "3052499"

    print(f"\n1. SEARCHING NCBI for taxid {taxid}...")

    # Search nucleotide database
    query = f"txid{taxid}[Organism:exp]"
    print(f"Query: {query}")

    try:
        nucleotide_ids = search_ncbi(query, db="nucleotide", retmax=1000)
        print(f"  Nucleotide database: {len(nucleotide_ids)} records")

        protein_ids = search_ncbi(query, db="protein", retmax=1000)
        print(f"  Protein database: {len(protein_ids)} records")

        # Focus on nucleotide records for segment extraction
        if not nucleotide_ids:
            print("  ❌ No nucleotide records found!")
            return 1

        print(f"\n2. DOWNLOADING {len(nucleotide_ids)} nucleotide records...")

        all_records = {}

        for records, failed in fetch_batch(
            nucleotide_ids,
            db="nucleotide",
            batch_size=config["ncbi"]["batch_size"],
            sleep_time=config["ncbi"]["sleep_between_batches"],
        ):
            for record in records:
                all_records[record.id] = record

        print(f"  ✓ Downloaded {len(all_records)} records successfully")

        print("\n3. ANALYZING DOWNLOADED RECORDS...")

        # Check what we actually got
        organism_counts = {}
        for record in all_records.values():
            organism = record.annotations.get('organism', 'Unknown')
            organism_counts[organism] = organism_counts.get(organism, 0) + 1

        print("  Organisms found:")
        for org, count in sorted(organism_counts.items()):
            print(f"    {org}: {count} records")

        print("\n4. TESTING SEGMENT EXTRACTION...")

        extractors = {
            'S': extract_s_segment_protein,
            'M': extract_m_segment_protein,
            'L': extract_l_segment_protein
        }

        segment_results = {seg: [] for seg in ['S', 'M', 'L']}

        for record_id, record in all_records.items():
            # Classify species from description
            species, clade = classify_species_from_description(record.description)

            for segment, extractor in extractors.items():
                protein_data = extractor(record)
                if protein_data:
                    isolate, strain = extract_isolate_strain(record)
                    segment_results[segment].append({
                        'accession': record.id,
                        'organism': record.annotations.get('organism', 'Unknown'),
                        'description': record.description[:100] + "...",
                        'species_classified': species,
                        'isolate': isolate,
                        'strain': strain,
                        'sequence_length': len(protein_data['sequence']),
                        'method': protein_data['method']
                    })

        print("\n5. EXTRACTION RESULTS:")
        total_extracted = 0
        for segment in ['S', 'M', 'L']:
            count = len(segment_results[segment])
            total_extracted += count
            print(f"  {segment}-segment: {count} proteins extracted")

            if count > 0:
                print(f"    Sample records:")
                for i, result in enumerate(segment_results[segment][:3]):  # Show first 3
                    print(f"      {i+1}. {result['accession']} - {result['species_classified']} - {result['sequence_length']}aa")

        print(f"\n6. SUMMARY:")
        print(f"  Total records downloaded: {len(all_records)}")
        print(f"  Total proteins extracted: {total_extracted}")
        print(f"  Sin Nombre sequences: {sum(1 for seg_results in segment_results.values() for r in seg_results if 'Sin Nombre' in r['species_classified'])}")

        # Check for expected organisms
        expected_organisms = [
            "Orthohantavirus sinnombreense",
            "Sin Nombre virus",
            "Sin nombre virus"
        ]

        found_expected = any(any(exp in org for exp in expected_organisms)
                           for org in organism_counts.keys())

        if found_expected:
            print(f"  ✅ Found expected Sin Nombre organisms")
        else:
            print(f"  ❌ No expected Sin Nombre organisms found")
            print(f"  Organisms retrieved: {list(organism_counts.keys())}")

        return 0

    except Exception as e:
        print(f"  ERROR: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())