#!/usr/bin/env python3
"""
Step 4: Fixed validation with correct sequence-to-embedding mapping.
"""

import numpy as np
import json
import re
from Bio import SeqIO
from pathlib import Path

def get_segment_from_extraction_method(header):
    """Extract segment using extraction method as ground truth."""
    SEGMENT_MAP = {
        'n_protein_full': 'S',
        'gpc_full': 'M', 'gpc_truncated_1022': 'M', 'gn_domain': 'M', 'gpc_nterm': 'M',
        'nucleotide_tier1_gpc': 'M', 'nucleotide_tier2_gpc': 'M', 'nucleotide_tier2_large_cds': 'M',
        'rdrp_full': 'L',
        'protein_db_direct': None, 'protein_db': None,
    }

    if re.search(r'_[SML]\s', header):
        return re.search(r'_([SML])\s', header).group(1)

    for method, segment in SEGMENT_MAP.items():
        if method in header:
            return segment
    return None

def extract_base_accession(record_id):
    """Fixed accession extraction that handles periods correctly."""
    # Remove segment suffix first (_S, _M, _L)
    if '_' in record_id and record_id.split('_')[-1] in ['S', 'M', 'L']:
        base_with_version = record_id.rsplit('_', 1)[0]
    else:
        base_with_version = record_id

    # Keep the version number (everything before first space, if any)
    return base_with_version.split()[0]

def extract_species_from_header(header):
    """Extract species name from header."""
    match = re.match(r'^[^\s]+\s+([^|]+)', header)
    if match:
        return match.group(1).strip()
    return 'unknown'

def main():
    print("=" * 70)
    print("STEP 4: FIXED VALIDATION - SEQUENCE-TO-EMBEDDING MAPPING")
    print("=" * 70)

    # Load FASTA sequences
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    sequences = []  # List of all sequences with their info
    species_sequence_counts = {}  # species -> {S: count, M: count, L: count}

    print(f"\n📋 Loading FASTA sequences...")
    for record in SeqIO.parse(visualization_file, "fasta"):
        segment = get_segment_from_extraction_method(record.description)
        if segment is None:
            continue

        base_accession = extract_base_accession(record.id)
        species = extract_species_from_header(record.description)

        sequences.append({
            'record_id': record.id,
            'base_accession': base_accession,
            'segment': segment,
            'species': species
        })

        if species not in species_sequence_counts:
            species_sequence_counts[species] = {'S': 0, 'M': 0, 'L': 0}
        species_sequence_counts[species][segment] += 1

    print(f"  Loaded {len(sequences)} sequences from FASTA")

    # Load embedding indices
    print(f"\n📋 Loading embedding indices...")
    embedded_accessions = {'S': set(), 'M': set(), 'L': set()}

    for segment in ['S', 'M', 'L']:
        index_file = Path(f"results/embeddings/embeddings_{segment}_index.json")
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)
            if 'accession_ids' in index:
                embedded_accessions[segment] = set(index['accession_ids'])

    total_embedded_accessions = sum(len(accs) for accs in embedded_accessions.values())
    print(f"  Loaded {total_embedded_accessions} embedding accessions")

    # Map sequences to embeddings
    print(f"\n🔍 Mapping sequences to embeddings...")

    species_embedded_counts = {}  # species -> {S: count, M: count, L: count}
    missing_sequences = []
    total_mapped = 0

    for seq_info in sequences:
        accession = seq_info['base_accession']
        segment = seq_info['segment']
        species = seq_info['species']

        if species not in species_embedded_counts:
            species_embedded_counts[species] = {'S': 0, 'M': 0, 'L': 0}

        if accession in embedded_accessions[segment]:
            species_embedded_counts[species][segment] += 1
            total_mapped += 1
        else:
            missing_sequences.append(seq_info)

    mapping_rate = (total_mapped / len(sequences)) * 100
    print(f"  Mapped {total_mapped}/{len(sequences)} sequences ({mapping_rate:.1f}%)")

    # Check for zero vectors
    print(f"\n🔍 Checking for zero vectors...")
    total_embeddings = 0
    zero_vectors = 0

    for segment in ['S', 'M', 'L']:
        embedding_file = Path(f"results/embeddings/embeddings_{segment}.npy")
        if embedding_file.exists():
            embeddings = np.load(embedding_file)
            segment_zeros = sum(1 for emb in embeddings if np.linalg.norm(emb) < 1e-6)
            zero_vectors += segment_zeros
            total_embeddings += len(embeddings)

    zero_percentage = (zero_vectors / total_embeddings) * 100 if total_embeddings > 0 else 0
    print(f"  Zero vectors: {zero_vectors}/{total_embeddings} ({zero_percentage:.1f}%)")

    # Generate coverage report
    print(f"\n📊 SPECIES COVERAGE REPORT:")
    print(f"{'Species':<15} | {'FASTA':<5} | {'Embedded':<8} | {'Missing':<7} | {'Coverage':<8}")
    print("-" * 65)

    target_species = ['Seoul', 'Andes', 'Puumala', 'Prospect Hill', 'Choclo', 'Sin Nombre']

    all_species_pass = True
    overall_fasta = 0
    overall_embedded = 0

    for species in target_species:
        fasta_count = sum(species_sequence_counts.get(species, {}).values())
        embedded_count = sum(species_embedded_counts.get(species, {}).values())
        missing_count = fasta_count - embedded_count
        coverage_pct = (embedded_count / fasta_count) * 100 if fasta_count > 0 else 0

        # Pass criteria: >95% coverage for species with N≥5
        pass_criteria = coverage_pct > 95.0 if fasta_count >= 5 else True
        status = "✓" if pass_criteria else "❌"

        if not pass_criteria:
            all_species_pass = False

        print(f"{species:<15} | {fasta_count:<5} | {embedded_count:<8} | {missing_count:<7} | {coverage_pct:<6.1f}% {status}")

        overall_fasta += fasta_count
        overall_embedded += embedded_count

    # Other species summary
    other_fasta = 0
    other_embedded = 0
    for species, counts in species_sequence_counts.items():
        if species not in target_species:
            fasta_count = sum(counts.values())
            embedded_count = sum(species_embedded_counts.get(species, {}).values())
            other_fasta += fasta_count
            other_embedded += embedded_count

    if other_fasta > 0:
        other_coverage = (other_embedded / other_fasta) * 100
        print(f"{'Other species':<15} | {other_fasta:<5} | {other_embedded:<8} | {other_fasta-other_embedded:<7} | {other_coverage:<6.1f}%")
        overall_fasta += other_fasta
        overall_embedded += other_embedded

    overall_coverage = (overall_embedded / overall_fasta) * 100 if overall_fasta > 0 else 0

    print("-" * 65)
    print(f"{'TOTAL':<15} | {overall_fasta:<5} | {overall_embedded:<8} | {overall_fasta-overall_embedded:<7} | {overall_coverage:<6.1f}%")

    # Show missing sequence examples
    if missing_sequences:
        print(f"\n❌ Sample missing sequences:")
        missing_by_species = {}
        for seq in missing_sequences:
            species = seq['species']
            if species not in missing_by_species:
                missing_by_species[species] = []
            missing_by_species[species].append(f"{seq['base_accession']}({seq['segment']})")

        for species, missing_list in list(missing_by_species.items())[:5]:
            print(f"  {species}: {missing_list[:3]}")

    # Final assessment
    print(f"\n🎯 FINAL ASSESSMENT:")
    print(f"Overall coverage: {overall_coverage:.1f}%")
    print(f"Zero vectors: {zero_percentage:.1f}%")

    # Pass/Fail criteria
    overall_pass = (overall_coverage > 95.0 and
                   zero_percentage == 0.0 and
                   all_species_pass)

    if overall_pass:
        print(f"\n✅ PASS - Ready for figure generation!")
    else:
        print(f"\n❌ FAIL - Issues detected:")
        if overall_coverage <= 95.0:
            print(f"  - Coverage: {overall_coverage:.1f}% (need >95%)")
        if zero_percentage > 0:
            print(f"  - Zero vectors: {zero_percentage:.1f}% (need 0%)")
        if not all_species_pass:
            print(f"  - Some species below threshold")

    return 0 if overall_pass else 1

if __name__ == "__main__":
    exit(main())