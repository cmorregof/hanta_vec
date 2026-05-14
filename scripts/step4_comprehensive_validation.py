#!/usr/bin/env python3
"""
Step 4: Comprehensive validation - every FASTA accession has embedding, coverage per species.
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

def extract_species_from_header(header):
    """Extract species name from header."""
    match = re.match(r'^[^\s]+\s+([^|]+)', header)
    if match:
        return match.group(1).strip()
    return 'unknown'

def load_embedding_indices():
    """Load all embedding indices."""
    embeddings_dir = Path("results/embeddings")

    embedded_accessions = {'S': set(), 'M': set(), 'L': set()}

    for segment in ['S', 'M', 'L']:
        index_file = embeddings_dir / f"embeddings_{segment}_index.json"
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)

            if 'accession_ids' in index:
                embedded_accessions[segment] = set(index['accession_ids'])
            else:
                print(f"  Warning: No accession_ids found in {segment} index")

    return embedded_accessions

def validate_zero_vectors():
    """Check for zero vectors in all embeddings."""
    embeddings_dir = Path("results/embeddings")

    total_embeddings = 0
    total_zero_vectors = 0

    for segment in ['S', 'M', 'L']:
        embedding_file = embeddings_dir / f"embeddings_{segment}.npy"
        if embedding_file.exists():
            embeddings = np.load(embedding_file)
            n_sequences = embeddings.shape[0]

            zero_count = 0
            for embedding in embeddings:
                if np.linalg.norm(embedding) < 1e-6:
                    zero_count += 1

            total_embeddings += n_sequences
            total_zero_vectors += zero_count

            print(f"  {segment}: {zero_count}/{n_sequences} zero vectors")

    zero_percentage = (total_zero_vectors / total_embeddings) * 100 if total_embeddings > 0 else 0
    print(f"  Overall: {total_zero_vectors}/{total_embeddings} ({zero_percentage:.1f}%) zero vectors")

    return zero_percentage

def main():
    print("=" * 70)
    print("STEP 4: COMPREHENSIVE VALIDATION")
    print("=" * 70)

    # Load FASTA sequences
    print("\n📋 Loading FASTA sequences...")
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    fasta_sequences = {}  # accession -> {segment, species, etc.}
    species_segment_counts = {}  # species -> {S: count, M: count, L: count}

    for record in SeqIO.parse(visualization_file, "fasta"):
        segment = get_segment_from_extraction_method(record.description)
        if segment is None:
            continue

        base_accession = record.id.split('_')[0] if '_' in record.id else record.id
        species = extract_species_from_header(record.description)

        fasta_sequences[base_accession] = {
            'segment': segment,
            'species': species,
            'record_id': record.id
        }

        if species not in species_segment_counts:
            species_segment_counts[species] = {'S': 0, 'M': 0, 'L': 0}
        species_segment_counts[species][segment] += 1

    print(f"  Loaded {len(fasta_sequences)} sequences from FASTA")

    # Load embedding indices
    print("\n📋 Loading embedding indices...")
    embedded_accessions = load_embedding_indices()

    total_embedded = sum(len(accs) for accs in embedded_accessions.values())
    print(f"  Loaded {total_embedded} embedding accessions")

    # Validate coverage
    print("\n🔍 Validating coverage...")

    species_embedded_counts = {}  # species -> {S: count, M: count, L: count}
    total_missing = 0
    missing_by_species = {}

    for accession, info in fasta_sequences.items():
        segment = info['segment']
        species = info['species']

        if species not in species_embedded_counts:
            species_embedded_counts[species] = {'S': 0, 'M': 0, 'L': 0}

        if accession in embedded_accessions[segment]:
            species_embedded_counts[species][segment] += 1
        else:
            total_missing += 1
            if species not in missing_by_species:
                missing_by_species[species] = []
            missing_by_species[species].append(f"{accession}({segment})")

    # Validate zero vectors
    print("\n🔍 Validating zero vectors...")
    zero_percentage = validate_zero_vectors()

    # Generate report
    print(f"\n📊 COVERAGE VALIDATION REPORT:")
    print(f"{'Species':<15} | {'FASTA':<5} | {'Embedded':<8} | {'Missing':<7} | {'Coverage':<8}")
    print("-" * 65)

    target_species = ['Seoul', 'Andes', 'Puumala', 'Prospect Hill', 'Choclo', 'Sin Nombre']

    all_species_pass = True
    total_fasta_seqs = 0
    total_embedded_seqs = 0

    for species in target_species:
        fasta_count = sum(species_segment_counts.get(species, {}).values())
        embedded_count = sum(species_embedded_counts.get(species, {}).values())
        missing_count = fasta_count - embedded_count
        coverage_pct = (embedded_count / fasta_count) * 100 if fasta_count > 0 else 0

        # Pass criteria: >95% coverage for species with N≥5
        pass_criteria = coverage_pct > 95.0 if fasta_count >= 5 else True
        status = "✓" if pass_criteria else "❌"

        if not pass_criteria:
            all_species_pass = False

        print(f"{species:<15} | {fasta_count:<5} | {embedded_count:<8} | {missing_count:<7} | {coverage_pct:<6.1f}% {status}")

        total_fasta_seqs += fasta_count
        total_embedded_seqs += embedded_count

    # Show other species
    other_fasta = 0
    other_embedded = 0
    for species, counts in species_segment_counts.items():
        if species not in target_species:
            fasta_count = sum(counts.values())
            embedded_count = sum(species_embedded_counts.get(species, {}).values())
            other_fasta += fasta_count
            other_embedded += embedded_count

    if other_fasta > 0:
        print(f"{'Other species':<15} | {other_fasta:<5} | {other_embedded:<8} | {other_fasta-other_embedded:<7} | {(other_embedded/other_fasta)*100:<6.1f}%")
        total_fasta_seqs += other_fasta
        total_embedded_seqs += other_embedded

    overall_coverage = (total_embedded_seqs / total_fasta_seqs) * 100 if total_fasta_seqs > 0 else 0

    print("-" * 65)
    print(f"{'TOTAL':<15} | {total_fasta_seqs:<5} | {total_embedded_seqs:<8} | {total_fasta_seqs-total_embedded_seqs:<7} | {overall_coverage:<6.1f}%")

    # Final assessment
    print(f"\n🎯 FINAL ASSESSMENT:")
    print(f"Overall coverage: {overall_coverage:.1f}%")
    print(f"Zero vectors: {zero_percentage:.1f}%")

    # Show missing sequences if any
    if missing_by_species:
        print(f"\n❌ Missing embeddings by species:")
        for species, missing_list in missing_by_species.items():
            print(f"  {species}: {missing_list[:3]}{'...' if len(missing_list) > 3 else ''}")

    # Pass/Fail determination
    overall_pass = (overall_coverage > 95.0 and
                   zero_percentage == 0.0 and
                   all_species_pass)

    if overall_pass:
        print(f"\n✅ PASS - Validation successful!")
        print(f"Ready for figure generation")
    else:
        print(f"\n❌ FAIL - Validation issues detected:")
        if overall_coverage <= 95.0:
            print(f"  - Overall coverage too low: {overall_coverage:.1f}% (need >95%)")
        if zero_percentage > 0:
            print(f"  - Zero vectors detected: {zero_percentage:.1f}%")
        if not all_species_pass:
            print(f"  - Some species below 95% coverage threshold")

    return 0 if overall_pass else 1

if __name__ == "__main__":
    exit(main())