#!/usr/bin/env python3
"""
Diagnostic: Analyze embedding coverage gap (787 sequences vs 429 embeddings).
"""

import numpy as np
import json
from Bio import SeqIO
from pathlib import Path
import re

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

def analyze_fasta_sequences():
    """Analyze sequences by species and segment from FASTA."""
    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    species_data = {}

    for record in SeqIO.parse(visualization_file, "fasta"):
        # Extract species
        description = record.description
        match = re.match(r'^[^\s]+\s+([^|]+)', description)
        if match:
            species = match.group(1).strip()
        else:
            species = 'unknown'

        # Extract segment
        segment = get_segment_from_extraction_method(description)
        if segment is None:
            continue

        # Get base accession
        base_acc = record.id.split('_')[0] if '_' in record.id else record.id

        if species not in species_data:
            species_data[species] = {'S': set(), 'M': set(), 'L': set()}

        species_data[species][segment].add(base_acc)

    return species_data

def analyze_embedding_coverage():
    """Analyze which sequences have embeddings."""
    embeddings_dir = Path("results/embeddings")

    embedded_accessions = {'S': set(), 'M': set(), 'L': set()}
    embedding_counts = {'S': 0, 'M': 0, 'L': 0}

    for segment in ['S', 'M', 'L']:
        # Check final embedding files
        embedding_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        index_file = embeddings_dir / f"embeddings_{segment}_final_index.json"

        if embedding_file.exists() and index_file.exists():
            embeddings = np.load(embedding_file)
            with open(index_file, 'r') as f:
                index = json.load(f)

            # Get actual accession list (skip metadata keys)
            accessions = [k for k in index.keys() if not k.startswith(('accession_ids', 'embedding_shape', 'clean_count'))]

            embedded_accessions[segment].update(accessions)
            embedding_counts[segment] = len(accessions)

            print(f"  {segment} segment: {len(accessions)} embeddings in final file")

    return embedded_accessions, embedding_counts

def check_file_sizes():
    """Check current embedding file sizes."""
    embeddings_dir = Path("results/embeddings")

    print("\n📊 EMBEDDING FILE SIZES:")
    print(f"{'File':<30} {'Size':<15} {'Shape':<15}")
    print("-" * 60)

    for segment in ['S', 'M', 'L']:
        # Final files
        final_file = embeddings_dir / f"embeddings_{segment}_final.npy"
        if final_file.exists():
            embeddings = np.load(final_file)
            file_size = f"{final_file.stat().st_size / (1024*1024):.1f} MB"
            print(f"embeddings_{segment}_final.npy{'':<7} {file_size:<15} {embeddings.shape}")

        # New files (from recent computation)
        new_file = embeddings_dir / f"embeddings_{segment}_new.npy"
        if new_file.exists():
            embeddings = np.load(new_file)
            file_size = f"{new_file.stat().st_size / (1024*1024):.1f} MB"
            print(f"embeddings_{segment}_new.npy{'':<8} {file_size:<15} {embeddings.shape}")

def main():
    print("=" * 70)
    print("DIAGNOSTIC: EMBEDDING COVERAGE GAP ANALYSIS")
    print("=" * 70)

    # Analyze sequences in FASTA
    print("\n🔍 ANALYZING FASTA SEQUENCES...")
    fasta_species = analyze_fasta_sequences()

    # Analyze embeddings
    print("\n🔍 ANALYZING EMBEDDINGS...")
    embedded_accessions, embedding_counts = analyze_embedding_coverage()

    # Check file sizes
    check_file_sizes()

    # Generate coverage report
    print(f"\n📋 SPECIES EMBEDDING COVERAGE REPORT:")
    print(f"{'Species':<15} | {'FASTA':<5} | {'Embedded':<8} | {'Missing':<7} | {'Coverage':<8}")
    print("-" * 65)

    target_species = ['Seoul', 'Andes', 'Puumala', 'Prospect Hill', 'Choclo', 'Sin Nombre']

    total_fasta = 0
    total_embedded = 0

    for species in target_species:
        if species in fasta_species:
            # Count sequences in FASTA
            fasta_accessions = set()
            for segment in ['S', 'M', 'L']:
                fasta_accessions.update(fasta_species[species][segment])
            fasta_count = len(fasta_accessions)

            # Count embeddings
            embedded_count = 0
            for segment in ['S', 'M', 'L']:
                # Find intersection with embedded accessions
                embedded_for_species = embedded_accessions[segment].intersection(fasta_species[species][segment])
                embedded_count += len(embedded_for_species)

            missing_count = fasta_count - embedded_count
            coverage_pct = (embedded_count / fasta_count) * 100 if fasta_count > 0 else 0

            status = "✓" if coverage_pct > 85 or fasta_count < 5 else "❌"

            print(f"{species:<15} | {fasta_count:<5} | {embedded_count:<8} | {missing_count:<7} | {coverage_pct:<6.1f}% {status}")

            total_fasta += fasta_count
            total_embedded += embedded_count
        else:
            print(f"{species:<15} | {'0':<5} | {'0':<8} | {'0':<7} | {'0.0%':<8}")

    # Overall coverage
    overall_coverage = (total_embedded / total_fasta) * 100 if total_fasta > 0 else 0

    print("-" * 65)
    print(f"{'TOTAL':<15} | {total_fasta:<5} | {total_embedded:<8} | {total_fasta-total_embedded:<7} | {overall_coverage:<6.1f}%")

    print(f"\n🎯 COVERAGE ANALYSIS:")
    print(f"Overall coverage: {overall_coverage:.1f}%")
    print(f"Coverage gap: {100-overall_coverage:.1f}%")

    if overall_coverage < 85:
        print(f"❌ INSUFFICIENT COVERAGE for figure generation")
        print(f"Need >85% coverage for reliable analysis")
    else:
        print(f"✅ SUFFICIENT COVERAGE for figure generation")

    # Segment breakdown
    print(f"\n📊 SEGMENT BREAKDOWN:")
    for segment in ['S', 'M', 'L']:
        segment_fasta = sum(len(fasta_species.get(species, {}).get(segment, set()))
                          for species in target_species if species in fasta_species)
        segment_embedded = len(embedded_accessions[segment])
        segment_coverage = (segment_embedded / segment_fasta) * 100 if segment_fasta > 0 else 0
        print(f"{segment}: {segment_embedded}/{segment_fasta} ({segment_coverage:.1f}%)")

    return 0

if __name__ == "__main__":
    exit(main())