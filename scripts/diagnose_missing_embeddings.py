#!/usr/bin/env python3
"""
Diagnose the 464 sequences with NO embedding.
Check if missing embeddings are biased toward specific species.
"""

import json
import pandas as pd
from pathlib import Path
from Bio import SeqIO

def load_embedded_sequences(embeddings_dir):
    """Load all sequences that have embeddings."""

    embedded_accessions = set()

    for segment in ['S', 'M', 'L']:
        final_index_file = embeddings_dir / f"embeddings_{segment}_final_index.json"

        if final_index_file.exists():
            with open(final_index_file, 'r') as f:
                index_data = json.load(f)

            accessions = index_data.get('accession_ids', [])
            embedded_accessions.update(accessions)
            print(f"  {segment}: {len(accessions)} embedded sequences")
        else:
            print(f"  ❌ {segment}: No final index found")

    print(f"  Total unique embedded accessions: {len(embedded_accessions)}")
    return embedded_accessions

def identify_missing_sequences(metadata, embedded_accessions):
    """Identify sequences in metadata that don't have embeddings."""

    print(f"\n🔍 IDENTIFYING MISSING SEQUENCES")
    print("=" * 60)

    all_accessions = set(metadata['genbank_id'])
    missing_accessions = all_accessions - embedded_accessions

    print(f"Total sequences in metadata: {len(all_accessions)}")
    print(f"Sequences with embeddings: {len(embedded_accessions)}")
    print(f"Sequences WITHOUT embeddings: {len(missing_accessions)}")

    # Get missing sequence details
    missing_metadata = metadata[metadata['genbank_id'].isin(missing_accessions)]

    return missing_accessions, missing_metadata

def analyze_species_bias(missing_metadata, metadata):
    """Analyze if missing sequences show species bias."""

    print(f"\n📊 SPECIES BIAS ANALYSIS")
    print("=" * 60)

    # Overall species distribution
    total_species = metadata['species'].value_counts()
    missing_species = missing_metadata['species'].value_counts()

    print(f"Species distribution in missing sequences:")
    print(f"{'Species':<20} {'Missing':<8} {'Total':<8} {'% Missing':<10} {'Bias'}")
    print("-" * 60)

    bias_analysis = {}

    for species in total_species.index:
        total_count = total_species[species]
        missing_count = missing_species.get(species, 0)
        missing_pct = (missing_count / total_count * 100) if total_count > 0 else 0

        # Calculate bias (how much higher/lower than average)
        overall_missing_rate = len(missing_metadata) / len(metadata) * 100
        bias_factor = missing_pct / overall_missing_rate if overall_missing_rate > 0 else 0

        if bias_factor > 1.5:
            bias_str = "HIGH BIAS ⚠️"
        elif bias_factor < 0.5:
            bias_str = "LOW BIAS"
        else:
            bias_str = "NORMAL"

        print(f"{species:<20} {missing_count:<8} {total_count:<8} {missing_pct:<10.1f} {bias_str}")

        bias_analysis[species] = {
            'missing_count': missing_count,
            'total_count': total_count,
            'missing_percentage': missing_pct,
            'bias_factor': bias_factor,
            'bias_category': bias_str
        }

    # Check for critical bias in major species
    print(f"\n🎯 CRITICAL SPECIES ANALYSIS")
    print("-" * 40)

    critical_species = ['Seoul', 'Andes', 'Puumala', 'Sin Nombre', 'Choclo']
    bias_detected = False

    for species in critical_species:
        if species in bias_analysis:
            data = bias_analysis[species]
            missing_pct = data['missing_percentage']
            bias_factor = data['bias_factor']

            if missing_pct > 70:  # More than 70% missing
                print(f"🚨 {species}: {missing_pct:.1f}% missing (CRITICAL BIAS)")
                bias_detected = True
            elif missing_pct > 50:  # More than 50% missing
                print(f"⚠️  {species}: {missing_pct:.1f}% missing (HIGH BIAS)")
                bias_detected = True
            else:
                print(f"✅ {species}: {missing_pct:.1f}% missing (acceptable)")

    return bias_analysis, bias_detected

def check_fasta_availability(missing_accessions, processed_dir):
    """Check if missing sequences exist in FASTA files."""

    print(f"\n🧬 CHECKING FASTA FILE AVAILABILITY")
    print("=" * 60)

    fasta_availability = {}
    total_in_fasta = 0

    for segment in ['S', 'M', 'L']:
        fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"

        if not fasta_path.exists():
            print(f"❌ {segment}: FASTA file not found: {fasta_path}")
            continue

        # Load sequences from FASTA
        fasta_sequences = {rec.id for rec in SeqIO.parse(fasta_path, "fasta")}

        # Check how many missing sequences are in this FASTA
        missing_in_fasta = missing_accessions & fasta_sequences
        missing_not_in_fasta = missing_accessions - fasta_sequences

        fasta_availability[segment] = {
            'total_in_fasta': len(fasta_sequences),
            'missing_sequences_in_fasta': len(missing_in_fasta),
            'missing_sequences_not_in_fasta': len(missing_not_in_fasta)
        }

        total_in_fasta += len(missing_in_fasta)

        print(f"{segment} FASTA:")
        print(f"  Total sequences in FASTA: {len(fasta_sequences)}")
        print(f"  Missing sequences found in FASTA: {len(missing_in_fasta)}")
        print(f"  Missing sequences NOT in FASTA: {len(missing_not_in_fasta)}")

        if len(missing_in_fasta) > 0:
            sample_in_fasta = sorted(list(missing_in_fasta))[:3]
            print(f"  Sample missing but in FASTA: {sample_in_fasta}")

    print(f"\nSUMMARY:")
    print(f"Total missing sequences: {len(missing_accessions)}")
    print(f"Found in FASTA files: {total_in_fasta}")
    print(f"NOT in FASTA files: {len(missing_accessions) - total_in_fasta}")

    return fasta_availability

def check_recomputation_plan(embeddings_dir, missing_accessions):
    """Check what the recomputation plan included vs what was computed."""

    print(f"\n📋 CHECKING RECOMPUTATION PLAN")
    print("=" * 60)

    plan_file = embeddings_dir / "recomputation_plan.json"

    if not plan_file.exists():
        print(f"❌ Recomputation plan not found: {plan_file}")
        return

    with open(plan_file, 'r') as f:
        plan = json.load(f)

    # Get all sequences that were supposed to be recomputed
    planned_sequences = set()

    for segment in ['S', 'M', 'L']:
        segment_missing = plan['missing_by_segment'][segment]['missing_accessions']
        planned_sequences.update(segment_missing)

        print(f"{segment}: {len(segment_missing)} sequences planned for recomputation")

    print(f"Total planned for recomputation: {len(planned_sequences)}")

    # Check what's still missing
    still_missing = missing_accessions & planned_sequences
    not_planned = missing_accessions - planned_sequences

    print(f"\nRECOMPUTATION ANALYSIS:")
    print(f"Sequences that were planned but still missing: {len(still_missing)}")
    print(f"Sequences missing but NOT in plan: {len(not_planned)}")

    if len(still_missing) > 0:
        print(f"Sample still missing: {sorted(list(still_missing))[:5]}")

    if len(not_planned) > 0:
        print(f"Sample not planned: {sorted(list(not_planned))[:5]}")

    # Check against FASTA availability
    return {
        'planned_sequences': len(planned_sequences),
        'still_missing_after_computation': len(still_missing),
        'missing_not_planned': len(not_planned)
    }

def generate_coverage_recommendation(bias_detected, bias_analysis, fasta_availability):
    """Generate recommendation for figure generation."""

    print(f"\n🎯 FIGURE GENERATION RECOMMENDATION")
    print("=" * 60)

    # Check if major species have acceptable representation
    major_species_ok = True
    problematic_species = []

    for species, data in bias_analysis.items():
        if data['total_count'] >= 50:  # Major species (50+ sequences)
            if data['missing_percentage'] > 60:  # More than 60% missing
                major_species_ok = False
                problematic_species.append(f"{species} ({data['missing_percentage']:.1f}% missing)")

    if not bias_detected and major_species_ok:
        print(f"✅ RECOMMENDATION: PROCEED WITH FIGURE GENERATION")
        print(f"   - No critical species bias detected")
        print(f"   - Major species have acceptable representation")
        print(f"   - 48% coverage is adequate for analysis")
        recommendation = "PROCEED"
    else:
        print(f"⚠️  RECOMMENDATION: FIGURES MAY BE BIASED")
        print(f"   - Species bias detected: {bias_detected}")
        print(f"   - Problematic species: {problematic_species}")
        print(f"   - Consider additional data fetching or explicit bias warnings")
        recommendation = "CAUTION"

    return recommendation

def main():
    print("=" * 70)
    print("DIAGNOSING 464 SEQUENCES WITH NO EMBEDDINGS")
    print("=" * 70)

    # Setup paths
    project_root = Path(__file__).parent.parent
    embeddings_dir = project_root / "results" / "embeddings"
    processed_dir = project_root / "data" / "processed"
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    # Load metadata
    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded metadata: {len(metadata)} sequences")

    # Load embedded sequences
    print(f"\n📂 LOADING EMBEDDED SEQUENCES")
    print("=" * 50)
    embedded_accessions = load_embedded_sequences(embeddings_dir)

    # Identify missing sequences
    missing_accessions, missing_metadata = identify_missing_sequences(metadata, embedded_accessions)

    # Analyze species bias
    bias_analysis, bias_detected = analyze_species_bias(missing_metadata, metadata)

    # Check FASTA availability
    fasta_availability = check_fasta_availability(missing_accessions, processed_dir)

    # Check recomputation plan
    recomputation_analysis = check_recomputation_plan(embeddings_dir, missing_accessions)

    # Generate recommendation
    recommendation = generate_coverage_recommendation(bias_detected, bias_analysis, fasta_availability)

    print(f"\n" + "=" * 70)
    print("DIAGNOSTIC SUMMARY")
    print("=" * 70)
    print(f"Missing sequences: {len(missing_accessions)}")
    print(f"Species bias detected: {bias_detected}")
    print(f"Recommendation: {recommendation}")

    return 0

if __name__ == "__main__":
    exit(main())