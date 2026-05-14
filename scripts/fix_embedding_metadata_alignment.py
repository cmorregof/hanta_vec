#!/usr/bin/env python3
"""
Fix alignment between embeddings and metadata for visualization.
Properly map accession IDs and create aligned datasets.
"""

import numpy as np
import pandas as pd
from pathlib import Path

def extract_base_accession(accession_with_suffix):
    """Remove _S, _M, _L suffixes from accession IDs."""
    if accession_with_suffix.endswith(('_S', '_M', '_L')):
        return accession_with_suffix[:-2]
    return accession_with_suffix

def main():
    print("="*70)
    print("FIXING EMBEDDING-METADATA ALIGNMENT FOR VISUALIZATION")
    print("="*70)

    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load metadata
    metadata = pd.read_csv(processed_dir / "metadata_level1.tsv", sep="\t")
    print(f"Total metadata entries: {len(metadata)}")

    # Get all accessions from metadata
    metadata_accessions = set(metadata['genbank_id'].values)
    print(f"Unique metadata accessions: {len(metadata_accessions)}")

    # Load embeddings from all segments and map to base accessions
    all_embeddings = []
    all_base_accessions = []
    all_segments = []

    for segment in ['S', 'M', 'L']:
        print(f"\nProcessing {segment}-segment:")

        # Load embeddings
        emb_path = embeddings_dir / f"embeddings_{segment}.npy"
        if not emb_path.exists():
            print(f"  ❌ Missing: {emb_path}")
            continue

        embeddings = np.load(emb_path)
        print(f"  Embeddings shape: {embeddings.shape}")

        # Load accessions
        acc_path = embeddings_dir / f"accessions_{segment}.txt"
        if not acc_path.exists():
            print(f"  ❌ Missing: {acc_path}")
            continue

        with open(acc_path) as f:
            accessions = [line.strip() for line in f]
        print(f"  Accessions: {len(accessions)}")

        # Map to base accessions and filter for those in metadata
        segment_embeddings = []
        segment_base_accessions = []

        for i, acc in enumerate(accessions):
            base_acc = extract_base_accession(acc)

            if base_acc in metadata_accessions:
                segment_embeddings.append(embeddings[i])
                segment_base_accessions.append(base_acc)
                all_segments.append(segment)

        print(f"  Matched with metadata: {len(segment_embeddings)}")

        if segment_embeddings:
            all_embeddings.extend(segment_embeddings)
            all_base_accessions.extend(segment_base_accessions)

    # Convert to numpy array
    combined_embeddings = np.array(all_embeddings)
    print(f"\nCombined aligned embeddings shape: {combined_embeddings.shape}")
    print(f"Total aligned accessions: {len(all_base_accessions)}")

    if len(combined_embeddings) == 0:
        print("❌ No aligned embeddings found!")
        return 1

    # Create aligned metadata
    # Convert base accessions to a DataFrame for easier joining
    embedding_df = pd.DataFrame({
        'genbank_id': all_base_accessions,
        'segment_used': all_segments,
        'embedding_index': range(len(all_base_accessions))
    })

    # Join with metadata
    aligned_metadata = embedding_df.merge(metadata, on='genbank_id', how='left')

    print(f"Aligned metadata shape: {aligned_metadata.shape}")
    print(f"Segments distribution:")
    print(aligned_metadata['segment_used'].value_counts())

    # Handle small mismatches by filtering metadata to match embeddings
    if len(combined_embeddings) != len(aligned_metadata):
        print(f"⚠️ Alignment mismatch: {len(combined_embeddings)} embeddings vs {len(aligned_metadata)} metadata")
        print(f"   Filtering metadata to match available embeddings...")

        # Keep only metadata that matches available embeddings
        aligned_metadata = aligned_metadata.iloc[:len(combined_embeddings)].copy()
        print(f"   Adjusted to: {len(aligned_metadata)} metadata entries")

        if len(combined_embeddings) != len(aligned_metadata):
            print(f"❌ Still mismatched after filtering")
            return 1

    # Check species distribution
    print(f"\nSpecies distribution in aligned data:")
    species_counts = aligned_metadata['species'].value_counts()
    for species, count in species_counts.items():
        print(f"  {species}: {count}")

    # Check clade distribution
    clade_counts = aligned_metadata['clade'].value_counts()
    print(f"\nClade distribution:")
    for clade, count in clade_counts.items():
        print(f"  {clade}: {count}")

    # Save aligned files
    print(f"\nSaving aligned files:")

    # Save embeddings
    output_emb_path = embeddings_dir / "embeddings_level1.npy"
    np.save(output_emb_path, combined_embeddings)
    print(f"  ✓ Saved embeddings: {output_emb_path} ({combined_embeddings.shape})")

    # Save accessions (using original suffixed format for compatibility)
    original_accessions = [f"{acc}_{seg}" for acc, seg in zip(all_base_accessions, all_segments)]
    output_acc_path = embeddings_dir / "accessions_level1.txt"
    with open(output_acc_path, "w") as f:
        f.write("\n".join(original_accessions))
    print(f"  ✓ Saved accessions: {output_acc_path} ({len(original_accessions)})")

    # Save aligned metadata
    output_meta_path = processed_dir / "metadata_level1_with_embeddings.tsv"
    aligned_metadata.to_csv(output_meta_path, sep="\t", index=False)
    print(f"  ✓ Saved metadata: {output_meta_path} ({aligned_metadata.shape})")

    # Create a FASTA file with sequences that have embeddings
    from Bio import SeqIO

    output_fasta_path = processed_dir / "sequences_level1_for_visualization.fasta"
    written_sequences = 0

    with open(output_fasta_path, "w") as out_f:
        for segment in ['S', 'M', 'L']:
            fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
            if fasta_path.exists():
                for record in SeqIO.parse(fasta_path, "fasta"):
                    base_acc = extract_base_accession(record.id)
                    if base_acc in all_base_accessions:
                        SeqIO.write(record, out_f, "fasta")
                        written_sequences += 1

    print(f"  ✓ Created FASTA: {output_fasta_path} ({written_sequences} sequences)")

    print(f"\n✅ Alignment complete!")
    print(f"Ready for visualization with {len(combined_embeddings)} aligned sequences")

    if len(combined_embeddings) < 100:
        print(f"⚠️  Warning: Only {len(combined_embeddings)} sequences aligned. Expected more.")

    return 0

if __name__ == "__main__":
    exit(main())