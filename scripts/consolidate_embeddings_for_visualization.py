#!/usr/bin/env python3
"""
Consolidate three-segment embeddings for visualization.
Creates combined embedding matrix and accession list expected by figure script.
"""

import numpy as np
import pandas as pd
from pathlib import Path

def main():
    print("="*60)
    print("CONSOLIDATING THREE-SEGMENT EMBEDDINGS FOR VISUALIZATION")
    print("="*60)

    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load metadata to get the correct order
    metadata = pd.read_csv(processed_dir / "metadata_level1.tsv", sep="\t")
    print(f"Total metadata entries: {len(metadata)}")

    # Load each segment's embeddings and accessions
    all_embeddings = []
    all_accessions = []

    for segment in ['S', 'M', 'L']:
        print(f"\nLoading {segment}-segment:")

        # Load embeddings
        emb_path = embeddings_dir / f"embeddings_{segment}.npy"
        embeddings = np.load(emb_path)
        print(f"  Embeddings shape: {embeddings.shape}")

        # Load accessions
        acc_path = embeddings_dir / f"accessions_{segment}.txt"
        with open(acc_path) as f:
            accessions = [line.strip() for line in f]
        print(f"  Accessions: {len(accessions)}")

        # Verify consistency
        if len(embeddings) != len(accessions):
            print(f"  ❌ Mismatch: {len(embeddings)} embeddings vs {len(accessions)} accessions")
            return 1

        all_embeddings.append(embeddings)
        all_accessions.extend(accessions)

    # Combine embeddings vertically
    combined_embeddings = np.vstack(all_embeddings)
    print(f"\nCombined embeddings shape: {combined_embeddings.shape}")
    print(f"Total accessions: {len(all_accessions)}")

    # Verify against metadata
    metadata_accessions = set(metadata['genbank_id'].values)
    embedding_accessions = set(all_accessions)

    missing_in_embeddings = metadata_accessions - embedding_accessions
    extra_in_embeddings = embedding_accessions - metadata_accessions

    print(f"\nConsistency check:")
    print(f"  Metadata accessions: {len(metadata_accessions)}")
    print(f"  Embedding accessions: {len(embedding_accessions)}")
    print(f"  Missing in embeddings: {len(missing_in_embeddings)}")
    print(f"  Extra in embeddings: {len(extra_in_embeddings)}")

    if missing_in_embeddings:
        print(f"  Missing examples: {list(missing_in_embeddings)[:5]}")
    if extra_in_embeddings:
        print(f"  Extra examples: {list(extra_in_embeddings)[:5]}")

    # Save consolidated files for visualization
    print(f"\nSaving consolidated files:")

    # Save embeddings
    output_emb_path = embeddings_dir / "embeddings_level1.npy"
    np.save(output_emb_path, combined_embeddings)
    print(f"  ✓ Saved embeddings to {output_emb_path}")

    # Save accessions
    output_acc_path = embeddings_dir / "accessions_level1.txt"
    with open(output_acc_path, "w") as f:
        f.write("\n".join(all_accessions))
    print(f"  ✓ Saved accessions to {output_acc_path}")

    # Filter metadata to match available embeddings
    metadata_filtered = metadata[metadata['genbank_id'].isin(embedding_accessions)].copy()

    # Sort metadata to match accession order
    accession_to_idx = {acc: i for i, acc in enumerate(all_accessions)}
    metadata_filtered['embedding_order'] = metadata_filtered['genbank_id'].map(accession_to_idx)
    metadata_filtered = metadata_filtered.sort_values('embedding_order').drop('embedding_order', axis=1)

    print(f"\nFiltered metadata shape: {metadata_filtered.shape}")
    print(f"Matches embedding count: {len(metadata_filtered) == len(combined_embeddings)}")

    # Save metadata subset that matches embeddings
    output_meta_path = processed_dir / "metadata_level1_with_embeddings.tsv"
    metadata_filtered.to_csv(output_meta_path, sep="\t", index=False)
    print(f"  ✓ Saved filtered metadata to {output_meta_path}")

    # Create a consolidated FASTA for visualization (use first available sequence per accession)
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    output_fasta_path = processed_dir / "sequences_level1_for_visualization.fasta"
    with open(output_fasta_path, "w") as f:
        for segment in ['S', 'M', 'L']:
            fasta_path = processed_dir / f"{segment}_sequences_level1.fasta"
            for record in SeqIO.parse(fasta_path, "fasta"):
                if record.id in embedding_accessions:
                    SeqIO.write(record, f, "fasta")

    print(f"  ✓ Created visualization FASTA: {output_fasta_path}")

    print(f"\n✅ Consolidation complete!")
    print(f"Ready for visualization with {len(combined_embeddings)} sequences")

    return 0

if __name__ == "__main__":
    exit(main())