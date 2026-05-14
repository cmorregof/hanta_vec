#!/usr/bin/env python3
"""
Clean embedding recomputation from scratch using single visualization FASTA.
"""

import numpy as np
import json
import torch
import re
from Bio import SeqIO
from pathlib import Path
from transformers import EsmModel, EsmTokenizer
from tqdm import tqdm

def get_segment_from_extraction_method(header):
    """Extract segment using extraction method as ground truth."""

    SEGMENT_MAP = {
        # S-segment: nucleocapsid
        'n_protein_full': 'S',
        # M-segment: glycoprotein (all variants)
        'gpc_full': 'M',
        'gpc_truncated_1022': 'M',
        'gn_domain': 'M',
        'gpc_nterm': 'M',
        'nucleotide_tier1_gpc': 'M',
        'nucleotide_tier2_gpc': 'M',
        'nucleotide_tier2_large_cds': 'M',
        # L-segment: polymerase
        'rdrp_full': 'L',
        # Ambiguous — exclude entirely
        'protein_db_direct': None,
        'protein_db': None,
    }

    # First try accession suffix (_S, _M, _L)
    if re.search(r'_[SML]\s', header):
        return re.search(r'_([SML])\s', header).group(1)

    # Then use extraction method mapping (ground truth)
    for method, segment in SEGMENT_MAP.items():
        if method in header:
            return segment  # None for excluded methods

    return None

def extract_species_from_header(header):
    """Extract species name from header."""
    # Format: >ACCESSION[_SEGMENT] Species | isolate | method
    match = re.match(r'^[^\s]+\s+([^|]+)', header)
    if match:
        return match.group(1).strip()
    return 'unknown'

def load_esm_model():
    """Load ESM-2 model for embedding computation."""

    print("🤖 Loading ESM-2 model...")

    model_name = "facebook/esm2_t12_35M_UR50D"
    tokenizer = EsmTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)

    # Use MPS if available
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    model = model.to(device)
    model.eval()

    print(f"✓ Model loaded: {model_name}")
    print(f"✓ Device: {device}")
    print(f"✓ Max sequence length: 1022")

    return model, tokenizer, device

def compute_embedding(sequence, model, tokenizer, device, max_length=1022):
    """Compute ESM-2 embedding for a single sequence."""

    # Truncate if necessary
    if len(sequence) > max_length:
        sequence = sequence[:max_length]

    # Tokenize
    inputs = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=max_length+2)
    inputs = {k: v.to(device) for k, v in inputs.items()}

    # Compute embedding
    with torch.no_grad():
        outputs = model(**inputs)
        # Use mean pooling over sequence length (excluding special tokens)
        embedding = outputs.last_hidden_state[0, 1:-1, :].mean(dim=0)

    return embedding.cpu().numpy()

def process_sequences_by_segment():
    """Process all sequences from visualization FASTA, organized by segment."""

    print("📋 Processing sequences from visualization FASTA...")

    visualization_file = Path("data/processed/sequences_level1_for_visualization.fasta")

    sequences_by_segment = {'S': [], 'M': [], 'L': []}
    excluded_count = 0

    for record in SeqIO.parse(visualization_file, "fasta"):
        # Extract segment
        segment = get_segment_from_extraction_method(record.description)

        if segment is None:
            excluded_count += 1
            continue

        # Extract base accession
        base_accession = record.id.split('_')[0] if '_' in record.id else record.id

        # Extract species
        species = extract_species_from_header(record.description)

        # Store sequence info
        seq_info = {
            'accession': base_accession,
            'record_id': record.id,
            'species': species,
            'sequence': str(record.seq),
            'description': record.description
        }

        sequences_by_segment[segment].append(seq_info)

    print(f"✓ Processed {visualization_file}")
    for segment in ['S', 'M', 'L']:
        print(f"  {segment}: {len(sequences_by_segment[segment])} sequences")
    print(f"  Excluded (ambiguous methods): {excluded_count}")

    return sequences_by_segment

def compute_embeddings_for_segment(sequences, segment, model, tokenizer, device):
    """Compute embeddings for all sequences in a segment."""

    print(f"\n🧬 Computing {segment} segment embeddings...")

    embeddings = []
    accessions = []
    metadata = []
    failed_count = 0
    truncated_count = 0

    progress_bar = tqdm(sequences, desc=f"{segment} embeddings", unit="seq")

    for seq_info in progress_bar:
        try:
            sequence = seq_info['sequence']
            accession = seq_info['accession']

            # Track truncations
            if len(sequence) > 1022:
                truncated_count += 1

            # Compute embedding
            embedding = compute_embedding(sequence, model, tokenizer, device)

            # Validate embedding (no zero vectors)
            embedding_norm = np.linalg.norm(embedding)
            if embedding_norm < 1e-6:
                print(f"    Warning: Near-zero embedding for {accession} (norm: {embedding_norm})")
                failed_count += 1
                continue

            # Store results
            embeddings.append(embedding)
            accessions.append(accession)
            metadata.append(seq_info)

        except Exception as e:
            print(f"    Error processing {seq_info['accession']}: {e}")
            failed_count += 1
            continue

    progress_bar.close()

    # Convert to numpy array
    if embeddings:
        embeddings_array = np.stack(embeddings)
    else:
        embeddings_array = np.empty((0, 480))  # Empty array with correct shape

    print(f"✓ {segment} Results:")
    print(f"  Successful embeddings: {len(embeddings)}")
    print(f"  Failed embeddings: {failed_count}")
    print(f"  Truncated sequences: {truncated_count}")
    print(f"  Final shape: {embeddings_array.shape}")
    print(f"  Mean norm: {np.mean([np.linalg.norm(emb) for emb in embeddings]):.6f}")

    return embeddings_array, accessions, metadata

def save_segment_embeddings(embeddings, accessions, metadata, segment):
    """Save embeddings and index for a segment."""

    embeddings_dir = Path("results/embeddings")

    # Save embeddings
    embedding_file = embeddings_dir / f"embeddings_{segment}.npy"
    np.save(embedding_file, embeddings)

    # Create index
    index = {
        'accession_ids': accessions,
        'embedding_shape': embeddings.shape,
        'metadata': metadata,
        'segment': segment
    }

    # Save index
    index_file = embeddings_dir / f"embeddings_{segment}_index.json"
    with open(index_file, 'w') as f:
        json.dump(index, f, indent=2, default=str)

    print(f"💾 Saved {segment} embeddings:")
    print(f"  Embeddings: {embedding_file}")
    print(f"  Index: {index_file}")

    return len(accessions)

def main():
    print("=" * 70)
    print("CLEAN EMBEDDING RECOMPUTATION FROM SCRATCH")
    print("=" * 70)

    # Process sequences from single FASTA
    sequences_by_segment = process_sequences_by_segment()

    # Load ESM-2 model
    model, tokenizer, device = load_esm_model()

    total_computed = 0
    total_failed = 0

    # Process each segment
    for segment in ['S', 'M', 'L']:
        if sequences_by_segment[segment]:
            embeddings, accessions, metadata = compute_embeddings_for_segment(
                sequences_by_segment[segment], segment, model, tokenizer, device
            )

            # Save results
            segment_count = save_segment_embeddings(embeddings, accessions, metadata, segment)
            total_computed += segment_count
        else:
            print(f"\n🧬 No sequences found for {segment} segment")

    # Final summary
    print(f"\n{'='*50}")
    print("CLEAN RECOMPUTATION SUMMARY")
    print(f"{'='*50}")

    total_input = sum(len(sequences_by_segment[seg]) for seg in ['S', 'M', 'L'])
    success_rate = (total_computed / total_input) * 100 if total_input > 0 else 0

    print(f"Total input sequences: {total_input}")
    print(f"Total embeddings computed: {total_computed}")
    print(f"Success rate: {success_rate:.1f}%")

    print(f"\n✅ Clean recomputation complete!")
    print(f"Ready for validation step.")

    return 0

if __name__ == "__main__":
    exit(main())