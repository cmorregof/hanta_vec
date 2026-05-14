#!/usr/bin/env python3
"""
Fix F8 with corrected embedding lookup and normalization.
Address zero-norm embeddings and index mismatch issues.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def find_matched_isolates_fixed(metadata):
    """Find isolates with sequences in multiple segments - FIXED VERSION."""

    print("  Finding matched isolates...")

    # Group by isolate and species to find cross-segment matches
    isolate_segment_map = {}

    for _, row in metadata.iterrows():
        isolate = row.get('isolate', '')
        species = row['species']
        segment = row['segment_used']
        genbank_id = row['genbank_id']
        record_id = row.get('record_id', '')

        # Handle NaN isolates by using 'empty' as placeholder
        if pd.isna(isolate) or isolate == '':
            isolate = 'empty'

        # Create a key for matching
        key = f"{species}_{isolate}"

        if key not in isolate_segment_map:
            isolate_segment_map[key] = {}

        if segment not in isolate_segment_map[key]:
            isolate_segment_map[key][segment] = []

        isolate_segment_map[key][segment].append({
            'genbank_id': genbank_id,
            'embedding_index': row.name,  # pandas index
            'record_id': record_id,
            'row_data': row
        })

    # Find isolates with multiple segments
    matched_isolates = {}

    for key, segments in isolate_segment_map.items():
        if len(segments) >= 2:  # Has at least 2 segments
            # Extract species from key
            species = key.split('_')[0]

            if species == 'Unknown':
                continue  # Skip Unknown species for cleaner analysis

            if species not in matched_isolates:
                matched_isolates[species] = []

            matched_isolates[species].append({
                'isolate_key': key,
                'segments': segments
            })

    # Debug output
    print(f"  Found matched isolates across segments:")
    total_isolates = 0
    for species, isolates in matched_isolates.items():
        total_isolates += len(isolates)
        segment_combinations = {}
        for isolate in isolates:
            segs = sorted(isolate['segments'].keys())
            combo = '+'.join(segs)
            segment_combinations[combo] = segment_combinations.get(combo, 0) + 1

        print(f"    {species}: {len(isolates)} isolates")
        for combo, count in segment_combinations.items():
            print(f"      {combo}: {count} isolates")

    print(f"  Total matched isolates: {total_isolates}")
    return matched_isolates

def create_segment_index_mapping(metadata):
    """Create mapping from record_id to segment-specific index."""

    print("  Creating segment index mapping...")

    segment_mappings = {}

    for segment in ['S', 'M', 'L']:
        segment_data = metadata[metadata['segment_used'] == segment].copy()
        segment_data = segment_data.reset_index(drop=True)  # Create new 0-based index

        # Map record_id to new segment index
        mapping = {}
        for new_idx, (_, row) in enumerate(segment_data.iterrows()):
            record_id = row.get('record_id', '')
            genbank_id = row['genbank_id']

            # Try both record_id and genbank_id as keys
            if record_id:
                mapping[record_id] = new_idx
            mapping[genbank_id] = new_idx

        segment_mappings[segment] = mapping
        print(f"    {segment}: {len(mapping)} mappings")

    return segment_mappings

def make_f8_fixed_embeddings(metadata, embeddings_dir, output_dir):
    """Generate F8 with proper embedding normalization and index mapping."""

    print("  Computing cross-segment embedding similarities with fixed lookup...")

    # Load embeddings
    embeddings_files = {
        'S': embeddings_dir / "embeddings_S.npy",
        'M': embeddings_dir / "embeddings_M.npy",
        'L': embeddings_dir / "embeddings_L.npy"
    }

    single_embedding_file = embeddings_dir / "embeddings_level1.npy"

    embeddings = {}

    # Load single embedding file (concatenated)
    if single_embedding_file.exists():
        print(f"    Loading single embedding file: {single_embedding_file}")
        single_embeddings = np.load(single_embedding_file)
        print(f"    Single embeddings shape: {single_embeddings.shape}")

        # Check if embeddings are zero vectors
        norms = np.linalg.norm(single_embeddings, axis=1)
        zero_count = np.sum(norms == 0)
        print(f"    Zero-norm embeddings: {zero_count}/{len(single_embeddings)}")

        if zero_count > len(single_embeddings) * 0.5:
            print(f"    ⚠️ Too many zero embeddings, trying segment-specific files...")
        else:
            embeddings['single'] = single_embeddings

    # Load segment-specific embeddings
    for seg, path in embeddings_files.items():
        if path.exists():
            seg_embeddings = np.load(path)
            print(f"    Loaded {seg}: shape {seg_embeddings.shape}")

            # Check norms
            norms = np.linalg.norm(seg_embeddings, axis=1)
            zero_count = np.sum(norms == 0)
            print(f"    {seg} zero-norm embeddings: {zero_count}/{len(seg_embeddings)}")

            embeddings[seg] = seg_embeddings

    if not embeddings:
        print("    ❌ No valid embeddings found!")
        return

    # Create segment index mapping
    segment_mappings = create_segment_index_mapping(metadata)

    # Find matched isolates
    matched_isolates = find_matched_isolates_fixed(metadata)

    if not matched_isolates:
        print("  ❌ No matched isolates found across segments")
        return

    # Compute pairwise similarities with fixed lookup
    similarity_data = []

    for species, isolates in matched_isolates.items():
        print(f"  Processing {species}: {len(isolates)} isolates...")

        for isolate_info in isolates:
            segments = isolate_info['segments']
            segment_pairs = [('S', 'M'), ('S', 'L'), ('M', 'L')]

            for seg1, seg2 in segment_pairs:
                if seg1 in segments and seg2 in segments:
                    try:
                        # Get first sequence from each segment for this isolate
                        seq1_info = segments[seg1][0]
                        seq2_info = segments[seg2][0]

                        # Get embeddings using segment-specific approach
                        emb1 = None
                        emb2 = None

                        # Try segment-specific embeddings first
                        if seg1 in embeddings and seg2 in embeddings:
                            # Map record_id/genbank_id to segment index
                            record1 = seq1_info.get('record_id', seq1_info['genbank_id'])
                            record2 = seq2_info.get('record_id', seq2_info['genbank_id'])

                            if record1 in segment_mappings[seg1] and record2 in segment_mappings[seg2]:
                                idx1 = segment_mappings[seg1][record1]
                                idx2 = segment_mappings[seg2][record2]

                                if idx1 < len(embeddings[seg1]) and idx2 < len(embeddings[seg2]):
                                    emb1 = embeddings[seg1][idx1]
                                    emb2 = embeddings[seg2][idx2]

                        # Fallback to single embeddings if available
                        if (emb1 is None or emb2 is None) and 'single' in embeddings:
                            emb1_idx = seq1_info['embedding_index']
                            emb2_idx = seq2_info['embedding_index']

                            if emb1_idx < len(embeddings['single']) and emb2_idx < len(embeddings['single']):
                                emb1 = embeddings['single'][emb1_idx]
                                emb2 = embeddings['single'][emb2_idx]

                        if emb1 is not None and emb2 is not None:
                            # Check if embeddings are valid (non-zero)
                            norm1 = np.linalg.norm(emb1)
                            norm2 = np.linalg.norm(emb2)

                            if norm1 > 1e-10 and norm2 > 1e-10:  # Avoid zero vectors
                                # Compute cosine similarity
                                cos_sim = np.dot(emb1, emb2) / (norm1 * norm2)

                                if np.isfinite(cos_sim):  # Check for valid similarity
                                    similarity_data.append({
                                        'species': species,
                                        'segment_pair': f'{seg1}-{seg2}',
                                        'cosine_similarity': cos_sim,
                                        'isolate': isolate_info['isolate_key']
                                    })

                                    # Debug successful computation
                                    if len(similarity_data) <= 10:  # First few
                                        print(f"      ✓ {species} {seg1}-{seg2}: {cos_sim:.6f}")
                            else:
                                print(f"      ❌ Zero norm: {species} {seg1}-{seg2} (norms: {norm1:.3e}, {norm2:.3e})")

                    except Exception as e:
                        print(f"    Error computing similarity for {species} {seg1}-{seg2}: {e}")
                        continue

    if not similarity_data:
        print("  ❌ No similarity data computed")
        return

    df_sim = pd.DataFrame(similarity_data)

    print(f"  Computed {len(similarity_data)} pairwise similarities")
    print(f"  Species breakdown:")
    for species, group in df_sim.groupby('species'):
        print(f"    {species}: {len(group)} pairs")

    # Create figure with proper species handling
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    segment_pairs = ['S-M', 'S-L', 'M-L']
    colors = plt.cm.Set2(np.linspace(0, 1, len(df_sim['species'].unique())))
    species_colors = dict(zip(df_sim['species'].unique(), colors))

    for i, pair in enumerate(segment_pairs):
        ax = axes[i]

        pair_data = df_sim[df_sim['segment_pair'] == pair].copy()

        if len(pair_data) > 0:
            print(f"  Creating plot for {pair}: {len(pair_data)} data points")

            # Get species with data for this segment pair
            species_with_data = pair_data['species'].unique()

            if len(species_with_data) > 0:
                try:
                    # Create violin plot
                    sns.violinplot(data=pair_data, x='species', y='cosine_similarity',
                                 ax=ax, palette='Set2')

                    # Add mean points and labels
                    for j, species in enumerate(species_with_data):
                        species_sims = pair_data[pair_data['species'] == species]['cosine_similarity']
                        if len(species_sims) > 0:
                            mean_sim = species_sims.mean()
                            if np.isfinite(mean_sim):
                                ax.scatter(j, mean_sim, color='red', s=50, zorder=10)
                                ax.text(j, mean_sim + 0.01, f'{mean_sim:.3f}',
                                       ha='center', va='bottom', fontsize=9, fontweight='bold')

                    ax.set_title(f'Segment {pair} Similarity\n(n={len(pair_data)} pairs)')
                    ax.set_ylabel('Cosine Similarity')
                    ax.set_ylim(0.90, 1.0)  # Focus on high similarity range
                    ax.grid(True, alpha=0.3)
                    ax.tick_params(axis='x', rotation=45)

                except Exception as e:
                    print(f"    Error creating violin plot for {pair}: {e}")
                    # Fallback to scatter plot
                    for species in species_with_data:
                        species_data = pair_data[pair_data['species'] == species]
                        ax.scatter([species] * len(species_data), species_data['cosine_similarity'],
                                 alpha=0.7, label=species)

                    ax.set_title(f'Segment {pair} Similarity\n(n={len(pair_data)} pairs)')
                    ax.set_ylabel('Cosine Similarity')
                    ax.set_ylim(0.90, 1.0)
                    ax.grid(True, alpha=0.3)

            else:
                ax.text(0.5, 0.5, f'No species data for\n{pair}', ha='center', va='center',
                       transform=ax.transAxes, fontsize=14)
                ax.set_title(f'Segment {pair} Similarity')

        else:
            ax.text(0.5, 0.5, f'No data for\n{pair}', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'Segment {pair} Similarity')

    plt.suptitle('Cross-Segment Embedding Similarity for Matched Isolates\n' +
                 'Higher similarity indicates consistent evolutionary signals across segments',
                 fontsize=16, y=1.02)
    plt.tight_layout()

    # Save figure
    output_path = output_dir / "small" / "F8_cross_segment_similarity.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved F8: {output_path}")

    # Print summary statistics
    print(f"\n  Cross-segment similarity summary:")
    for pair in segment_pairs:
        pair_data = df_sim[df_sim['segment_pair'] == pair]
        if len(pair_data) > 0:
            mean_sim = pair_data['cosine_similarity'].mean()
            std_sim = pair_data['cosine_similarity'].std()
            print(f"    {pair}: {mean_sim:.3f} ± {std_sim:.3f} (n={len(pair_data)})")

def main():
    print("="*70)
    print("FIXING F8 WITH CORRECTED EMBEDDING LOOKUP")
    print("="*70)

    # Load data
    embeddings_dir = Path(__file__).parent.parent / "results" / "embeddings"
    figures_dir = Path(__file__).parent.parent / "results" / "figures"
    processed_dir = Path(__file__).parent.parent / "data" / "processed"

    # Load aligned data
    metadata_path = processed_dir / "metadata_level1_with_embeddings.tsv"

    if not metadata_path.exists():
        print("❌ Required files not found!")
        return 1

    metadata = pd.read_csv(metadata_path, sep="\t")
    print(f"Loaded {len(metadata)} metadata entries")

    # Generate fixed F8
    print("\nGenerating Fixed F8: Cross-Segment Similarity Analysis...")
    make_f8_fixed_embeddings(metadata, embeddings_dir, figures_dir)

    print("\n✅ Fixed F8 generation complete with proper embedding lookup!")
    return 0

if __name__ == "__main__":
    exit(main())