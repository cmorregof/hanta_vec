"""
Create Level 0 (smoke test) and Level 1 (MVP) dataset splits.

Handles:
- Stratified sampling by species
- Capping per-species to avoid bias
- FASTA export
- Manifest generation (SHA256, accession, metadata)
"""

from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def stratified_sample_by_species(
    records_df: pd.DataFrame,
    records_dict: Dict[str, Dict],
    max_per_species: int = 50,
    seed: int = 42,
) -> Tuple[pd.DataFrame, Dict[str, Dict]]:
    """
    Sample stratified by species, with max per species.

    Args:
        records_df: metadata DataFrame
        records_dict: accession → {seq, metadata}
        max_per_species: max records per species
        seed: random seed

    Returns:
        (sampled_df, sampled_records_dict)
    """
    # Ensure species_clean exists
    if "species_clean" not in records_df.columns:
        records_df["species_clean"] = records_df.get("organism", "unknown").str.replace(" ", "_")

    sampled_df = records_df.groupby("species_clean", group_keys=False).apply(
        lambda x: x.sample(min(len(x), max_per_species), random_state=seed),
        include_groups=False
    )

    # Ensure the sampled_df has all required columns
    if "species_clean" not in sampled_df.columns:
        sampled_df["species_clean"] = sampled_df.get("organism", "unknown").str.replace(" ", "_")

    sampled_accessions = set(sampled_df["accession"].values)
    sampled_records = {
        acc: records_dict[acc]
        for acc in sampled_accessions
        if acc in records_dict
    }

    return sampled_df, sampled_records


def export_fasta(
    records_dict: Dict[str, Dict],
    output_path: Path,
):
    """Export sequences to FASTA file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seq_records = []
    for accession, data in records_dict.items():
        seq_record = SeqRecord(
            Seq(data["seq"]),
            id=accession,
            description=f"{data['metadata'].get('organism', '')}",
        )
        seq_records.append(seq_record)

    with open(output_path, "w") as f:
        SeqIO.write(seq_records, f, "fasta")

    print(f"✓ Exported {len(seq_records)} sequences to {output_path}")


def export_manifest(
    records_df: pd.DataFrame,
    records_dict: Dict[str, Dict],
    output_path: Path,
):
    """
    Export dataset manifest (SHA256, accession, species, metadata).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for _, row in records_df.iterrows():
        accession = row["accession"]
        if accession not in records_dict:
            continue

        data = records_dict[accession]
        manifest_rows.append({
            "accession": accession,
            "sha256_sequence": data.get("seq_hash", ""),
            "length": len(data["seq"]),
            "species": row["species_clean"],
            "country": row["country_norm"],
            "old_new_world": row["old_new_world"],
            "year": row["year"],
            "qc_status": "passed" if row.get("qc_passed") else "failed",
            "extraction_method": row.get("extraction_method", "unknown"),
        })

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(output_path, sep="\t", index=False)

    print(f"✓ Exported manifest with {len(manifest_df)} entries to {output_path}")


def build_level_0_and_1(
    records_dict: Dict[str, Dict],
    records_df: pd.DataFrame,
    output_dir: Path,
    level_0_count: Tuple[int, int] = (10, 20),  # (min, max)
    level_1_count: Tuple[int, int] = (100, 300),  # (min, max)
    max_per_species: int = 60,
) -> Tuple[Dict, Dict, pd.DataFrame, pd.DataFrame]:
    """
    Create Level 0 and Level 1 splits.

    Args:
        records_dict: accession → {seq, metadata}
        records_df: metadata DataFrame
        output_dir: where to save outputs
        level_0_count: (min, max) for Level 0
        level_1_count: (min, max) for Level 1
        max_per_species: cap per species before QC

    Returns:
        (level0_records, level1_records, level0_df, level1_df)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    total_records = len(records_dict)
    print(f"\nBuilding Level 0 and Level 1 from {total_records} QC-passed records")

    # Level 0: stratified, capped at max_per_species, target 10-20
    level0_df, level0_records = stratified_sample_by_species(
        records_df, records_dict, max_per_species=max_per_species // 2, seed=42
    )

    # Adjust to target range
    if len(level0_records) < level_0_count[0]:
        print(f"⚠ Level 0: only {len(level0_records)} records (target {level_0_count[0]}-{level_0_count[1]})")
    else:
        print(f"✓ Level 0: {len(level0_records)} records (target {level_0_count[0]}-{level_0_count[1]})")

    # Level 1: all remaining records, capped per species
    level1_df, level1_records = stratified_sample_by_species(
        records_df, records_dict, max_per_species=max_per_species, seed=42
    )

    if len(level1_records) < level_1_count[0]:
        print(f"⚠ Level 1: only {len(level1_records)} records (target {level_1_count[0]}-{level_1_count[1]})")
    elif len(level1_records) > level_1_count[1]:
        print(f"⚠ Level 1: {len(level1_records)} exceeds max {level_1_count[1]} (keeping all stratified)")
    else:
        print(f"✓ Level 1: {len(level1_records)} records (target {level_1_count[0]}-{level_1_count[1]})")

    # Export FASTA
    export_fasta(level0_records, output_dir / "gn_sequences_level0.fasta")
    export_fasta(level1_records, output_dir / "gn_sequences_level1.fasta")

    # Export metadata
    level0_df.to_csv(output_dir / "metadata_level0.tsv", sep="\t", index=False)
    level1_df.to_csv(output_dir / "metadata_level1.tsv", sep="\t", index=False)

    # Export manifests
    export_manifest(level0_df, level0_records, output_dir / "../manifests/dataset_manifest_level0.tsv")
    export_manifest(level1_df, level1_records, output_dir / "../manifests/dataset_manifest_level1.tsv")

    return level0_records, level1_records, level0_df, level1_df
