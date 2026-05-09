#!/usr/bin/env python3
"""
Phase 1: Build curated dataset from NCBI sequences.

Steps:
1. Fetch Gn sequences from NCBI for all taxa
2. Extract Gn protein and metadata
3. Apply QC filters (length, ambiguous AA, stops, duplicates)
4. Create Level 0 (10-20) and Level 1 (100-300) splits
5. Generate reports

Outputs:
- data/processed/gn_sequences_level0.fasta
- data/processed/gn_sequences_level1.fasta
- data/processed/metadata_level0.tsv
- data/processed/metadata_level1.tsv
- results/manifests/dataset_manifest_level{0,1}.tsv
- results/manifests/qc_report.json
"""

import json
import sys
import logging
from pathlib import Path
from typing import Dict

import yaml
from Bio import Entrez, SeqIO

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from data.fetch import (
    load_config,
    setup_ncbi,
    search_ncbi,
    fetch_batch,
    extract_gn_from_record,
)
from data.metadata import (
    extract_metadata_from_record,
    build_metadata_df,
)
from data.qc import (
    apply_qc_filters,
    remove_exact_duplicates,
    remove_near_duplicates,
    generate_qc_report,
)
from data.splits import build_level_0_and_1


def setup_logging(log_dir: Path) -> logging.Logger:
    """Setup logging."""
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "phase1.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def fetch_and_extract_gn(
    config: Dict,
    logger: logging.Logger,
    num_per_taxon: int = 80,
) -> tuple:
    """
    Fetch and extract Gn sequences for all taxa.

    Returns:
        (records_with_seqs, total_fetched)
    """
    setup_ncbi(config["ncbi"].get("email"), config["ncbi"].get("api_key"))

    records_with_seqs = {}
    total_fetched = 0

    for species_name, taxon_id in config["taxa"].items():
        if species_name == "genus_level":
            continue

        logger.info(f"{species_name} (taxid={taxon_id}):")

        query = f"txid{taxon_id}[Organism:exp] AND (Gn[Title] OR GPC[Title] OR 'glycoprotein precursor'[Title])"

        try:
            ids = search_ncbi(query, retmax=num_per_taxon)
            logger.info(f"  Found {len(ids)} records")

            if not ids:
                continue

            # Fetch in batches
            for records, failed in fetch_batch(
                ids,
                batch_size=config["ncbi"]["batch_size"],
                sleep_time=config["ncbi"]["sleep_between_batches"],
            ):
                for record in records:
                    gn_seq = extract_gn_from_record(record)
                    if gn_seq and len(gn_seq) > 300:  # Rough min check
                        metadata = extract_metadata_from_record(
                            record, gn_seq, species_name
                        )
                        records_with_seqs[record.id] = {
                            "seq": gn_seq,
                            "metadata": metadata,
                        }
                        total_fetched += 1

        except Exception as e:
            logger.error(f"  Error: {e}")

    logger.info(f"✓ Total Gn sequences extracted: {len(records_with_seqs)}")
    return records_with_seqs, total_fetched


def main():
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    log_dir = Path(__file__).parent.parent / "logs"
    logger = setup_logging(log_dir)

    logger.info("="*80)
    logger.info("PHASE 1: BUILD DATASET")
    logger.info("="*80)

    config = load_config(config_path)

    # Step 1: Fetch and extract
    logger.info("\n1. FETCH & EXTRACT Gn sequences from NCBI")
    records_with_seqs, total_fetched = fetch_and_extract_gn(config, logger, num_per_taxon=80)

    if not records_with_seqs:
        logger.error("No sequences extracted. Aborting.")
        return 1

    # Step 2: Apply QC
    logger.info("\n2. APPLY QC FILTERS")
    passed_qc, failed_qc = apply_qc_filters(records_with_seqs, config["qc"])
    logger.info(f"  Passed QC: {len(passed_qc)}")
    logger.info(f"  Failed QC: {len(failed_qc)}")

    # Step 3: Remove exact duplicates
    logger.info("\n3. REMOVE EXACT DUPLICATES")
    unique_seqs, exact_dups = remove_exact_duplicates(passed_qc)
    logger.info(f"  Unique: {len(unique_seqs)}")
    logger.info(f"  Exact duplicates removed: {len(exact_dups)}")

    # Step 4: Remove near-duplicates
    logger.info("\n4. REMOVE NEAR-DUPLICATES (≥99% identity)")
    final_seqs, near_dups = remove_near_duplicates(unique_seqs, identity_threshold=0.99)
    logger.info(f"  Final: {len(final_seqs)}")
    logger.info(f"  Near-duplicates removed: {len(near_dups)}")

    # Step 5: Build metadata dataframe
    logger.info("\n5. BUILD METADATA DATAFRAME")
    metadata_df = build_metadata_df(final_seqs)
    logger.info(f"  Metadata shape: {metadata_df.shape}")

    # Step 6: Generate QC report
    logger.info("\n6. GENERATE QC REPORT")
    qc_report = generate_qc_report(
        total_fetched,
        len(final_seqs),
        failed_qc,
        exact_dups,
        near_dups,
        metadata_df,
    )

    manifests_dir = Path(__file__).parent.parent / "results" / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    with open(manifests_dir / "qc_report.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    logger.info(json.dumps(qc_report, indent=2))

    # Step 7: Build Level 0 and Level 1
    logger.info("\n7. BUILD LEVEL 0 & LEVEL 1 SPLITS")
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    build_level_0_and_1(
        final_seqs,
        metadata_df,
        processed_dir,
        level_0_count=(10, 20),
        level_1_count=(100, 300),
        max_per_species=60,
    )

    logger.info("\n" + "="*80)
    logger.info("✓ PHASE 1 COMPLETE")
    logger.info("="*80)

    # Check if passed acceptance criteria
    if qc_report["passed_qc"] < 100:
        logger.warning("⚠ Level 1 has < 100 records. Recommend manual augmentation.")
        return 1

    if qc_report["old_world_fraction"] > 0.8 or qc_report["new_world_fraction"] > 0.8:
        logger.warning("⚠ Heavy bias toward Old/New World. Per-species cap may have failed.")
        return 1

    logger.info("✓ Dataset ready for Phase 2 (embeddings)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
