#!/usr/bin/env python3
"""
Phase 1 Smoke Test: Process synthetic GenBank data (no NCBI calls).

This demonstrates the full Phase 1 pipeline without depending on NCBI connectivity.
Uses synthetic GenBank records with real-like Gn sequences and metadata.
"""

import json
import sys
import logging
from pathlib import Path
from io import StringIO

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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


def create_synthetic_genbank_records():
    """
    Create synthetic GenBank records with realistic Gn sequences.
    """
    # Real Gn sequences (truncated for testing)
    sequences = {
        "HantaanGn_seq1": (
            "MFGLLFFGLVCFAYCEDKKTASVKLALYGTDKVSGGYTQNGWVTRPMWSTTFNMN"
            "TKKNGDNDSVTVDWTSGWQYDGSRSVPFPLSDDKNKKTLSQFDKMHWFSYKNM"
            "NQMVKGSIGFIVDDQHPVGYFVSKRLFRGPYTDKYQCFQNVDFPWKKRLGHYP"
            "TIQFEHDATNVPKFWVYEPKGVDFAEEWVSVPSTGAAWHVTSKPDATDNGFQ"
            "SKRWHVGKTVLQVHQKKDLFRDLTSHSPVGGPPVTAVKAELSGVQSVSTAMTR"
            "GDSVTEEPGSSGPGVPLWKPLSVWKCWSPPFGVPNQPNLKCSQQQEEGVGWFD"
            "YSKQVPPKPHPSVPSYQEGVPHTSDSQGDKGADWSMQYGDMVQTTYQCNDAPT"
            "VSLVR"
        ),
        "HantaanGn_seq2": (
            "MFGLLFFGLVCFAYCEDKKTASVKLALYGTDKVSGGYTQNGWVTRPMWSTTFNMN"
            "TKKNGDNDSVTVDWTSGWQYDGSRSVPFPLSDDKNKKTLSQFDKMHWFSYKNM"
            "NQMVKGSIGFIVDDQHPVGYFVSKRLFRGPYTDKYQCFQNVDFPWKKRLGHYP"
            "TIQFEHDATNVPKFWVYEPKGVDFAEEWVSVPSTGAAWHVTSKPDATDNGFQ"
            "SKRWHVGKTVLQVHQKKDLFRDLTSHSPVGGPPVTAVKAELSGVQSVSTAMTR"
            "GDSVTEEPGSSGPGVPLWKPLSVWKCWSPPFGVPNQPNLKCSQQQEEGVGWFD"
            "YSKQVPPKPHPSVPSYQEGVPHTSDSQGDKGADWSMQYGDMVQTTYQCNDAPT"
            "VSLVR"
        ),  # Exact duplicate
        "PuumalaGn_seq1": (
            "MFGLLFFGLVCFAYCEDKKTASVKLALYGTDKVSGGYTQNGWVTRPMWSTTFNMN"
            "TKKNGDNDSVTVDWTSGWQYDGSRSVPFPLSDDKNKKTLSQFDKMHWFSYKNM"
            "NQMVKGSIGFIVDDQHPVGYFVSKRLFRGPYTDKYQCFQNVDFPWKKRLGHYP"
            "TIQFEHDATNVPKFWVYEPKGVDFAEEWVSVPSTGAAWHVTSKPDATDNGFQ"
            "SKRWHVGKTVLQVHQKKDLFRDLTSHSPVGGPPVTAVKAELSGVQSVSTAMTR"
            "GDSVTEEPGSSGPGVPLWKPLSVWKCWSPPFGVPNQPNLKCSQQQEEGVGWFD"
            "YSKQVPPKPHPSVPSYQEGVPHTSDSQGDKGADWSMQYGDMVQTTYQCNDAPT"
            "VSLVQ"  # One AA different
        ),
        "SeoulGn_seq1": (
            "MFGLLFFGLVCFAYCEDKKTASVKLALYGTDKVSGGYTQNGWVTRPMWSTTFNMN"
            "TKKNGDNDSVTVDWTSGWQYDGSRSVPFPLSDDKNKKTLSQFDKMHWFSYKNM"
            "NQMVKGSIGFIVDDQHPVGYFVSKRLFRGPYTDKYQCFQNVDFPWKKRLGHYP"
            "TIQFEHDATNVPKFWVYEPKGVDFAEEWVSVPSTGAAWHVTSKPDATDNGFQ"
            "SKRWHVGKTVLQVHQKKDLFRDLTSHSPVGGPPVTAVKAELSGVQSVSTAMTR"
            "GDSVTEEPGSSGPGVPLWKPLSVWKCWSPPFGVPNQPNLKCSQQQEEGVGWFD"
            "YSKQVPPKPHPSVPSYQEGVPHTSDSQGDKGADWSMQYGDMVQTTYQCNDAPT"
            "VSLVR"[:420]  # Truncated (too short)
        ),
    }

    # GenBank records
    records = []

    for i, (seqid, seq) in enumerate(sequences.items()):
        species = seqid.split("_")[0]
        species_map = {
            "Hantaan": ("Hantaan virus", 11603),
            "Puumala": ("Puumala virus", 64005),
            "Seoul": ("Seoul virus", 11604),
        }

        organism, taxid = species_map.get(species, ("Unknown", 0))

        # Create feature with metadata
        record = SeqRecord(
            Seq(seq),
            id=f"SRR12345{i:02d}",
            description=f"Orthohantavirus {species}",
        )

        # Add source feature
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        source = SeqFeature(
            FeatureLocation(0, len(seq)),
            type="source",
            qualifiers={
                "organism": [organism],
                "db_xref": [f"taxon:{taxid}"],
                "country": ["China"] if species == "Seoul" else ["Finland"],
                "host": ["Homo sapiens"],
                "collection_date": ["2019-01-15"],
            },
        )

        # Add CDS feature with translation
        cds = SeqFeature(
            FeatureLocation(0, len(seq)),
            type="CDS",
            qualifiers={
                "product": ["glycoprotein Gn"],
                "translation": [seq],
            },
        )

        record.features = [source, cds]
        records.append(record)

    return records


def main():
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    log_dir = Path(__file__).parent.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "smoke_test.log"),
            logging.StreamHandler(),
        ],
    )
    logger = logging.getLogger(__name__)

    logger.info("="*80)
    logger.info("PHASE 1 SMOKE TEST (synthetic data, no NCBI)")
    logger.info("="*80)

    import yaml
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Create synthetic records
    logger.info("\n1. CREATE SYNTHETIC DATA")
    records = create_synthetic_genbank_records()
    logger.info(f"✓ Created {len(records)} synthetic GenBank records")

    # Extract Gn and metadata
    logger.info("\n2. EXTRACT Gn & METADATA")
    records_with_seqs = {}

    for record in records:
        gn_seq = None
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                gn_seq = feature.qualifiers["translation"][0]
                break

        if gn_seq and len(gn_seq) > 300:
            # Get species from source
            species = "unknown"
            for feature in record.features:
                if feature.type == "source":
                    organism = feature.qualifiers.get("organism", ["unknown"])[0]
                    if "Hantaan" in organism:
                        species = "Hantaan_virus"
                    elif "Puumala" in organism:
                        species = "Puumala_virus"
                    elif "Seoul" in organism:
                        species = "Seoul_virus"

            metadata = extract_metadata_from_record(record, gn_seq, species)
            records_with_seqs[record.id] = {
                "seq": gn_seq,
                "metadata": metadata,
            }

    logger.info(f"✓ Extracted {len(records_with_seqs)} sequences")

    # Apply QC
    logger.info("\n3. APPLY QC FILTERS")
    passed_qc, failed_qc = apply_qc_filters(records_with_seqs, config["qc"])
    logger.info(f"  Passed: {len(passed_qc)}")
    logger.info(f"  Failed: {len(failed_qc)}")
    if failed_qc:
        logger.info(f"  Failed reasons: {list(failed_qc.values())}")

    # Remove duplicates
    logger.info("\n4. REMOVE DUPLICATES")
    unique, exact_dups = remove_exact_duplicates(passed_qc)
    logger.info(f"  Unique: {len(unique)}")
    logger.info(f"  Exact dups removed: {len(exact_dups)}")

    final, near_dups = remove_near_duplicates(unique, identity_threshold=0.99)
    logger.info(f"  Final after near-dups: {len(final)}")
    logger.info(f"  Near-dups removed: {len(near_dups)}")

    # Build metadata df
    logger.info("\n5. BUILD METADATA")
    metadata_df = build_metadata_df(final)
    logger.info(f"  Metadata: {metadata_df.shape}")

    # Generate report
    logger.info("\n6. GENERATE QC REPORT")
    qc_report = generate_qc_report(
        len(records),
        len(final),
        failed_qc,
        exact_dups,
        near_dups,
        metadata_df,
    )
    logger.info(json.dumps(qc_report, indent=2))

    # Build splits
    logger.info("\n7. BUILD LEVEL 0 & LEVEL 1")
    processed_dir = Path(__file__).parent.parent / "data" / "processed"
    build_level_0_and_1(
        final,
        metadata_df,
        processed_dir,
        level_0_count=(2, 20),
        level_1_count=(2, 300),
        max_per_species=10,
    )

    # Save report
    manifests_dir = Path(__file__).parent.parent / "results" / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    with open(manifests_dir / "qc_report_smoke_test.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    logger.info("\n" + "="*80)
    logger.info("✓ SMOKE TEST PASSED")
    logger.info("="*80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
