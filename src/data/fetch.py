"""
Fetch Orthohantavirus Gn sequences from NCBI.

Handles:
- Entrez esearch/efetch queries
- Batching and retries with exponential backoff
- Checkpoint system (resume on interrupt)
- GenBank parsing
"""

import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml
from Bio import Entrez, SeqIO


def load_config(config_path: Path) -> Dict:
    """Load config.yaml."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def setup_logging(log_dir: Path) -> logging.Logger:
    """Setup logging to file and console."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"fetch_{int(time.time())}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def setup_ncbi(email: str, api_key: Optional[str] = None):
    """Configure NCBI Entrez."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def load_checkpoint(checkpoint_path: Path) -> Dict[str, str]:
    """Load fetch status checkpoint."""
    if checkpoint_path.exists():
        with open(checkpoint_path) as f:
            return json.load(f)
    return {}


def save_checkpoint(checkpoint_path: Path, status: Dict[str, str]):
    """Save fetch status checkpoint."""
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    with open(checkpoint_path, "w") as f:
        json.dump(status, f, indent=2)


def search_ncbi(
    query: str,
    db: str = "protein",
    retmax: int = 10000,
) -> List[str]:
    """
    Search NCBI and return list of IDs.

    Args:
        query: Entrez query string
        db: database (protein, nucleotide)
        retmax: max results to return

    Returns:
        List of sequence IDs
    """
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        record = Entrez.read(handle)
        return record.get("IdList", [])
    except Exception as e:
        raise RuntimeError(f"NCBI search failed: {e}")


def fetch_batch(
    ids: List[str],
    db: str = "protein",
    batch_size: int = 200,
    sleep_time: float = 0.4,
    max_retries: int = 3,
    backoff: List[int] = [1, 4, 16],
) -> Tuple[List, List]:
    """
    Fetch sequence batch from NCBI with retries.

    Args:
        ids: List of sequence IDs
        db: database
        batch_size: records per batch
        sleep_time: sleep between batches
        max_retries: retry attempts
        backoff: seconds to wait (exponential)

    Yields:
        (successful_records, failed_ids)
    """
    failed_ids = []

    for i in range(0, len(ids), batch_size):
        batch_ids = ids[i : i + batch_size]
        batch_num = i // batch_size + 1

        retry_count = 0
        while retry_count < max_retries:
            try:
                print(
                    f"  Batch {batch_num}: fetching {len(batch_ids)} records...",
                    end=" ",
                    flush=True,
                )
                handle = Entrez.efetch(
                    db=db,
                    id=",".join(batch_ids),
                    rettype="gb",
                    retmode="text",
                )
                records = list(SeqIO.parse(handle, "genbank"))
                print(f"✓ {len(records)} records")
                time.sleep(sleep_time)
                yield records, []
                break

            except Exception as e:
                retry_count += 1
                if retry_count < max_retries:
                    wait = backoff[retry_count - 1]
                    print(f"retry in {wait}s ({e})")
                    time.sleep(wait)
                else:
                    print(f"FAILED after {max_retries} retries")
                    failed_ids.extend(batch_ids)
                    yield [], batch_ids


def extract_gn_from_record(record) -> Optional[str]:
    """
    Extract Gn (glycoprotein N-terminal, ~480 aa) from GenBank record (nucleotide).

    Priority (Tier 1 → Tier 2 → Tier 3):

    Tier 1: Explicitly annotated Gn/G1 with translation
      - "Gn", "gn", "G1", "g1"
      - "glycoprotein Gn", "glycoprotein G1"
      - "envelope glycoprotein Gn", "surface glycoprotein Gn"

    Tier 2: GPC / glycoprotein precursor with translation → truncate to first 480 aa
      - "GPC", "glycoprotein precursor", "M segment"
      - Returns translation[0:480]

    Tier 3: Generic "envelope glycoprotein" → truncate to 480 aa

    Returns:
        Gn sequence (str, ~350-600 aa) or None if not found
    """
    if not record.features:
        return None

    # Tier 1: Explicitly annotated Gn/G1
    gn_qualifiers = [
        "gn", "g1",
        "glycoprotein gn", "glycoprotein g1",
        "envelope glycoprotein gn",
        "surface glycoprotein gn",
    ]

    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            product = feature.qualifiers.get("product", [""])[0].lower()

            if any(qual in product for qual in gn_qualifiers):
                seq = feature.qualifiers["translation"][0]
                if 300 < len(seq) < 800:
                    return seq

    # Tier 2: GPC / glycoprotein precursor
    gpc_seq = None
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            product = feature.qualifiers.get("product", [""])[0].lower()
            translation = feature.qualifiers["translation"][0]

            if any(term in product for term in ["gpc", "glycoprotein precursor"]):
                # If right size for GPC (~1000-1300 aa), truncate to Gn
                if 900 < len(translation) < 1400:
                    return translation[:480]
                gpc_seq = translation

            # Also check for generic "envelope glycoprotein"
            elif "envelope glycoprotein" in product and not any(q in product for q in ["gc", "g2"]):
                if 900 < len(translation) < 1400:
                    return translation[:480]

    # Tier 3: Fallback to GPC if size is suitable
    if gpc_seq:
        if 900 < len(gpc_seq) < 1400:
            return gpc_seq[:480]
        elif len(gpc_seq) > 480:
            return gpc_seq[:480]

    return None


def fetch_gn_sequences(
    config_path: Path,
    output_dir: Path,
    num_per_taxon: int = 80,  # Download 80, keep 60 post-QC
) -> Tuple[List, Dict]:
    """
    Fetch Gn sequences from NCBI for all configured taxa.

    Searches nucleotide database for M-segment (genome segment M) records,
    extracts the GPC (glycoprotein precursor) protein via translation.

    Args:
        config_path: Path to config.yaml
        output_dir: Where to save .gb files
        num_per_taxon: Max sequences to download per taxon

    Returns:
        (list_of_records, metadata_dict)
    """
    config = load_config(config_path)
    setup_ncbi(
        config["ncbi"].get("email"),
        config["ncbi"].get("api_key"),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    all_records = []
    metadata = {}

    for species_name, taxon_id in config["taxa"].items():
        if species_name == "genus_level":
            continue

        print(f"\n{species_name} (taxid={taxon_id}):")

        # Query: taxon + glycoprotein precursor in nucleotide database
        # M-segment encodes the glycoprotein precursor (GPC, ~1140 aa)
        # We extract Gn (first ~480 aa) from the translated sequence
        query = f"txid{taxon_id}[Organism:exp] AND ('glycoprotein precursor' OR 'envelope glycoprotein' OR Gn)"

        try:
            ids = search_ncbi(query, db="nucleotide", retmax=num_per_taxon)
            print(f"  Found {len(ids)} records")

            if not ids:
                continue

            # Fetch in batches (nucleotide database)
            batch_count = 0
            for records, failed in fetch_batch(
                ids,
                db="nucleotide",
                batch_size=config["ncbi"]["batch_size"],
                sleep_time=config["ncbi"]["sleep_between_batches"],
            ):
                batch_count += 1
                all_records.extend(records)

                for record in records:
                    gn_seq = extract_gn_from_record(record)
                    if gn_seq and 350 <= len(gn_seq) <= 600:
                        # Store both seq and metadata for later processing
                        metadata[record.id] = {
                            "seq": gn_seq,
                            "metadata": {
                                "organism": species_name,
                                "taxon_id": taxon_id,
                                "description": record.description,
                                "seq_length": len(gn_seq),
                            }
                        }

        except Exception as e:
            print(f"  Error: {e}")

    print(f"\n✓ Total Gn sequences extracted: {len(metadata)}")
    return metadata
