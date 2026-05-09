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
    Extract Gn (glycoprotein N-terminal) from GenBank record.

    Priority:
    1. Explicitly annotated "Gn" CDS
    2. Explicitly annotated "glycoprotein" CDS (assume first 480 aa)
    3. Explicitly annotated "GPC" / "glycoprotein precursor" (truncate to 480 aa)

    Returns:
        Gn sequence (str) or None if not found
    """
    # Try explicit Gn
    for feature in record.features:
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", [""])[0].lower()
            if "gn" in product or "g1" in product:
                try:
                    seq_str = str(
                        feature.extract(record.seq)
                    )  # This is DNA; we need protein
                    # If this is a protein feature, it's already translated
                    if len(seq_str) > 1000:
                        return None  # Too long for Gn, likely DNA
                    return seq_str
                except:
                    pass

    # Try explicit Gc
    gn_seq = None
    gpc_seq = None

    for feature in record.features:
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", [""])[0].lower()

            # Try to get Gn from translation
            if "translation" in feature.qualifiers:
                translation = feature.qualifiers["translation"][0]

                if "gn" in product or "g1" in product:
                    return translation

                # Collect GPC for fallback
                if "gpc" in product or "glycoprotein precursor" in product:
                    gpc_seq = translation

    # Fallback: truncate GPC to first 480 aa (Gn region)
    if gpc_seq and len(gpc_seq) > 350:
        return gpc_seq[:480]  # Truncate to Gn

    return gn_seq


def fetch_gn_sequences(
    config_path: Path,
    output_dir: Path,
    num_per_taxon: int = 80,  # Download 80, keep 60 post-QC
) -> Tuple[List, Dict]:
    """
    Fetch Gn sequences from NCBI for all configured taxa.

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

        # Query: taxon + Gn/GPC keywords
        query = f"txid{taxon_id}[Organism:exp] AND (Gn[Title] OR GPC[Title] OR 'glycoprotein precursor'[Title])"

        try:
            ids = search_ncbi(query, retmax=num_per_taxon)
            print(f"  Found {len(ids)} records")

            if not ids:
                continue

            # Fetch in batches
            batch_count = 0
            for records, failed in fetch_batch(
                ids,
                batch_size=config["ncbi"]["batch_size"],
                sleep_time=config["ncbi"]["sleep_between_batches"],
            ):
                batch_count += 1
                all_records.extend(records)

                for record in records:
                    gn_seq = extract_gn_from_record(record)
                    if gn_seq:
                        metadata[record.id] = {
                            "organism": species_name,
                            "taxon_id": taxon_id,
                            "description": record.description,
                            "seq_length": len(gn_seq),
                        }

        except Exception as e:
            print(f"  Error: {e}")

    print(f"\n✓ Total records fetched: {len(all_records)}")
    return all_records, metadata
