"""
Extract and normalize metadata from GenBank records.

Handles:
- Parsing organism, country, host, collection_date
- Normalizing country → region → Old/New World
- Parsing collection dates to years
- Categorizing hosts (human, rodent, unknown)
"""

import re
from typing import Dict, Optional

import pandas as pd


# Taxonomic name cleaning
SPECIES_MAP = {
    "Andes virus": "Andes_virus",
    "Hantaan virus": "Hantaan_virus",
    "Puumala virus": "Puumala_virus",
    "Sin Nombre virus": "Sin_Nombre_virus",
    "Seoul virus": "Seoul_virus",
    "Dobrava virus": "Dobrava_virus",
    "Maporal virus": "Maporal_virus",
    "Tula virus": "Tula_virus",
    "Prospect Hill virus": "Prospect_Hill_virus",
}

# Country → Region → Old/New World mapping
COUNTRY_REGION_MAP = {
    # North America
    "United States": {"region": "North America", "world": "New World"},
    "USA": {"region": "North America", "world": "New World"},
    "Canada": {"region": "North America", "world": "New World"},
    "Mexico": {"region": "Central America", "world": "New World"},
    # Central/South America
    "Argentina": {"region": "South America", "world": "New World"},
    "Chile": {"region": "South America", "world": "New World"},
    "Brazil": {"region": "South America", "world": "New World"},
    "Peru": {"region": "South America", "world": "New World"},
    "Colombia": {"region": "South America", "world": "New World"},
    # Europe
    "Finland": {"region": "Scandinavia", "world": "Old World"},
    "Sweden": {"region": "Scandinavia", "world": "Old World"},
    "Norway": {"region": "Scandinavia", "world": "Old World"},
    "Germany": {"region": "Central Europe", "world": "Old World"},
    "France": {"region": "Western Europe", "world": "Old World"},
    "Belgium": {"region": "Western Europe", "world": "Old World"},
    "United Kingdom": {"region": "Western Europe", "world": "Old World"},
    "UK": {"region": "Western Europe", "world": "Old World"},
    "Spain": {"region": "Southern Europe", "world": "Old World"},
    "Serbia": {"region": "Eastern Europe", "world": "Old World"},
    "Slovenia": {"region": "Eastern Europe", "world": "Old World"},
    "Croatia": {"region": "Eastern Europe", "world": "Old World"},
    "Russia": {"region": "Eastern Europe", "world": "Old World"},
    # Asia
    "China": {"region": "East Asia", "world": "Old World"},
    "Japan": {"region": "East Asia", "world": "Old World"},
    "South Korea": {"region": "East Asia", "world": "Old World"},
    "Korea": {"region": "East Asia", "world": "Old World"},
    "Thailand": {"region": "Southeast Asia", "world": "Old World"},
    "Vietnam": {"region": "Southeast Asia", "world": "Old World"},
    # Africa
    "Kenya": {"region": "East Africa", "world": "Old World"},
    # Oceania
    "Australia": {"region": "Oceania", "world": "Old World"},
}

# Host categorization
HOST_MAP = {
    "human": ["Homo sapiens", "human"],
    "rodent": [
        "Peromyscus",
        "Microtus",
        "Rattus",
        "Myodes",
        "Apodemus",
        "Lemmus",
        "Neotoma",
        "Sigmodon",
    ],
    "unknown": ["unknown", "not determined", ""],
}


def clean_species_name(organism_str: str) -> str:
    """Normalize species name from organism string."""
    organism_str = organism_str.strip()
    for key, val in SPECIES_MAP.items():
        if key.lower() in organism_str.lower():
            return val
    return organism_str.replace(" ", "_")


def parse_country(country_str: Optional[str]) -> tuple:
    """
    Parse country from GenBank country field.
    Format: "Country: Argentina: Patagonia" or "Argentina"
    Returns: (country_raw, country_norm, region, old_new_world)
    """
    if not country_str or country_str.lower() == "unknown":
        return "unknown", "unknown", "unknown", "unknown"

    # Remove "Country: " prefix if present
    country_str = country_str.replace("Country: ", "").strip()

    # Extract main country (before first colon if exists)
    parts = country_str.split(":")
    main_country = parts[0].strip()

    # Try to match in map
    for key, mapping in COUNTRY_REGION_MAP.items():
        if key.lower() in main_country.lower() or main_country.lower() in key.lower():
            return country_str, main_country, mapping["region"], mapping["world"]

    # Fallback
    return country_str, main_country, "unknown", "unknown"


def parse_year(date_str: Optional[str]) -> Optional[int]:
    """
    Parse year from collection date.
    Formats: "2018", "2018-03", "2018-03-15", "2018-09"
    Returns: year (int) or None if invalid
    """
    if not date_str:
        return None

    date_str = date_str.strip()

    # Try YYYY
    if re.match(r"^\d{4}$", date_str):
        year = int(date_str)
        if 1970 <= year <= 2025:
            return year

    # Try YYYY-MM or YYYY-MM-DD
    match = re.match(r"^(\d{4})-", date_str)
    if match:
        year = int(match.group(1))
        if 1970 <= year <= 2025:
            return year

    return None


def categorize_host(host_str: Optional[str]) -> str:
    """Categorize host from host field."""
    if not host_str:
        return "unknown"

    host_lower = host_str.lower()

    for category, keywords in HOST_MAP.items():
        for keyword in keywords:
            if keyword.lower() in host_lower:
                return category

    return "unknown"


def extract_metadata_from_record(record, seq_str: str, species_name: str) -> Dict:
    """
    Extract all metadata fields from a GenBank record.

    Args:
        record: Bio.SeqIO record
        seq_str: Gn sequence string
        species_name: Species name (e.g., "Andes_virus")

    Returns:
        dict with all metadata fields
    """
    # Extract organism from features
    organism = species_name
    taxon_id = None
    country_raw = None
    host_raw = None
    collection_date = None

    for feature in record.features:
        if feature.type == "source":
            organism = (
                feature.qualifiers.get("organism", [species_name])[0]
                or species_name
            )
            # Get taxon ID
            db_xref = feature.qualifiers.get("db_xref", [])
            for ref in db_xref:
                if "taxon" in ref.lower():
                    taxon_id = ref.split(":")[-1]
            country_raw = feature.qualifiers.get("country", [None])[0]
            host_raw = feature.qualifiers.get("host", [None])[0]
            collection_date = feature.qualifiers.get("collection_date", [None])[0]

    # Normalize fields
    country_raw, country_norm, region, old_new_world = parse_country(country_raw)
    year = parse_year(collection_date)
    host_category = categorize_host(host_raw)
    species_clean = clean_species_name(organism)

    return {
        "accession": record.id,
        "description": record.description,
        "organism": organism,
        "taxon_id": taxon_id,
        "species_clean": species_clean,
        "country_raw": country_raw,
        "country_norm": country_norm,
        "region": region,
        "old_new_world": old_new_world,
        "host_raw": host_raw or "unknown",
        "host_category": host_category,
        "collection_date": collection_date,
        "year": year,
        "seq_length": len(seq_str),
        "extraction_method": "gn_annotated",  # Will be updated during QC
    }


def build_metadata_df(records_with_sequences: Dict[str, Dict]) -> pd.DataFrame:
    """
    Build metadata DataFrame from dict of {accession: {seq, metadata}}.

    Args:
        records_with_sequences: dict with accession → {seq, metadata}

    Returns:
        pandas DataFrame
    """
    rows = []
    for accession, data in records_with_sequences.items():
        row = data["metadata"].copy()
        rows.append(row)

    df = pd.DataFrame(rows)

    # Ensure required columns exist with defaults
    required_cols = [
        "accession",
        "organism",
        "taxon_id",
        "species_clean",
        "country_raw",
        "country_norm",
        "region",
        "old_new_world",
        "host_raw",
        "host_category",
        "year",
        "seq_length",
        "extraction_method",
    ]

    for col in required_cols:
        if col not in df.columns:
            if col == "species_clean":
                df[col] = df.get("organism", "unknown").str.replace(" ", "_")
            elif col == "region":
                df[col] = "unknown"
            elif col == "old_new_world":
                df[col] = "unknown"
            else:
                df[col] = None

    return df[required_cols]
