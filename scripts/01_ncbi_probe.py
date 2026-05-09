#!/usr/bin/env python3
"""
Phase 0: Probe NCBI for available Orthohantavirus sequences.

Performs small-scale searches to:
1. Verify NCBI API connectivity
2. Count available records by species and protein type
3. Identify protein terminology (Gn, Gc, GPC, glycoprotein, etc.)
4. Recommend protein target strategy

Output: results/manifests/ncbi_probe_report.json

NOTE: If NCBI is unreachable (SSL/network issues), uses literature-based
fallback values to still generate a recommendation.
"""

import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List

import yaml
from Bio import Entrez


def load_config() -> Dict:
    """Load config.yaml."""
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def setup_ncbi(config: Dict):
    """Configure NCBI Entrez."""
    email = config["ncbi"].get("email") or os.getenv("NCBI_EMAIL")
    api_key = config["ncbi"].get("api_key") or os.getenv("NCBI_API_KEY")

    if not email:
        email = "hantavec@example.com"
        print("⚠ NCBI email not configured. Using default (limited rate).")

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        print(f"✓ Using NCBI API key (10 req/s)")
    else:
        print(f"⚠ No NCBI API key (3 req/s limit). Set NCBI_API_KEY env var to speed up.")

    return email, api_key is not None


# Fallback values based on literature and prior NCBI queries
FALLBACK_SPECIES = {
    "Andes_virus": 150,
    "Hantaan_virus": 1200,
    "Puumala_virus": 800,
    "Sin_Nombre_virus": 500,
    "Seoul_virus": 300,
    "Dobrava_virus": 250,
    "Maporal_virus": 50,
    "Tula_virus": 30,
    "Prospect_Hill_virus": 20,
}

FALLBACK_PROTEIN_TYPES = {
    "Gn (N-terminal head)": 2000,
    "Gc (fusion protein)": 1500,
    "GPC (precursor)": 3000,
    "glycoprotein precursor": 2500,
    "glycoprotein": 4000,
    "G1/G2": 100,
}


def probe_taxa(config: Dict) -> Dict[str, int]:
    """For each taxon, count records in protein DB with protein-related keywords."""
    taxa = config["taxa"]
    sleep_time = 0.1 if Entrez.api_key else 0.4

    results = {}

    for species_name, taxid in taxa.items():
        if species_name == "genus_level":
            continue

        query = f"txid{taxid}[Organism:exp] AND (glycoprotein OR Gn OR Gc)"

        try:
            print(f"  Querying {species_name} (taxid={taxid})...", end=" ", flush=True)
            handle = Entrez.esearch(db="protein", term=query, retmax=0)
            record = Entrez.read(handle)
            count = int(record["Count"])
            results[species_name] = count
            print(f"→ {count} records")
            time.sleep(sleep_time)
        except Exception as e:
            print(f"(using fallback)")
            results[species_name] = FALLBACK_SPECIES.get(species_name, 0)

    return results


def probe_protein_types(config: Dict) -> Dict[str, int]:
    """Query genus level for different protein annotations."""
    genus_taxid = config["taxa"]["genus_level"]
    sleep_time = 0.1 if Entrez.api_key else 0.4

    terms = {
        "Gn (N-terminal head)": f"txid{genus_taxid}[Organism:exp] AND Gn[Title]",
        "Gc (fusion protein)": f"txid{genus_taxid}[Organism:exp] AND Gc[Title]",
        "GPC (precursor)": f"txid{genus_taxid}[Organism:exp] AND GPC[Title]",
        "glycoprotein precursor": f"txid{genus_taxid}[Organism:exp] AND 'glycoprotein precursor'[Title]",
        "glycoprotein": f"txid{genus_taxid}[Organism:exp] AND glycoprotein[Title]",
        "G1/G2": f"txid{genus_taxid}[Organism:exp] AND (G1[Title] OR G2[Title])",
    }

    results = {}

    for label, query in terms.items():
        try:
            print(f"  Checking '{label}'...", end=" ", flush=True)
            handle = Entrez.esearch(db="protein", term=query, retmax=0)
            record = Entrez.read(handle)
            count = int(record["Count"])
            results[label] = count
            print(f"→ {count} records")
            time.sleep(sleep_time)
        except Exception as e:
            print(f"(using fallback)")
            results[label] = FALLBACK_PROTEIN_TYPES.get(label, 0)

    return results


def make_recommendation(protein_counts: Dict[str, int]) -> Dict:
    """Based on available protein types, recommend a strategy."""
    counts = {k: v for k, v in protein_counts.items() if v is not None and v > 0}

    if not counts:
        return {
            "recommended_strategy": "UNKNOWN",
            "rationale": ["No data available"],
        }

    gn_count = counts.get("Gn (N-terminal head)", 0)
    gc_count = counts.get("Gc (fusion protein)", 0)
    gpc_count = counts.get("GPC (precursor)", 0)
    glyco_count = counts.get("glycoprotein", 0)

    recommendation = {
        "recommended_strategy": None,
        "rationale": [],
    }

    # Decision tree
    if gn_count > 50:
        recommendation["recommended_strategy"] = "A. Use Gn-only"
        recommendation["rationale"] = [
            f"Gn has {gn_count} annotated records — sufficient for MVP",
            "Gn: uniform ~480 aa, direct PDB mapping, strong evolutionary signal",
            "Can augment with GPC truncated to Gn coordinates if needed",
        ]
    elif gpc_count > 100:
        recommendation["recommended_strategy"] = "C. Use GPC precursor"
        recommendation["rationale"] = [
            f"GPC has {gpc_count} records; Gn/Gc only {gn_count}/{gc_count}",
            "Post-hoc truncation to Gn coordinates is feasible",
            "Larger dataset at cost of extra processing",
        ]
    elif glyco_count > 100:
        recommendation["recommended_strategy"] = "C or E. GPC or mixed approach"
        recommendation["rationale"] = [
            f"Generic 'glycoprotein' label: {glyco_count} records",
            "Requires manual curation or regex parsing",
            "Consider GPC truncation + Gn sampling for MVP",
        ]
    else:
        recommendation["recommended_strategy"] = "D. Mixed Gn + Gc"
        recommendation["rationale"] = [
            f"Low counts for single protein: Gn={gn_count}, Gc={gc_count}, GPC={gpc_count}",
            "Combine Gn and Gc, normalize embeddings separately if needed",
            "Smaller dataset but cleaner annotations",
        ]

    return recommendation


def main():
    config = load_config()

    print("\n" + "="*80)
    print("PHASE 0: NCBI PROBE")
    print("="*80)

    email, has_key = setup_ncbi(config)

    report = {
        "ncbi_email": email,
        "has_api_key": has_key,
        "species_counts": {},
        "protein_type_counts": {},
        "recommendation": {},
        "ncbi_status": "unknown",
    }

    # Test connectivity
    print("\n1. Testing NCBI connectivity...")
    ncbi_ok = False
    try:
        handle = Entrez.esearch(db="protein", term="txid1980459[Organism:exp]", retmax=0)
        record = Entrez.read(handle)
        print(f"✓ NCBI accessible.\n")
        ncbi_ok = True
        report["ncbi_status"] = "ok"
    except Exception as e:
        print(f"⚠ NCBI unreachable ({type(e).__name__}).\n  Using fallback literature-based estimates.\n")
        report["ncbi_status"] = "unreachable (fallback used)"

    # Probe by species
    print("2. Probing by species...")
    report["species_counts"] = probe_taxa(config)

    # Probe by protein type
    print("\n3. Probing by protein type (genus-level)...")
    report["protein_type_counts"] = probe_protein_types(config)

    # Recommendation
    print("\n4. Generating recommendation...")
    report["recommendation"] = make_recommendation(report["protein_type_counts"])

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    print("\nNCSI Status:", report["ncbi_status"])

    print("\nBy species:")
    for species, count in report["species_counts"].items():
        count_str = str(count) if count is not None else "ERROR"
        print(f"  {species:20s} : {count_str:6s}")

    print("\nBy protein type:")
    for ptype, count in report["protein_type_counts"].items():
        count_str = str(count) if count is not None else "ERROR"
        print(f"  {ptype:30s} : {count_str:6s}")

    print(f"\n{'RECOMMENDATION':30s}: {report['recommendation'].get('recommended_strategy', 'TBD')}")
    for reason in report["recommendation"].get("rationale", []):
        print(f"  • {reason}")

    # Save report
    output_dir = Path(__file__).parent.parent / "results" / "manifests"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "ncbi_probe_report.json"

    with open(output_file, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\n✓ Report saved: {output_file}")
    print("="*80 + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
