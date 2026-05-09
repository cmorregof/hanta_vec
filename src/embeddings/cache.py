"""
Cache management for embeddings.

Builds auditable cache index for reproducibility.
"""

import json
from pathlib import Path


def build_cache_index(config: dict, output_path: Path) -> dict:
    """
    Build cache index from cached embedding files.

    Scans results/embeddings/cache/ and generates auditable index.

    Args:
        config: config dict
        output_path: where to save cache_index.json

    Returns:
        cache index dict
    """
    cache_dir = Path(config["paths"]["embeddings"]) / "cache"

    if not cache_dir.exists():
        return {
            "model": config["embeddings"]["model"],
            "total_cached": 0,
            "entries": [],
        }

    entries = []
    for subdir in cache_dir.glob("**/"):
        for npy_file in subdir.glob("*.npy"):
            try:
                import numpy as np
                emb = np.load(npy_file)
                entries.append({
                    "cache_key": npy_file.stem,
                    "shape": list(emb.shape),
                })
            except Exception as e:
                print(f"Warning: failed to read {npy_file}: {e}")

    index = {
        "model": config["embeddings"]["model"],
        "total_cached": len(entries),
        "entries": entries,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(index, f, indent=2)

    return index
