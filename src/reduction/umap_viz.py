"""
UMAP dimensionality reduction for ESM-2 embeddings.
"""

import numpy as np


def run_umap(coords_50d: np.ndarray, config: dict) -> np.ndarray:
    """
    Run UMAP on first 50 PCs.

    Args:
        coords_50d: (N, 50) PCA coordinates
        config: config dict with reduction parameters

    Returns:
        coords_2d: (N, 2) UMAP coordinates
    """
    import umap

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=config["reduction"]["umap_n_neighbors"],
        min_dist=config["reduction"]["umap_min_dist"],
        metric="cosine",
        random_state=config["reduction"]["random_state"],
        low_memory=False,
    )
    return reducer.fit_transform(coords_50d)
