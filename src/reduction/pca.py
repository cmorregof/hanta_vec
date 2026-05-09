"""
PCA reduction for ESM-2 embeddings.
"""

import numpy as np
from sklearn.decomposition import PCA


def run_pca(embeddings: np.ndarray, config: dict) -> tuple:
    """
    Run PCA on embeddings to 2D and 50D.

    Args:
        embeddings: (N, 480) array of ESM-2 embeddings
        config: config dict

    Returns:
        coords_2d: (N, 2) PCA coordinates
        coords_50d: (N, 50) PCA coordinates
        explained_variance_ratio: (50,) variance explained per component
        pca_object: fitted PCA object for reference
    """
    n_components_full = min(50, embeddings.shape[0] - 1, embeddings.shape[1])

    pca_full = PCA(n_components=n_components_full, random_state=42)
    coords_50d = pca_full.fit_transform(embeddings)

    pca_2d = PCA(n_components=2, random_state=42)
    coords_2d = pca_2d.fit_transform(embeddings)

    return coords_2d, coords_50d, pca_full.explained_variance_ratio_, pca_full
