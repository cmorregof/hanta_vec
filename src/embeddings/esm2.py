"""
ESM-2 embedding computation with caching and mean pooling.

Handles:
- Model loading (MPS/CUDA/CPU device selection)
- Mean pooling over residues (excluding special tokens)
- Embedding caching with SHA256 keys
"""

import hashlib
import numpy as np
import torch
from pathlib import Path
from transformers import EsmTokenizer, EsmModel


def load_model(config: dict) -> tuple:
    """
    Load ESM-2 model and tokenizer with device selection.

    Priority: MPS (Apple Silicon) > CUDA > CPU

    Returns:
        (model, tokenizer, device)
    """
    model_name = config["embeddings"]["model"]

    tokenizer = EsmTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)

    # Device selection
    if torch.backends.mps.is_available():
        device = torch.device("mps")
        print("✓ Using MPS (Apple Silicon)")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"✓ Using CUDA (GPU: {torch.cuda.get_device_name(0)})")
    else:
        device = torch.device("cpu")
        print("✓ Using CPU")

    model = model.to(device)
    model.eval()

    return model, tokenizer, device


def mean_pool(outputs, attention_mask) -> np.ndarray:
    """
    Mean pooling over residue tokens, excluding special tokens.

    Excludes [BOS] (position 0) and [EOS] (last position) from pooling.

    Args:
        outputs: model output with last_hidden_state
        attention_mask: (batch, seq_len) attention mask

    Returns:
        (batch, hidden_dim) numpy array on CPU
    """
    token_embeddings = outputs.last_hidden_state  # (batch, seq_len, hidden)

    # Exclude first ([BOS]) and last ([EOS]) positions
    embeddings = token_embeddings[:, 1:-1, :]  # (batch, seq_len-2, hidden)
    mask = attention_mask[:, 1:-1].unsqueeze(-1).float()  # (batch, seq_len-2, 1)

    # Masked mean
    pooled = (embeddings * mask).sum(1) / mask.sum(1).clamp(min=1e-9)

    return pooled.cpu().numpy()


def get_or_compute_embedding(
    seq: str,
    accession: str,
    model,
    tokenizer,
    device,
    config: dict,
) -> np.ndarray:
    """
    Get embedding from cache or compute and cache it.

    Cache key: SHA256(seq | model_name | "mean_pool")

    Args:
        seq: protein sequence
        accession: sequence accession (for logging only)
        model: ESM-2 model
        tokenizer: ESM tokenizer
        device: torch device
        config: config dict with paths and embeddings settings

    Returns:
        (hidden_dim,) numpy array
    """
    # Cache key
    cache_key = hashlib.sha256(
        f"{seq}|{model.config.name_or_path}|mean_pool".encode()
    ).hexdigest()

    cache_dir = Path(config["paths"]["embeddings"]) / "cache"
    cache_path = cache_dir / cache_key[:2] / f"{cache_key}.npy"

    # Try cache
    if cache_path.exists():
        return np.load(cache_path)

    # Truncate if needed
    max_len = config["embeddings"]["max_length"]  # 1022
    if len(seq) > max_len:
        seq = seq[:max_len]

    # Tokenize
    inputs = tokenizer(
        seq,
        return_tensors="pt",
        padding=False,
        truncation=True,
        max_length=max_len + 2,  # +2 for special tokens
    )
    inputs = {k: v.to(device) for k, v in inputs.items()}

    # Forward pass
    with torch.no_grad():
        outputs = model(**inputs)

    # Pool
    embedding = mean_pool(outputs, inputs["attention_mask"])[0]  # (hidden_dim,)

    # Cache
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache_path, embedding)

    return embedding
