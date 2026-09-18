"""Bounded-memory file hashing shared by ATAC production scripts."""

import hashlib


DEFAULT_CHUNK_SIZE = 1024 * 1024


def sha256_file(path, chunk_size=DEFAULT_CHUNK_SIZE):
    """Return a SHA-256 hex digest without reading the whole file into memory."""
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()
