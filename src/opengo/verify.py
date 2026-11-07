"""Input integrity check against the supplied FILE_MANIFEST.tsv (SHA-256 + size)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from opengo.config import PROJECT_ROOT

MANIFEST_PATH = PROJECT_ROOT / "FILE_MANIFEST.tsv"
CHUNK = 1 << 20  # 1 MiB

# Supplied files relocated (byte-identical) to free their path for repository files.
# The supplied README.md moved so the repository README can live at the root.
RELOCATED = {"README.md": "docs/supplied/README.md"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def verify_inputs(root: Path = PROJECT_ROOT, manifest: Path = MANIFEST_PATH) -> pd.DataFrame:
    """One row per manifest entry; status is ok | missing | size_mismatch | hash_mismatch."""
    expected = pd.read_csv(manifest, sep="\t")
    rows = []
    for rec in expected.itertuples(index=False):
        path = root / RELOCATED.get(rec.path, rec.path)
        if not path.exists():
            rows.append((rec.path, rec.bytes, None, "missing"))
            continue
        size = path.stat().st_size
        if size != rec.bytes:
            rows.append((rec.path, rec.bytes, size, "size_mismatch"))
            continue
        status = "ok" if sha256(path) == rec.sha256 else "hash_mismatch"
        rows.append((rec.path, rec.bytes, size, status))
    return pd.DataFrame(rows, columns=["path", "expected_bytes", "observed_bytes", "status"])
