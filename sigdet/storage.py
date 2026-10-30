"""Q6. Where will figures and intermediate results be stored?

    data/raw/        downloaded TCGA/PCAWG files      NOT committed (gitignored)
    data/processed/  96-context catalogs              committed if small
    data/sims/       simulated train/val/test sets    NOT committed; rebuilt from seeds
    splits/          saved patient splits             committed
    results/         metrics as CSV                   committed
    figures/         plots made by scripts            committed

In Git we replace the raw data with data/manifest.csv: every file's size and
SHA-256 checksum, so anyone can re-download and confirm identical files.
Our simulated sets carry a .sha256 sidecar so we detect a changed test set (Q5).
"""
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .config import ROOT

# Written only if the project has no .gitignore yet.
GITIGNORE = """# Raw and large generated data (rebuilt from manifest and seeds)
data/raw/
data/sims/
__pycache__/
*.pyc
.pytest_cache/
"""

# GDC file UUIDs, which we pick out of download paths for the manifest.
_UUID_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")


def project_paths(root=ROOT):
    """Every folder our scripts read from or write to, in one place."""
    root = Path(root)
    return {
        "raw": root / "data" / "raw",
        "processed": root / "data" / "processed",
        "sims": root / "data" / "sims",
        "manifest": root / "data" / "manifest.csv",
        "splits": root / "splits",
        "results": root / "results",
        "figures": root / "figures",
    }


def ensure_dirs(root=ROOT):
    """We create any missing folders (and a default .gitignore) and return the paths."""
    paths = project_paths(root)
    for key, p in paths.items():
        if key != "manifest":  # the manifest is a file, not a folder
            p.mkdir(parents=True, exist_ok=True)
    gi = Path(root) / ".gitignore"
    if not gi.exists():
        gi.write_text(GITIGNORE)
    return paths


def sha256_file(path, chunk=1 << 20):
    """We hash a file in 1 MB chunks, so large files never load fully into memory."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def write_manifest(raw_dir, out_csv):
    """We list every raw file with size, checksum and GDC file UUID (if the path has one)."""
    raw_dir = Path(raw_dir)
    rows = []
    for p in sorted(raw_dir.rglob("*")):
        if p.is_file():
            rel = p.relative_to(raw_dir).as_posix()
            m = _UUID_RE.search(rel.lower())
            rows.append({
                "path": rel,
                "bytes": p.stat().st_size,
                "sha256": sha256_file(p),
                "gdc_file_uuid": m.group(0) if m else "",
            })
    df = pd.DataFrame(rows, columns=["path", "bytes", "sha256", "gdc_file_uuid"])
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def verify_manifest(manifest_csv, raw_dir):
    """We return a list of problems (missing or changed files). Empty = all good."""
    raw_dir = Path(raw_dir)
    problems = []
    for row in pd.read_csv(manifest_csv, dtype=str, keep_default_na=False).itertuples():
        p = raw_dir / row.path
        if not p.exists():
            problems.append(f"missing: {row.path}")
        elif sha256_file(p) != row.sha256:
            problems.append(f"changed: {row.path}")
    return problems


def save_sim_set(path, meta, overwrite=False, **arrays):
    """We save a simulated set once, with a checksum sidecar."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists. Load it with load_sim_set; do not regenerate.")
    path.parent.mkdir(parents=True, exist_ok=True)
    # We store metadata as a JSON string so loading never needs pickle.
    np.savez_compressed(path, meta=np.array(json.dumps(meta)), **arrays)
    digest = sha256_file(path)
    path.with_suffix(".sha256").write_text(digest + "\n")
    return digest


def load_sim_set(path):
    """We load a simulated set and verify it matches its checksum."""
    path = Path(path)
    expected = path.with_suffix(".sha256").read_text().strip()
    if sha256_file(path) != expected:
        raise ValueError(f"{path} changed since it was saved (checksum mismatch)")
    with np.load(path, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files if k != "meta"}
        meta = json.loads(str(z["meta"]))
    return data, meta


def save_results(df, name, paths, index=False):
    """We write a results table to results/<name>.csv."""
    out = paths["results"] / f"{name}.csv"
    df.to_csv(out, index=index)
    return out


def save_figure(fig, name, paths, dpi=150):
    """We write a matplotlib figure to figures/<name>.png."""
    out = paths["figures"] / f"{name}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    return out
