"""Q2. Can identifiers leak across partitions?  Q5. Can every method reuse
exactly the same split?

Q2: we split by PATIENT, not by sample. A TCGA sample barcode looks like
TCGA-05-4244-01A-01D-1105-08; its first 12 characters (TCGA-05-4244)
identify the patient. All samples of one patient land in the same partition.
We remove PCAWG donors that are also TCGA patients from PCAWG.

Q5: we write the split to disk once, with its seed and a SHA-256
fingerprint (a checksum that changes if even one ID changes). Every method
loads the same file, and loading fails if the file was edited.
"""
import hashlib
import json
import re
from datetime import date
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from .config import SEED, SPLIT_FRACTIONS

# "TCGA-" + 2-character tissue source site + "-" + 4-character participant code.
_PATIENT_RE = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})")


class LeakageError(RuntimeError):
    """Raised when one patient appears in more than one partition."""


def patient_id(barcode):
    """We turn a TCGA sample or aliquot barcode into the 12-character patient ID."""
    m = _PATIENT_RE.match(str(barcode).strip().upper())
    if not m:
        raise ValueError(f"Not a TCGA barcode: {barcode!r}")
    return m.group(1)


def make_patient_split(barcodes, cancer_types, fractions=SPLIT_FRACTIONS, seed=SEED):
    """We split patients into train/val/test, stratified by cancer type.

    Stratified means we split each cancer type separately with the same
    fractions, so every partition contains every cancer type.
    """
    if abs(sum(fractions.values()) - 1.0) > 1e-9:
        raise ValueError("Split fractions must sum to 1")
    df = pd.DataFrame({"barcode": list(barcodes), "cancer": list(cancer_types)})
    df["patient"] = df["barcode"].map(patient_id)
    # A patient listed under two cancer types would make stratification ambiguous.
    n_types = df.groupby("patient")["cancer"].nunique()
    if (n_types > 1).any():
        raise ValueError(f"Patients with several cancer types: {list(n_types[n_types > 1].index[:5])}")

    # We sort before shuffling so the result depends only on the seed, not on input order.
    patients = df.drop_duplicates("patient").sort_values("patient")
    rng = np.random.default_rng(seed)
    names = list(fractions)
    split = {name: [] for name in names}
    for _, group in patients.groupby("cancer", sort=True):
        ids = group["patient"].to_numpy().copy()
        rng.shuffle(ids)
        n = len(ids)
        n_train = int(round(n * fractions[names[0]]))
        n_val = int(round(n * fractions[names[1]]))
        # Test gets whatever remains, so rounding never loses a patient.
        split[names[0]] += ids[:n_train].tolist()
        split[names[1]] += ids[n_train:n_train + n_val].tolist()
        split[names[2]] += ids[n_train + n_val:].tolist()
    split = {k: sorted(v) for k, v in split.items()}
    assert_no_patient_leak(split)
    return split


def assert_no_patient_leak(split):
    """We raise LeakageError if any patient is in two partitions."""
    for a, b in combinations(split, 2):
        shared = set(split[a]) & set(split[b])
        if shared:
            raise LeakageError(f"{len(shared)} patients in both {a} and {b}, e.g. {sorted(shared)[:3]}")
    return True


def assign_partition(barcodes, split):
    """We look up the partition of each sample barcode through its patient ID."""
    lookup = {p: name for name, ids in split.items() for p in ids}
    out = []
    for bc in barcodes:
        pid = patient_id(bc)
        if pid not in lookup:
            raise KeyError(f"Patient {pid} is not in the saved split")
        out.append(lookup[pid])
    return np.array(out)


def remove_overlapping_donors(pcawg_donor_ids, tcga_patient_ids):
    """We return a boolean mask: True for PCAWG donors to KEEP (not also TCGA patients).

    TCGA-derived PCAWG donors carry TCGA barcodes, so we match them on the
    12-character patient ID.
    """
    tcga = {patient_id(p) for p in tcga_patient_ids}
    keep = []
    for d in pcawg_donor_ids:
        try:
            keep.append(patient_id(d) not in tcga)
        except ValueError:
            keep.append(True)  # non-TCGA donor ID, cannot overlap
    return np.array(keep, dtype=bool)


def split_fingerprint(split):
    """We hash the sorted patient lists, so any added, removed or moved ID changes the hash."""
    payload = json.dumps({k: sorted(split[k]) for k in sorted(split)}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def save_split(split, path, seed=SEED, fractions=SPLIT_FRACTIONS, overwrite=False):
    """We write the split once and refuse to overwrite it unless overwrite=True."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Load it with load_split; do not regenerate.")
    assert_no_patient_leak(split)
    doc = {
        "seed": seed,
        "fractions": fractions,
        "created": date.today().isoformat(),
        "n_patients": {k: len(v) for k, v in split.items()},
        "sha256": split_fingerprint(split),
        "partitions": {k: sorted(v) for k, v in split.items()},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2))
    return doc["sha256"]


def load_split(path):
    """We load a saved split and verify its fingerprint and leak-freedom."""
    doc = json.loads(Path(path).read_text())
    split = doc["partitions"]
    if split_fingerprint(split) != doc["sha256"]:
        raise ValueError(f"{path} was modified after it was saved (fingerprint mismatch)")
    assert_no_patient_leak(split)
    return split, doc
