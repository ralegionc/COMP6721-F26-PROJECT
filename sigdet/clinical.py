"""Clinical data: we read smoking status and pack-years per TCGA patient.

The GDC returns each case with a list of "exposures" records. We use the
fields exposures.tobacco_smoking_status and exposures.pack_years_smoked, and
we turn values the GDC marks as missing ("--", "Unknown", "Not Reported")
into NaN.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .labels import smoking_labels

# Lower-cased strings that we treat as "no information".
_MISSING = {"", "--", "'--", "unknown", "not reported", "not allowed to collect", "none", "nan"}


def _first(values):
    """We take the first non-missing value (a case can have several exposure records)."""
    for v in values:
        if v is not None and str(v).strip().lower() not in _MISSING:
            return v
    return None


def flatten_gdc_cases(hits):
    """We turn GDC /cases JSON hits into one row per patient."""
    rows = []
    for h in hits:
        exps = h.get("exposures") or []
        rows.append({
            "patient": h.get("submitter_id"),
            "project": (h.get("project") or {}).get("project_id"),
            "smoking_raw": _first(e.get("tobacco_smoking_status") for e in exps),
            "pack_years": _first(e.get("pack_years_smoked") for e in exps),
        })
    df = pd.DataFrame(rows, columns=["patient", "project", "smoking_raw", "pack_years"])
    df["pack_years"] = pd.to_numeric(df["pack_years"], errors="coerce")
    return df.drop_duplicates("patient").set_index("patient")


def load_clinical(path):
    """We load clinical data saved by gdc_download.py (cases.json) or a CSV with
    columns patient, project, smoking_raw, pack_years, and add 'smoker' (Q3)."""
    path = Path(path)
    if path.suffix == ".json":
        df = flatten_gdc_cases(json.loads(path.read_text()))
    else:
        df = pd.read_csv(path, dtype={"smoking_raw": str}).set_index("patient")
        df["pack_years"] = pd.to_numeric(df["pack_years"], errors="coerce")
    # We pass None (not NaN) for missing values so smoking_labels treats them uniformly.
    df["smoker"] = smoking_labels(df["smoking_raw"].where(df["smoking_raw"].notna(), None)).values
    return df


def unmapped_smoking_values(df):
    """Raw values that ended up NaN although present: we check these by hand."""
    raw = df["smoking_raw"]
    mask = raw.notna() & df["smoker"].isna() & ~raw.astype(str).str.strip().str.lower().isin(_MISSING)
    return raw[mask].value_counts()
