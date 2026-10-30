"""Our real-data checks (directional evidence only; real tumors have no ground truth).

E5: in lung adenocarcinoma, is our estimated tobacco-signature (SBS4) weight
higher in smokers than in lifelong non-smokers? We use
- a one-sided Mann-Whitney U test (rank-based, no normality assumption)
- the Spearman correlation of SBS4 weight with pack-years, among smokers
Our negative controls should come out near zero if attribution is sound:
- SBS4 weight in melanoma (SKCM); UV signatures (SBS7*) in LUAD
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from .config import PRESENCE_THRESHOLD
from .splits import patient_id


def per_patient(exposures):
    """We average exposures over a patient's samples (index: sample barcodes),
    so a patient with two tumor samples counts once in each test."""
    return exposures.groupby(exposures.index.map(patient_id)).mean()


def smoking_test(luad_exposures, clinical, signature="SBS4"):
    """We compare SBS4 weight between smokers and lifelong non-smokers (E5)."""
    e = per_patient(luad_exposures)[signature]
    # Only patients with both mutation data and a known smoking label take part.
    df = clinical.join(e.rename("weight"), how="inner").dropna(subset=["smoker"])
    smokers, never = df.loc[df["smoker"] == 1, "weight"], df.loc[df["smoker"] == 0, "weight"]
    out = {"n_smokers": int(len(smokers)), "n_never": int(len(never)),
           "median_smokers": float(smokers.median()) if len(smokers) else np.nan,
           "median_never": float(never.median()) if len(never) else np.nan}
    if len(smokers) and len(never):
        # One-sided: our hypothesis is that smokers have MORE SBS4, not just a difference.
        res = mannwhitneyu(smokers, never, alternative="greater")
        out.update({"mannwhitney_U": float(res.statistic), "mannwhitney_p": float(res.pvalue)})
    # Dose-response check among smokers with recorded pack-years.
    py = df[(df["smoker"] == 1) & df["pack_years"].notna()]
    if len(py) >= 5:
        rho, p = spearmanr(py["weight"], py["pack_years"])
        out.update({"spearman_rho_pack_years": float(rho), "spearman_p": float(p), "n_pack_years": int(len(py))})
    return out


def negative_control(exposures, signatures, label):
    """We measure how much weight goes to signatures that should be absent."""
    cols = [c for c in exposures.columns if c in signatures]
    if not cols:
        return {"control": label, "note": "signature not in this cancer's subset (weight fixed at 0)"}
    w = per_patient(exposures)[cols].sum(axis=1)
    return {"control": label, "median_weight": float(w.median()),
            "share_present": float((w >= PRESENCE_THRESHOLD).mean()), "n": int(len(w))}
