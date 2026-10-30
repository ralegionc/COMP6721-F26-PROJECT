"""Q1. Is the target clearly defined?

Our target: for each tumor, an exposure vector. This is one number per
candidate signature of that tumor's cancer type, giving the share of the
tumor's mutations we attribute to that signature. The numbers are
non-negative and sum to 1. We know the true values only for simulated tumors.
"""
import numpy as np
import pandas as pd

from .config import ARTIFACT_SIGNATURES, CANCER_SIGNATURES, CONTEXTS, SEED

TARGET_DESCRIPTION = (
    "Exposure vector per tumor: non-negative shares over the cancer type's "
    "candidate COSMIC SBS signatures (artifacts excluded), summing to 1. "
    "Known only for simulated tumors; real tumors give directional checks only."
)


def load_signatures(path=None):
    """We load the signature matrix (96 contexts x signatures).

    With a path, we read a COSMIC SBS tab-separated file whose first column
    holds the context labels, for example "A[C>A]A". Without a path, we return
    synthetic demo signatures so the pipeline runs before real data exists.

    Returns (DataFrame indexed by context, is_demo flag).
    """
    if path is not None:
        df = pd.read_csv(path, sep="\t")
        df = df.set_index(df.columns[0])
        missing = set(CONTEXTS) - set(df.index)
        if missing:
            raise ValueError(f"Signature file lacks {len(missing)} of the 96 contexts")
        # We reorder rows to our standard context order so every matrix lines up.
        return df.loc[CONTEXTS].astype(float), False
    return demo_signatures(), True


def demo_signatures(seed=SEED):
    """Random but reproducible stand-in signatures. Not biologically real."""
    rng = np.random.default_rng(seed)
    names = sorted(set().union(*CANCER_SIGNATURES.values()))
    cols = {}
    for name in names:
        # We make clock-like and flat signatures broad and the rest peaky,
        # loosely imitating the shapes of the real ones.
        alpha = 5.0 if name in {"SBS1", "SBS5", "SBS40a"} else 0.3
        cols[name] = rng.dirichlet(np.full(len(CONTEXTS), alpha))
    return pd.DataFrame(cols, index=CONTEXTS)


def signature_subset(signatures, cancer):
    """We pick the candidate signatures for one cancer type, with artifacts removed."""
    if cancer not in CANCER_SIGNATURES:
        raise KeyError(f"No signature list configured for cancer type {cancer!r}")
    wanted = [s for s in CANCER_SIGNATURES[cancer] if s not in ARTIFACT_SIGNATURES]
    missing = [s for s in wanted if s not in signatures.columns]
    if missing:
        raise ValueError(f"{cancer}: signatures missing from matrix: {missing}")
    return signatures[wanted]


def validate_exposures(exposures, names=None, atol=1e-6):
    """We check that an exposure table matches our target definition.

    We raise ValueError on wrong columns, negative or non-finite values, or
    rows that do not sum to 1, and return True otherwise.
    """
    if isinstance(exposures, pd.DataFrame):
        if names is not None and list(exposures.columns) != list(names):
            raise ValueError("Exposure columns do not match the signature subset")
        values = exposures.to_numpy(dtype=float)
    else:
        values = np.asarray(exposures, dtype=float)
        if names is not None and values.shape[1] != len(names):
            raise ValueError("Exposure width does not match the signature subset")
    if values.ndim != 2:
        raise ValueError("Exposures must be a 2-D table: tumors x signatures")
    if not np.isfinite(values).all():
        raise ValueError("Exposures contain NaN or infinite values")
    # atol allows for tiny negative rounding errors from the solvers.
    if (values < -atol).any():
        raise ValueError("Exposures contain negative values")
    sums = values.sum(axis=1)
    bad = np.flatnonzero(np.abs(sums - 1.0) > atol)
    if bad.size:
        raise ValueError(f"{bad.size} tumors have exposures not summing to 1 (first row {bad[0]})")
    return True
