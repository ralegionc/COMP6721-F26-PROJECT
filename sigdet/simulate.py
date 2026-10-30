"""Our simulated tumors with known exposures (our only source of ground truth).

Q2 guard: calibrate() uses TRAINING patients only. Test patients never
influence what our simulated tumors look like.

We added realism knobs, so no method is tested only on its own assumptions:
- overdispersion: we draw counts with extra variability (Dirichlet-
  multinomial) instead of an idealized multinomial. Lower = noisier.
- perturb: we randomly distort each tumor's signatures around the
  reference profiles. Lower = larger distortions.
"""
import numpy as np
import pandas as pd

from .splits import patient_id

_EPS = 1e-9  # keeps Dirichlet parameters strictly positive


def calibrate(exposures, n_mutations, split):
    """We estimate mixture and mutation-count statistics from training tumors.

    exposures: DataFrame indexed by sample barcode, one column per signature
        (in practice, exposures we fitted to real TCGA catalogs).
    n_mutations: Series of total mutation counts, same index.
    split: dict of partition -> patient IDs.
    """
    train = set(split["train"])
    patients = exposures.index.map(patient_id)
    mask = np.asarray(patients.isin(train))
    if not mask.any():
        raise ValueError("No training tumors found; check the split and barcodes")
    E = exposures.to_numpy(dtype=float)[mask]
    n = np.asarray(n_mutations.loc[exposures.index], dtype=float)[mask]

    # We fit a Dirichlet by the method of moments: mean shares give the shape,
    # their spread gives the concentration (how similar tumors are to each other).
    mean = np.clip(E.mean(axis=0), 1e-4, None)
    mean = mean / mean.sum()
    var = E.var(axis=0, ddof=1) if len(E) > 1 else np.full(E.shape[1], 1e-3)
    ok = var > 1e-12
    # For a Dirichlet, var = mean(1 - mean) / (conc + 1); we solve for conc per
    # signature and take the median, which resists outliers.
    conc = np.median(mean[ok] * (1 - mean[ok]) / var[ok] - 1) if ok.any() else 10.0
    conc = float(np.clip(conc, 1.0, 1000.0))

    # We model mutation counts as log-normal, since they span several orders of magnitude.
    logn = np.log10(np.clip(n, 1, None))
    return {
        "signatures": list(exposures.columns),
        "alpha": (mean * conc).tolist(),
        "log10_n_mean": float(logn.mean()),
        "log10_n_sd": float(logn.std(ddof=1)) if len(logn) > 1 else 0.3,
        "n_train_tumors": int(mask.sum()),
        # We record who was used, so assert_calibrated_on_train() can audit it.
        "train_patients": sorted(set(patients[mask])),
    }


def assert_calibrated_on_train(stats, split):
    """We raise an error if calibration used any non-training patient (Q2)."""
    used = set(stats["train_patients"])
    outside = used - set(split["train"])
    if outside:
        raise RuntimeError(f"Simulator calibrated on non-training patients: {sorted(outside)[:3]}")
    return True


def simulate(stats, signatures, n_tumors, rng, n_mutations=None,
             overdispersion=200.0, perturb=500.0, n_bounds=(10, 20000)):
    """We draw simulated tumors.

    n_mutations: int for a fixed count (E2 grid), or None to sample counts
        from the calibrated distribution.
    Returns dict with counts (tumors x 96), exposures (tumors x K) and
    n_mutations (tumors,).
    """
    S = signatures[stats["signatures"]].to_numpy(dtype=float) + _EPS
    S = S / S.sum(axis=0, keepdims=True)
    alpha = np.asarray(stats["alpha"], dtype=float)
    # Step 1: the true exposures, one Dirichlet draw per tumor.
    E = rng.dirichlet(alpha, size=n_tumors)

    # Step 2: how many mutations each tumor has.
    if n_mutations is None:
        logn = rng.normal(stats["log10_n_mean"], stats["log10_n_sd"], size=n_tumors)
        N = np.clip(np.round(10 ** logn), *n_bounds).astype(int)
    else:
        N = np.full(n_tumors, int(n_mutations))

    # Step 3: the mutation counts over the 96 contexts.
    counts = np.empty((n_tumors, S.shape[0]), dtype=np.int64)
    for i in range(n_tumors):
        if perturb:
            # Each tumor gets its own slightly distorted copy of every signature.
            S_i = np.column_stack([rng.dirichlet(S[:, k] * perturb + _EPS) for k in range(S.shape[1])])
        else:
            S_i = S
        p = S_i @ E[i]  # the tumor's expected mutation profile
        p = p / p.sum()
        if overdispersion:
            # Extra tumor-to-tumor noise beyond plain multinomial sampling.
            p = rng.dirichlet(p * overdispersion + _EPS)
        counts[i] = rng.multinomial(N[i], p)
    return {"counts": counts, "exposures": E, "n_mutations": N}
