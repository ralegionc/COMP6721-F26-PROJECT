"""Our evaluation beyond the checklist run (used by run_evaluation.py).

Detection view: "is signature k present in this tumor?" is a detection
question, so besides F1 at the 5% threshold we report threshold-free and
operating-point measures:
- ROC-AUC and average precision (area under the precision-recall curve),
  with the estimated weight as the detection score, per signature then
  averaged (see detection_report for why not pooled);
- the hit rate (recall) at a fixed 5% false-alarm rate, where the threshold
  is set on VALIDATION tumors, never on the test tumors;
- the false discovery rate at that threshold;
- a limit of detection: hit rate against how many mutations the signature
  actually contributed (true weight x mutation count).

M3 robustness sets (simulated, saved once like the main test sets):
- "unseen signature": a share of each tumor's mutations comes from a
  signature that is NOT in the cancer type's fitting list (overlaps the
  team's E3, which removes a signature from the list instead);
- "inverted prior": exposures drawn with the calibrated Dirichlet parameters
  reversed, so common signatures become rare and rare ones common. A method
  that learned the usual mix (our M3) is exposed; NNLS and EM have no prior.

Paired Wilcoxon tests with Holm's correction are not repeated here: they
come from sigdet.metrics.compare_methods, run by run_checklist.py
(results/method_tests_<source>.csv), which already covers every method.
"""
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from .config import N_GRID, PRESENCE_THRESHOLD
from .simulate import simulate

# E3: the held-out signature per cancer type (from our Week 8 plan), and the
# share of each tumor's mutations it contributes.
E3_UNSEEN = {"LUAD": "SBS17a", "SKCM": "SBS17a", "COAD": "SBS10a"}
E3_UNSEEN_SHARE = 0.2
FALSE_ALARM_RATE = 0.05
LOD_BINS = [0, 2, 4, 8, 16, 32, 64, 128, 256, np.inf]


# ---------------------------------------------------------------------------
# Detection metrics
# ---------------------------------------------------------------------------

def detection_threshold(E_val_true, E_val_pred, false_alarm=FALSE_ALARM_RATE,
                        threshold=PRESENCE_THRESHOLD):
    """The weight above which we call a signature present, chosen so that
    `false_alarm` of the truly absent (signature, tumor) pairs in the
    VALIDATION set are called present."""
    absent = np.asarray(E_val_true) < threshold
    return float(np.quantile(np.asarray(E_val_pred)[absent], 1 - false_alarm))


def detection_report(E_true, E_pred, call_threshold, threshold=PRESENCE_THRESHOLD):
    """Detection metrics for one test set.

    ROC-AUC and average precision are computed PER SIGNATURE (across tumors)
    and then averaged ("macro"). Pooling all signatures into one curve would
    reward simply knowing which signatures are common: the B0 baseline, which
    gives every tumor the same answer, would score about 0.8. Per signature,
    B0 scores 0.5, as a method with no information about the tumor should.
    The hit rate and false discovery rate use the single validation-chosen
    threshold over all (tumor, signature) pairs, i.e. the operating point.
    """
    Y = np.asarray(E_true) >= threshold
    P = np.asarray(E_pred, dtype=float)
    aucs, aps = [], []
    for k in range(Y.shape[1]):
        if 0 < Y[:, k].sum() < len(Y):  # undefined when a signature is always or never present
            aucs.append(roc_auc_score(Y[:, k], P[:, k]))
            aps.append(average_precision_score(Y[:, k], P[:, k]))
    y, s = Y.ravel(), P.ravel()
    called = s > call_threshold
    return {
        "roc_auc_macro": float(np.mean(aucs)) if aucs else np.nan,
        "avg_precision_macro": float(np.mean(aps)) if aps else np.nan,
        "hit_rate_at_fa5": float(called[y].mean()) if y.any() else np.nan,
        "false_alarm_rate": float(called[~y].mean()) if (~y).any() else np.nan,
        "fdr": float((called & ~y).sum() / max(called.sum(), 1)),
    }


def limit_of_detection(E_true, E_pred, n_mutations, call_threshold, threshold=PRESENCE_THRESHOLD):
    """Hit rate of truly present signatures, binned by how many mutations the
    signature contributed (true weight x mutation count)."""
    E_true, E_pred = np.asarray(E_true), np.asarray(E_pred)
    contributed = E_true * np.asarray(n_mutations, dtype=float)[:, None]
    present = E_true >= threshold
    rows = []
    for lo, hi in zip(LOD_BINS[:-1], LOD_BINS[1:]):
        sel = present & (contributed >= lo) & (contributed < hi)
        if sel.sum():
            rows.append({"mutations_from_signature": f"{lo:g}-{hi:g}" if np.isfinite(hi) else f"{lo:g}+",
                         "bin_low": lo, "n_pairs": int(sel.sum()),
                         "hit_rate": float((E_pred[sel] > call_threshold).mean())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# E3 robustness sets
# ---------------------------------------------------------------------------

def simulate_unseen(stats, signatures, unseen, share, n_per_point, rng, n_grid=N_GRID):
    """Test tumors where `share` of the mutations come from signature `unseen`,
    which is not in the fitting list. The truth we score against is the
    exposure over the fitting list (the remaining 1 - share, rescaled to 1),
    which is what an ideal method restricted to that list should report."""
    one = {"signatures": [unseen], "alpha": [1.0]}
    chunks = []
    for n in n_grid:
        n_extra = int(round(share * n))
        base = simulate(stats, signatures, n_per_point, rng, n_mutations=n - n_extra)
        extra = simulate(one, signatures, n_per_point, rng, n_mutations=n_extra)
        chunks.append({"counts": base["counts"] + extra["counts"], "exposures": base["exposures"],
                       "n_mutations": np.full(n_per_point, n)})
    return {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}


def inverted_prior(stats):
    """The calibrated Dirichlet parameters with their order reversed: the most
    common signature gets the smallest parameter and vice versa."""
    alpha = np.asarray(stats["alpha"], dtype=float)
    order = np.argsort(alpha)
    flipped = np.empty_like(alpha)
    flipped[order] = alpha[order[::-1]]
    return {**stats, "alpha": flipped.tolist()}


def simulate_grid(stats, signatures, n_per_point, rng, n_grid=N_GRID):
    """Test tumors at every mutation count of the grid (like the main test set)."""
    chunks = [simulate(stats, signatures, n_per_point, rng, n_mutations=n) for n in n_grid]
    return {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
