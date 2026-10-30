"""Q4. Would majority-class accuracy overstate success?

Most signatures are absent in most tumors, so "always absent" scores high
accuracy while finding nothing. majority_class_check() measures exactly how
much accuracy would overstate success on a given test set. The metrics we
actually report do not suffer from that inflation:

- cosine similarity between true and estimated exposure vectors (1 = same)
- MAE (mean absolute error) on signatures that are truly active
- per-signature precision, recall, F1 and PR-AUC for presence detection
- breakdown point N*
"""
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, f1_score, precision_score, recall_score

from .config import BREAKDOWN_COSINE, PRESENCE_THRESHOLD
from .labels import presence_labels


def cosine_per_tumor(E_true, E_pred):
    """We compute the cosine similarity of each tumor's true and estimated vectors."""
    a = np.asarray(E_true, dtype=float)
    b = np.asarray(E_pred, dtype=float)
    num = (a * b).sum(axis=1)
    den = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    # An all-zero vector has no direction; we score it 0 instead of dividing by zero.
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def mae_active(E_true, E_pred, threshold=PRESENCE_THRESHOLD):
    """We average the absolute error over (tumor, signature) pairs that are truly present."""
    a = np.asarray(E_true, dtype=float)
    b = np.asarray(E_pred, dtype=float)
    mask = a >= threshold
    return float(np.abs(a - b)[mask].mean()) if mask.any() else float("nan")


def majority_class_check(E_true, threshold=PRESENCE_THRESHOLD):
    """What an 'always absent' predictor would score on this test set."""
    y = presence_labels(E_true, threshold)
    absent_rate = 1.0 - y.mean()
    return {
        "absent_label_rate": float(absent_rate),
        # Predicting "absent" everywhere is right exactly as often as labels are absent.
        "always_absent_accuracy": float(absent_rate),
        "always_absent_recall": 0.0,
        "always_absent_f1": 0.0,
        "accuracy_would_overstate": bool(absent_rate > 0.5),
    }


def presence_report(E_true, E_pred, names, threshold=PRESENCE_THRESHOLD):
    """Our per-signature presence metrics. Positive class: weight >= threshold.

    PR-AUC uses the estimated weight as a score, so it does not depend on
    the threshold. It is undefined (NaN) when a signature is always or
    never present in the test set.
    """
    y = presence_labels(E_true, threshold)
    yhat = presence_labels(E_pred, threshold)
    scores = np.asarray(E_pred, dtype=float)
    rows = []
    for k, name in enumerate(names):
        support = int(y[:, k].sum())  # number of tumors where signature k is truly present
        defined = 0 < support < len(y)
        rows.append({
            "signature": name,
            "support": support,
            "precision": precision_score(y[:, k], yhat[:, k], zero_division=0),
            "recall": recall_score(y[:, k], yhat[:, k], zero_division=0),
            "f1": f1_score(y[:, k], yhat[:, k], zero_division=0),
            "pr_auc": average_precision_score(y[:, k], scores[:, k]) if defined else np.nan,
        })
    table = pd.DataFrame(rows)
    summary = {
        # Macro F1 skips signatures that never occur, since their F1 is meaningless.
        "macro_f1": float(table.loc[table["support"] > 0, "f1"].mean()),
        "micro_f1": float(f1_score(y.ravel(), yhat.ravel(), zero_division=0)),
        "macro_pr_auc": float(table["pr_auc"].mean(skipna=True)),
        "accuracy": float((y == yhat).mean()),  # we report this only to show inflation
    }
    return table, summary


def confusion_matrix_signatures(E_true, E_pred, names):
    """Our signature confusion matrix: rows = true dominant signature,
    columns = estimated dominant signature, counts of tumors."""
    t = np.asarray(E_true).argmax(axis=1)
    p = np.asarray(E_pred).argmax(axis=1)
    m = np.zeros((len(names), len(names)), dtype=int)
    np.add.at(m, (t, p), 1)  # adds 1 at each (true, estimated) pair, repeats included
    return pd.DataFrame(m, index=[f"true_{n}" for n in names], columns=[f"est_{n}" for n in names])


def breakdown_point(per_tumor, threshold=BREAKDOWN_COSINE):
    """We find the smallest N whose median cosine reaches the threshold and
    stays there for every larger N. per_tumor needs columns: method, n, cosine."""
    out = {}
    for method, g in per_tumor.groupby("method"):
        med = g.groupby("n")["cosine"].median().sort_index()
        ok = med >= threshold
        n_star = np.nan
        # We walk down from the largest N and stop at the first failure,
        # so a lucky small N after a failure cannot count.
        for n in med.index[::-1]:
            if ok[n]:
                n_star = n
            else:
                break
        out[method] = n_star
    return pd.Series(out, name="breakdown_N")
