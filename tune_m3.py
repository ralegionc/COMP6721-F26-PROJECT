"""Our M3 hyper-parameter search, on VALIDATION data only (Q2).

    python tune_m3.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt

Needs the simulated sets from run_checklist.py (data/sims/real/). For each
cancer type and setting we train M3 and record the best validation error
(mean squared error of the exposures). The test tumors are never loaded. The winning setting is the
default in sigdet/methods.py (M3MLP). Writes results/m3_tuning_real.csv.
"""
import argparse
import itertools

import pandas as pd

from sigdet.config import CANCER_SIGNATURES
from sigdet.methods import M3MLP
from sigdet.storage import ensure_dirs, load_sim_set, save_results
from sigdet.target import load_signatures, signature_subset

GRID = {"hidden": [(64, 32), (128, 64)], "alpha": [1e-4, 1e-3, 1e-2, 3e-2, 1e-1], "n_augment": [8, 16]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signatures", required=True, help="COSMIC SBS signature file")
    args = ap.parse_args()
    paths = ensure_dirs()
    signatures, _ = load_signatures(args.signatures)
    rows = []
    for cancer in CANCER_SIGNATURES:
        S = signature_subset(signatures, cancer)
        # Only train and validation sets: the test set stays untouched.
        sets = {p: load_sim_set(paths["sims"] / "real" / f"{cancer}_{p}.npz")[0] for p in ("train", "val")}
        for hidden, alpha, aug in itertools.product(*GRID.values()):
            m = M3MLP(S, hidden=hidden, alpha=alpha, n_augment=aug).fit(sets["train"], sets["val"])
            rows.append({"cancer": cancer, "hidden": str(hidden), "alpha": alpha, "n_augment": aug,
                         "best_val_loss": m.best_val_loss_, "epochs_run": len(m.history_)})
        print(f"  {cancer}: done")
    table = pd.DataFrame(rows)
    save_results(table, "m3_tuning_real", paths)
    mean = table.groupby(["hidden", "alpha", "n_augment"])["best_val_loss"].mean().sort_values()
    print("\nMean best validation error over cancer types (lower is better):\n" + mean.round(4).to_string())


if __name__ == "__main__":
    main()
