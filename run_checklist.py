"""Our main experiment run: checklist Q1-Q6, baseline and stronger methods,
accuracy vs mutation count (E1, E2) and, with real data, the smoking test (E5).

    python run_checklist.py                                  # demo data
    python run_checklist.py --signatures <COSMIC SBS file>   # real run

We use real catalogs from data/processed/ when they exist for every cancer
type (made by preprocess_data.py); otherwise we build a synthetic demo cohort.
We store splits and simulated sets per data source (demo / real), save them
once, and reload them with fingerprint checks on every later run.
"""
import argparse
import json
import time

import matplotlib

matplotlib.use("Agg")  # draw to files only; we never open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sigdet.clinical import load_clinical
from sigdet.config import CANCER_SIGNATURES, CONTEXTS, N_GRID, SEED
from sigdet.labels import POSITIVE_CLASS, smoking_labels
from sigdet.methods import B0CohortMean, M0NNLS, M1ElasticNetNNLS, M2MultinomialEM
from sigdet.metrics import (breakdown_point, confusion_matrix_signatures, cosine_per_tumor,
                            mae_active, majority_class_check, presence_report)
from sigdet.realdata import negative_control, smoking_test
from sigdet.simulate import assert_calibrated_on_train, calibrate, simulate
from sigdet.splits import (assert_no_patient_leak, assign_partition, load_split,
                           make_patient_split, patient_id, save_split)
from sigdet.storage import (ensure_dirs, load_sim_set, save_figure, save_results,
                            save_sim_set, sha256_file)
from sigdet.target import TARGET_DESCRIPTION, load_signatures, signature_subset, validate_exposures

# Real TCGA tissue-source-site codes, so our demo barcodes look like real ones.
TSS = {"LUAD": "05", "SKCM": "EE", "COAD": "AA"}
# Simulated tumors per cancer type: training set, validation set, and test
# tumors at each mutation count in N_GRID.
N_TRAIN_SIM, N_VAL_SIM, N_PER_GRID_POINT = 1000, 200, 50


def make_demo_cohort(signatures, n_patients=60, seed=SEED):
    """Our synthetic stand-in for real catalogs; ~20% of patients have 2 samples."""
    rng = np.random.default_rng([seed, 999])  # a separate stream from the real-data seeds
    rows, cancers, smoking = [], [], {}
    smoking_values = ["Lifelong Non-Smoker", "Current Smoker", "Current Reformed Smoker for > 15 yrs",
                      "Current Reformed Smoker for < or = 15 yrs", "Not Reported"]
    for cancer, sigs in CANCER_SIGNATURES.items():
        S = signature_subset(signatures, cancer).to_numpy()
        for i in range(1, n_patients + 1):
            pid = f"TCGA-{TSS[cancer]}-{i:04d}"
            if cancer == "LUAD":
                smoking[pid] = rng.choice(smoking_values)
            # Every patient has a primary sample (01A); some also get a metastasis (06A),
            # which exercises our "keep a patient's samples together" logic.
            for sample in (["01A"] + (["06A"] if rng.random() < 0.2 else [])):
                e = rng.dirichlet(np.full(len(sigs), 0.5))
                n = int(np.clip(10 ** rng.normal(2.3, 0.5), 10, 20000))
                p = S @ e
                rows.append((f"{pid}-{sample}-01D-A000-08", rng.multinomial(n, p / p.sum())))
                cancers.append(cancer)
    catalogs = pd.DataFrame([r[1] for r in rows], index=[r[0] for r in rows], columns=CONTEXTS)
    clinical = pd.DataFrame({"smoking_raw": pd.Series(smoking)}).rename_axis("patient")
    clinical["pack_years"] = np.nan
    clinical["smoker"] = smoking_labels(clinical["smoking_raw"]).values
    return catalogs, pd.Series(cancers, index=catalogs.index), clinical


def load_real_cohort(paths):
    """We load our real catalogs and clinical table, or return None if any are missing."""
    files = {c: paths["processed"] / f"{c}_sbs96.csv" for c in CANCER_SIGNATURES}
    if not all(f.exists() for f in files.values()):
        return None
    parts, cancers = [], []
    for c, f in files.items():
        df = pd.read_csv(f, index_col=0)[CONTEXTS]
        df = df[df.sum(axis=1) > 0]  # a sample with no usable mutations tells us nothing
        parts.append(df)
        cancers += [c] * len(df)
    catalogs = pd.concat(parts)
    clin_path = paths["processed"] / "clinical.csv"
    clinical = load_clinical(clin_path) if clin_path.exists() else None
    return catalogs, pd.Series(cancers, index=catalogs.index), clinical


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signatures", help="COSMIC SBS signature file (tab-separated)")
    ap.add_argument("--demo", action="store_true", help="force the synthetic demo cohort")
    args = ap.parse_args()
    t0 = time.time()
    paths = ensure_dirs()
    report = {}

    # ---- Q1 ------------------------------------------------------------
    signatures, demo_sigs = load_signatures(args.signatures)
    subsets = {c: signature_subset(signatures, c) for c in CANCER_SIGNATURES}

    real = None if args.demo else load_real_cohort(paths)
    if real is not None:
        catalogs, cancer_of, clinical = real
        source = "real"
        if demo_sigs:
            print("WARNING: real catalogs with DEMO signatures. Pass --signatures <COSMIC file>; "
                  "real-data results below are not meaningful.")
    else:
        catalogs, cancer_of, clinical = make_demo_cohort(signatures)
        source = "demo"
    print(f"Data source: {source} ({len(catalogs)} samples)")
    report["Q1"] = (f"{TARGET_DESCRIPTION} Subsets: "
                    + ", ".join(f"{c}={s.shape[1]}" for c, s in subsets.items())
                    + (" (demo signatures)" if demo_sigs else ""))

    # ---- Q2 + Q5: patient split, saved once per data source -------------
    split_path = paths["splits"] / f"patient_split_{source}.json"
    if split_path.exists():
        split, doc = load_split(split_path)
        print(f"Loaded split {split_path.name} (sha256 {doc['sha256'][:12]})")
    else:
        split = make_patient_split(catalogs.index, cancer_of.values)
        save_split(split, split_path)
        print(f"Created split {split_path.name}")
    assert_no_patient_leak(split)
    partition = pd.Series(assign_partition(catalogs.index, split), index=catalogs.index)
    # Double check: every patient's samples must all sit in a single partition.
    assert (partition.groupby(catalogs.index.map(patient_id)).nunique() == 1).all()

    # ---- per cancer: calibrate on train, simulate, evaluate --------------
    sim_dir = paths["sims"] / source
    per_tumor, summary, presence, majority, fingerprints, real_fits = [], [], [], [], {}, {}
    m1_alphas = {}
    for ci, (cancer, S_df) in enumerate(subsets.items()):
        names, S = list(S_df.columns), S_df.to_numpy()
        mask = cancer_of == cancer
        # We fit NNLS to the real tumors; these fits feed the simulator and E5.
        fitted = pd.DataFrame(M0NNLS(S).predict(catalogs[mask].to_numpy()),
                              index=catalogs.index[mask], columns=names)
        real_fits[cancer] = fitted
        # calibrate() keeps only training patients (Q2), and we verify that it did.
        stats = calibrate(fitted, catalogs[mask].sum(axis=1), split)
        assert_calibrated_on_train(stats, split)

        sets = {}
        for pi, part in enumerate(["train", "val", "test"]):
            path = sim_dir / f"{cancer}_{part}.npz"
            # We simulate each set only once; later runs reload the saved file (Q5).
            if not path.exists():
                rng = np.random.default_rng([SEED, ci, pi])  # own seed per cancer type and partition
                if part == "test":
                    # The test set has a fixed number of tumors at each mutation count (E2).
                    chunks = [simulate(stats, signatures, N_PER_GRID_POINT, rng, n_mutations=n) for n in N_GRID]
                    sim = {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
                else:
                    sim = simulate(stats, signatures, N_TRAIN_SIM if part == "train" else N_VAL_SIM, rng)
                save_sim_set(path, {"cancer": cancer, "part": part, "signatures": names,
                                    "seed": [SEED, ci, pi]}, **sim)
            # Loading checks the checksum, so a changed file fails here.
            data, _ = load_sim_set(path)
            validate_exposures(data["exposures"], names)
            sets[part] = data
            fingerprints[f"{source}/{path.name}"] = sha256_file(path)

        test = sets["test"]
        methods = [B0CohortMean(), M0NNLS(S), M1ElasticNetNNLS(S), M2MultinomialEM(S)]
        majority.append({"cancer": cancer, **majority_class_check(test["exposures"])})
        # Every method sees exactly the same train, validation and test files.
        for m in methods:
            m.fit(sets["train"], val=sets["val"])
            if isinstance(m, M1ElasticNetNNLS):
                m1_alphas[cancer] = m.alpha_
            pred = m.predict(test["counts"])
            validate_exposures(pred, names)  # each method's output must meet the Q1 definition
            cos = cosine_per_tumor(test["exposures"], pred)
            per_tumor.append(pd.DataFrame({"cancer": cancer, "method": m.name,
                                           "n": test["n_mutations"], "cosine": cos}))
            table, summ = presence_report(test["exposures"], pred, names)
            presence.append(table.assign(cancer=cancer, method=m.name))
            summary.append({"cancer": cancer, "method": m.name, "median_cosine": float(np.median(cos)),
                            "mae_active": mae_active(test["exposures"], pred), **summ})
            if m.name in ("M0_nnls", "M2_multinomial_em"):
                save_results(confusion_matrix_signatures(test["exposures"], pred, names),
                             f"confusion_{cancer}_{m.name}_{source}", paths, index=True)
        print(f"  {cancer}: done ({time.time() - t0:.0f}s)")

    (paths["splits"] / f"sim_fingerprints_{source}.json").write_text(json.dumps(fingerprints, indent=2))

    # ---- gather results across cancer types and save them ----------------
    per_tumor = pd.concat(per_tumor, ignore_index=True)
    summary = pd.DataFrame(summary)
    majority = pd.DataFrame(majority)
    breakdown = pd.concat([breakdown_point(g).rename(c) for c, g in per_tumor.groupby("cancer")], axis=1)
    save_results(per_tumor, f"per_tumor_cosine_{source}", paths)
    save_results(summary, f"summary_{source}", paths)
    save_results(pd.concat(presence, ignore_index=True), f"presence_per_signature_{source}", paths)
    save_results(majority, f"majority_class_check_{source}", paths)
    save_results(breakdown, f"breakdown_point_{source}", paths, index=True)

    # One panel per cancer type: median cosine vs mutation count, with the 0.90 line for N*.
    fig, axes = plt.subplots(1, len(subsets), figsize=(4.2 * len(subsets), 3.4), sharey=True)
    for ax, (cancer, g) in zip(np.atleast_1d(axes), per_tumor.groupby("cancer")):
        for method, h in g.groupby("method"):
            med = h.groupby("n")["cosine"].median()
            ax.plot(med.index, med.values, marker="o", ms=3, label=method)
        ax.axhline(0.9, color="grey", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_title(cancer)
        ax.set_xlabel("mutations per tumor")
    np.atleast_1d(axes)[0].set_ylabel("median cosine (true vs estimated)")
    np.atleast_1d(axes)[0].legend(fontsize=7)
    save_figure(fig, f"cosine_vs_mutation_count_{source}", paths)

    # ---- E5: real-data smoking test + negative controls -------------------
    e5 = None
    if clinical is not None and "SBS4" in subsets["LUAD"].columns:
        e5 = smoking_test(real_fits["LUAD"], clinical)
        controls = []
        # Negative controls: we add a signature that should be absent and see how much weight it gets.
        for cancer, extra, label in [("SKCM", ["SBS4"], "SBS4 in SKCM"),
                                     ("LUAD", ["SBS7a", "SBS7b", "SBS7c", "SBS7d"], "SBS7 in LUAD")]:
            # The usual subset plus the control signatures, without duplicates.
            cols = list(dict.fromkeys(list(subsets[cancer].columns) + [s for s in extra if s in signatures.columns]))
            mask = cancer_of == cancer
            fit = pd.DataFrame(M0NNLS(signatures[cols].to_numpy()).predict(catalogs[mask].to_numpy()),
                               index=catalogs.index[mask], columns=cols)
            controls.append(negative_control(fit, set(extra), label))
        save_results(pd.DataFrame([e5]), f"e5_smoking_{source}", paths)
        save_results(pd.DataFrame(controls), f"e5_negative_controls_{source}", paths)

    # ---- checklist report ------------------------------------------------
    # For Q4 we quote the cancer type where "always absent" would look best.
    worst = majority.loc[majority["absent_label_rate"].idxmax()]
    verdict = ("accuracy WOULD overstate success" if worst["accuracy_would_overstate"]
               else "absent labels are not the majority in this run")
    m0 = summary[summary["method"] == "M0_nnls"]
    report["Q2"] = (f"Split by patient ({', '.join(f'{k}={len(v)}' for k, v in split.items())}); no overlap; "
                    "multi-sample patients kept together; simulator calibrated on training patients only.")
    sm = clinical["smoker"] if clinical is not None else pd.Series(dtype=float)
    if clinical is not None and "project" in clinical.columns:
        sm = sm[clinical["project"] == "TCGA-LUAD"]  # the smoking test uses LUAD only
    report["Q3"] = (f"Presence: {POSITIVE_CLASS['presence']}. Smoking: {POSITIVE_CLASS['smoking']}. "
                    f"LUAD labels: {int((sm == 1).sum())} positive, {int((sm == 0).sum())} negative, "
                    f"{int(sm.isna().sum())} excluded.")
    report["Q4"] = (f"Max 'absent' rate {worst['absent_label_rate']:.0%} ({worst['cancer']}); {verdict}. "
                    f"Reported: per-signature F1/PR-AUC, cosine, MAE vs B0 (M0 macro F1 "
                    f"{m0['macro_f1'].min():.2f}-{m0['macro_f1'].max():.2f}).")
    report["Q5"] = (f"Split + {len(fingerprints)} simulated sets saved once (seed {SEED}, SHA-256); "
                    f"all {len(methods)} methods scored on the same files.")
    report["Q6"] = "results/, figures/, splits/ committed; data/raw and data/sims gitignored (manifest + seeds)."

    # We write the same numbers as a Markdown summary for the progress report.
    lines = [f"# Run summary ({source} data, {time.strftime('%Y-%m-%d')})", "",
             "## Checklist"] + [f"- {q}: {report[q]}" for q in sorted(report)] + [
             "", "## Test-set results (simulated tumors, all mutation counts pooled)", "",
             "| cancer | method | median cosine | MAE (active) | macro F1 | PR-AUC |",
             "| --- | --- | --- | --- | --- | --- |"]
    for r in summary.itertuples():
        lines.append(f"| {r.cancer} | {r.method} | {r.median_cosine:.3f} | {r.mae_active:.3f} | "
                     f"{r.macro_f1:.3f} | {r.macro_pr_auc:.3f} |")
    lines += ["", "## Breakdown point N* (smallest N with median cosine >= 0.90)", "",
              "| method | " + " | ".join(breakdown.columns) + " |",
              "| --- " * (len(breakdown.columns) + 1) + "|"]
    for method, row in breakdown.iterrows():
        lines.append(f"| {method} | " + " | ".join("not reached" if pd.isna(v) else f"{int(v)}" for v in row) + " |")
    lines += ["", f"M1 penalty chosen on validation data: {m1_alphas}"]
    if e5:
        lines += ["", "## E5 smoking test (LUAD, real data)", "", json.dumps(e5, indent=1)]
    (paths["results"] / f"run_summary_{source}.md").write_text("\n".join(lines) + "\n")

    print("\nCHECKLIST")
    for q in sorted(report):
        print(f"  {q}: {report[q]}")
    print("\n" + summary[["cancer", "method", "median_cosine", "mae_active", "macro_f1", "macro_pr_auc"]]
          .round(3).to_string(index=False))
    print("\nBreakdown point N*:\n" + breakdown.to_string())
    if e5:
        print("\nE5:", e5)
    print(f"\nWrote results/run_summary_{source}.md ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
