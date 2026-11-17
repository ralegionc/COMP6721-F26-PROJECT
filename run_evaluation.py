"""Our evaluation run for M3 and the final figures (Bao's part).

    python run_checklist.py  --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt  # first
    python run_evaluation.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt

run_checklist.py must run first: we reuse its saved split and simulated sets
(data/sims/real/), so every number here comes from the same tumors. We only
WRITE new files (names below); we never touch splits/ or the checklist outputs.

  results/eval_detection_real.csv        per-signature ROC-AUC and average precision, hit rate at 5% false alarms, FDR
  results/eval_limit_of_detection_real.csv  hit rate vs mutations contributed by the signature
  results/eval_robustness_real.csv       M3 checks: unseen signature and inverted prior
  results/eval_e5_methods_real.csv       E5 smoking test for every method (all / test patients)
  results/eval_real_agreement_real.csv   real tumors: reconstruction and agreement with M2
  results/m3_training_history_real.csv   M3 loss curves
  results/e3_fingerprints_real.json      checksums of the robustness simulated sets
  results/evaluation_summary_real.md     the headline tables
  figures/eval_*.png                     final figures (listed in evaluation_summary_real.md)

The paired Wilcoxon tests (Holm-corrected) are not redone here: run_checklist.py already runs
sigdet.metrics.compare_methods on every pair of methods, M3 included, and writes
results/method_tests_real.csv. The summary quotes the rows that involve M3.
"""
import argparse
import json
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_checklist import load_real_cohort
from sigdet.config import CANCER_SIGNATURES, N_GRID, PRESENCE_THRESHOLD, SEED
from sigdet.evaluation import (E3_UNSEEN, E3_UNSEEN_SHARE, detection_report, detection_threshold,
                               inverted_prior, limit_of_detection, simulate_grid, simulate_unseen)
from sigdet.methods import B0CohortMean, M0NNLS, M1ElasticNetNNLS, M2MultinomialEM, M3MLP
from sigdet.metrics import cosine_per_tumor, mae_active, presence_report
from sigdet.realdata import smoking_test
from sigdet.simulate import calibrate
from sigdet.splits import assign_partition, load_split
from sigdet.storage import ensure_dirs, load_sim_set, save_figure, save_results, save_sim_set, sha256_file
from sigdet.target import load_signatures, signature_subset, validate_exposures

N_PER_GRID_POINT = 50  # same as the main test set
LABELS = {"B0_cohort_mean": "B0 cohort mean", "M0_nnls": "M0 NNLS", "M1_elastic_net": "M1 elastic-net NNLS",
          "M2_multinomial_em": "M2 multinomial EM", "M3_mlp": "M3 MLP", "M3_mlp_broad": "M3 MLP (broad)"}
COLORS = {"B0_cohort_mean": "#9e9e9e", "M0_nnls": "#222222", "M1_elastic_net": "#2ca02c",
          "M2_multinomial_em": "#1f77b4", "M3_mlp": "#d62728", "M3_mlp_broad": "#ff7f0e"}


def make_methods(S_df):
    S = S_df.to_numpy()
    return [B0CohortMean(), M0NNLS(S), M1ElasticNetNNLS(S), M2MultinomialEM(S),
            M3MLP(S_df), M3MLP(S_df, broad=True)]


def load_or_make(path, make, meta, fingerprints):
    """Our E3 sets are simulated once and reloaded with a checksum check, like the main sets."""
    if not path.exists():
        save_sim_set(path, meta, **make())
    data, _ = load_sim_set(path)
    fingerprints[f"real/{path.name}"] = sha256_file(path)
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signatures", required=True, help="COSMIC SBS signature file (tab-separated)")
    args = ap.parse_args()
    t0 = time.time()
    paths = ensure_dirs()
    signatures, demo = load_signatures(args.signatures)
    if demo:
        raise SystemExit("Pass the real COSMIC file with --signatures.")
    real = load_real_cohort(paths)
    if real is None:
        raise SystemExit("Real catalogs not found in data/processed/.")
    catalogs, cancer_of, clinical = real
    split, _ = load_split(paths["splits"] / "patient_split_real.json")
    partition = pd.Series(assign_partition(catalogs.index, split), index=catalogs.index)
    sim_dir = paths["sims"] / "real"

    per_tumor, detection, lod, robust, history, agreement, e5_rows = [], [], [], [], [], [], []
    confusion, recall_rows, fingerprints, real_fits = {}, [], {}, {}
    for ci, (cancer, names) in enumerate(CANCER_SIGNATURES.items()):
        S_df = signature_subset(signatures, cancer)
        sets = {}
        for part in ("train", "val", "test"):
            path = sim_dir / f"{cancer}_{part}.npz"
            if not path.exists():
                raise SystemExit(f"{path} is missing: run run_checklist.py --signatures ... first.")
            sets[part], _ = load_sim_set(path)

        # The simulator calibration, exactly as run_checklist.py does it (training patients only).
        mask = cancer_of == cancer
        fitted = pd.DataFrame(M0NNLS(S_df.to_numpy()).predict(catalogs[mask].to_numpy()),
                              index=catalogs.index[mask], columns=names)
        stats = calibrate(fitted, catalogs[mask].sum(axis=1), split)

        # E3 sets: own seed streams, separate from the main sets (which use [SEED, ci, 0..2]).
        e3 = {
            "unseen signature": load_or_make(
                sim_dir / f"e3_{cancer}_unseen.npz",
                lambda: simulate_unseen(stats, signatures, E3_UNSEEN[cancer], E3_UNSEEN_SHARE,
                                        N_PER_GRID_POINT, np.random.default_rng([SEED, ci, 10])),
                {"cancer": cancer, "kind": "unseen", "signature": E3_UNSEEN[cancer],
                 "share": E3_UNSEEN_SHARE, "seed": [SEED, ci, 10]}, fingerprints),
            "inverted prior": load_or_make(
                sim_dir / f"e3_{cancer}_inverted.npz",
                lambda: simulate_grid(inverted_prior(stats), signatures, N_PER_GRID_POINT,
                                      np.random.default_rng([SEED, ci, 11])),
                {"cancer": cancer, "kind": "inverted_prior", "seed": [SEED, ci, 11]}, fingerprints),
        }

        test = sets["test"]
        real_counts = catalogs[mask]
        for m in make_methods(S_df):
            t_fit = time.time()
            m.fit(sets["train"], val=sets["val"])
            if hasattr(m, "history_"):
                history.append(pd.DataFrame(m.history_).assign(cancer=cancer, method=m.name,
                                                               train_seconds=round(time.time() - t_fit, 1)))

            # ---- main test set: cosine, MAE, detection ----
            pred = m.predict(test["counts"])
            validate_exposures(pred, names)
            cos = cosine_per_tumor(test["exposures"], pred)
            err = np.abs(np.asarray(test["exposures"]) - pred)
            act = np.asarray(test["exposures"]) >= PRESENCE_THRESHOLD
            per_tumor.append(pd.DataFrame({"cancer": cancer, "method": m.name, "n": test["n_mutations"],
                                           "cosine": cos,
                                           "mae_active": (err * act).sum(1) / np.maximum(act.sum(1), 1)}))
            thr = detection_threshold(sets["val"]["exposures"], m.predict(sets["val"]["counts"]))
            for n in N_GRID:
                sel = test["n_mutations"] == n
                detection.append({"cancer": cancer, "method": m.name, "n": n, "call_threshold": thr,
                                  **detection_report(test["exposures"][sel], pred[sel], thr)})
            lod.append(limit_of_detection(test["exposures"], pred, test["n_mutations"], thr)
                       .assign(cancer=cancer, method=m.name))
            table, _ = presence_report(test["exposures"], pred, names)
            recall_rows.append(table.assign(cancer=cancer, method=m.name))
            if m.name in ("M0_nnls", "M2_multinomial_em", "M3_mlp", "M3_mlp_broad"):
                t, p = np.asarray(test["exposures"]).argmax(1), pred.argmax(1)
                cm = np.zeros((len(names), len(names)))
                np.add.at(cm, (t, p), 1)
                confusion[(cancer, m.name)] = cm  # raw counts; the figure normalizes rows

            # ---- E3 robustness ----
            for kind, data in [("in-distribution", test), *e3.items()]:
                pe = m.predict(data["counts"])
                c = cosine_per_tumor(data["exposures"], pe)
                for lo_, hi_, band in [(0, 100, "25-100"), (250, 10**9, "250-5000")]:
                    sel = (data["n_mutations"] >= lo_) & (data["n_mutations"] <= hi_)
                    robust.append({"cancer": cancer, "method": m.name, "test": kind, "mutations": band,
                                   "median_cosine": float(np.median(c[sel])),
                                   "mae_active": mae_active(data["exposures"][sel], pe[sel])})

            # ---- real tumors (no ground truth) ----
            est = pd.DataFrame(m.predict(real_counts.to_numpy()), index=real_counts.index, columns=names)
            real_fits[(cancer, m.name)] = est
            prof = real_counts.to_numpy(dtype=float)
            recon = est.to_numpy() @ (S_df.to_numpy() / S_df.to_numpy().sum(0)).T
            rc = cosine_per_tumor(prof / prof.sum(1, keepdims=True), recon)
            agreement.append(pd.DataFrame({"cancer": cancer, "method": m.name, "sample": real_counts.index,
                                           "partition": partition[real_counts.index].values,
                                           "n": prof.sum(1), "reconstruction_cosine": rc}))
            if cancer == "LUAD" and clinical is not None and m.name != "B0_cohort_mean":
                for who, idx in [("all patients", est.index),
                                 ("test patients", est.index[partition[est.index] == "test"])]:
                    e5_rows.append({"method": m.name, "patients": who, **smoking_test(est.loc[idx], clinical)})
        print(f"  {cancer}: done ({time.time() - t0:.0f}s)", flush=True)

    # ---------------------------------------------------------------- tables
    per_tumor = pd.concat(per_tumor, ignore_index=True)
    detection = pd.DataFrame(detection)
    lod = pd.concat(lod, ignore_index=True)
    robust = pd.DataFrame(robust)
    history = pd.concat(history, ignore_index=True)
    agreement = pd.concat(agreement, ignore_index=True)
    e5 = pd.DataFrame(e5_rows)
    recall = pd.concat(recall_rows, ignore_index=True)
    # Agreement with M2 on real tumors: cosine between each method's and M2's estimates.
    agree_m2 = []
    for (cancer, mname), est in real_fits.items():
        ref = real_fits[(cancer, "M2_multinomial_em")]
        agree_m2.append(pd.Series(cosine_per_tumor(ref.to_numpy(), est.to_numpy()), index=est.index)
                        .rename("cosine_vs_m2").to_frame().assign(cancer=cancer, method=mname))
    agree_m2 = pd.concat(agree_m2).rename_axis("sample").reset_index()
    agreement = agreement.merge(agree_m2, on=["cancer", "method", "sample"])

    # The team's paired tests (run_checklist.py), filtered to the comparisons that involve M3.
    tests_path = paths["results"] / "method_tests_real.csv"
    paired = pd.read_csv(tests_path)
    if not paired["method_a"].str.startswith("M3").any() and not paired["method_b"].str.startswith("M3").any():
        raise SystemExit(f"{tests_path.name} has no M3 rows: rerun run_checklist.py --signatures ... first.")
    save_results(detection, "eval_detection_real", paths)
    save_results(lod, "eval_limit_of_detection_real", paths)
    save_results(robust, "eval_robustness_real", paths)
    save_results(e5, "eval_e5_methods_real", paths)
    save_results(agreement, "eval_real_agreement_real", paths)
    save_results(history, "m3_training_history_real", paths)
    (paths["results"] / "e3_fingerprints_real.json").write_text(json.dumps(fingerprints, indent=2))

    make_figures(per_tumor, detection, lod, robust, recall, confusion, e5, history, agreement, real_fits, clinical, paths)
    write_summary(per_tumor, paired, detection, robust, recall, e5, agreement, paths)
    print(f"Wrote results/evaluation_summary_real.md ({time.time() - t0:.0f}s)")


# ---------------------------------------------------------------- figures
def _methods_in(df):
    return [m for m in LABELS if m in set(df["method"])]


def make_figures(per_tumor, detection, lod, robust, recall, confusion, e5, history, agreement, real_fits, clinical, paths):
    cancers = list(CANCER_SIGNATURES)

    # Figure 1: accuracy vs mutation count (median and interquartile range), cosine and MAE.
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.2), sharex=True)
    for j, cancer in enumerate(cancers):
        g = per_tumor[per_tumor.cancer == cancer]
        for m in _methods_in(g):
            h = g[g.method == m].groupby("n")
            for i, col in enumerate(["cosine", "mae_active"]):
                q = h[col].quantile([0.25, 0.5, 0.75]).unstack()
                axes[i, j].plot(q.index, q[0.5], marker="o", ms=3, color=COLORS[m], label=LABELS[m])
                axes[i, j].fill_between(q.index, q[0.25], q[0.75], color=COLORS[m], alpha=0.08)
        axes[0, j].axhline(0.9, color="grey", lw=0.8, ls="--")
        axes[0, j].set_title(cancer)
        axes[1, j].set_xlabel("mutations per tumor")
        for i in range(2):
            axes[i, j].set_xscale("log")
            axes[i, j].grid(alpha=0.3)
    axes[0, 0].set_ylabel("cosine, true vs estimated\n(median, IQR band; higher = better)")
    axes[1, 0].set_ylabel("MAE on present signatures\n(median, IQR band; lower = better)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(y=1.03, t="Accuracy vs mutation count on simulated test tumors (dashed line: N* threshold 0.90)")
    save_figure(fig, "eval_fig1_accuracy_vs_mutations_real", paths)
    plt.close(fig)

    # Figure 2: detection, average precision vs N, and the limit-of-detection curve.
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.8))
    for j, cancer in enumerate(cancers):
        g = detection[detection.cancer == cancer]
        for m in _methods_in(g):
            h = g[g.method == m]
            axes[j].plot(h.n, h.avg_precision_macro, marker="o", ms=3, color=COLORS[m], label=LABELS[m])
        axes[j].set_xscale("log")
        axes[j].set_title(f"{cancer}: average precision (per signature)")
        axes[j].set_xlabel("mutations per tumor")
        axes[j].grid(alpha=0.3)
    axes[0].set_ylabel("average precision (presence at 5%)")
    pooled = lod.groupby(["method", "bin_low", "mutations_from_signature"], as_index=False)[["n_pairs"]].sum()
    hits = lod.assign(h=lod.hit_rate * lod.n_pairs).groupby(["method", "bin_low"])["h"].sum()
    pooled["hit_rate"] = hits.values / pooled.n_pairs.values
    for m in _methods_in(pooled):
        h = pooled[pooled.method == m].sort_values("bin_low")
        axes[3].plot(range(len(h)), h.hit_rate, marker="o", ms=3, color=COLORS[m], label=LABELS[m])
    ticks = pooled[pooled.method == "M0_nnls"].sort_values("bin_low")["mutations_from_signature"]
    axes[3].set_xticks(range(len(ticks)), ticks, rotation=35, fontsize=8)
    axes[3].axhline(0.95, color="red", lw=0.8, ls="--")
    axes[3].set_title("Limit of detection (all types)")
    axes[3].set_xlabel("mutations contributed by the signature")
    axes[3].set_ylabel("hit rate at 5% false alarms")
    axes[3].grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    save_figure(fig, "eval_fig2_detection_real", paths)
    plt.close(fig)

    # Figure 3: recall per signature vs how common the signature is (the prior effect).
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.9), sharey=True)
    for j, cancer in enumerate(cancers):
        g = recall[recall.cancer == cancer]
        for m in ["M0_nnls", "M2_multinomial_em", "M3_mlp", "M3_mlp_broad"]:
            h = g[(g.method == m) & (g.support > 0)].sort_values("support")
            prev = h.support / 400
            axes[j].plot(prev, h.recall, marker="o", ms=4, color=COLORS[m], label=LABELS[m], alpha=0.85)
            if m == "M3_mlp":
                for x, y, s in zip(prev, h.recall, h.signature):
                    axes[j].annotate(s, (x, y), fontsize=6.5, xytext=(3, -8), textcoords="offset points")
        axes[j].set_title(cancer)
        axes[j].set_xlabel("share of test tumors where the signature is present")
        axes[j].grid(alpha=0.3)
    axes[0].set_ylabel("recall at the 5% threshold")
    axes[0].legend(fontsize=8)
    fig.suptitle(y=1.03, t="M3 finds common signatures and misses rare ones (it learned the usual mix)")
    save_figure(fig, "eval_fig3_recall_vs_prevalence_real", paths)
    plt.close(fig)

    # Figure 4: M3 robustness checks.
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharey="row")
    tests = ["in-distribution", "unseen signature", "inverted prior"]
    for i, band in enumerate(["25-100", "250-5000"]):
        for j, cancer in enumerate(cancers):
            g = robust[(robust.cancer == cancer) & (robust.mutations == band)]
            ms = _methods_in(g)
            w = 0.8 / len(ms)
            for k, m in enumerate(ms):
                vals = g[g.method == m].set_index("test").loc[tests, "median_cosine"]
                axes[i, j].bar(np.arange(3) + (k - (len(ms) - 1) / 2) * w, vals, w, color=COLORS[m], label=LABELS[m])
            axes[i, j].set_xticks(range(3), tests, fontsize=8)
            axes[i, j].set_ylim(0.0, 1.05)
            axes[i, j].set_title(f"{cancer}, {band} mutations", fontsize=10)
            axes[i, j].grid(axis="y", alpha=0.3)
        axes[i, 0].set_ylabel("median cosine")
    handles, labels_ = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="lower center", ncol=len(labels_), fontsize=8, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle(y=1.03, t=f"M3 robustness checks: {int(E3_UNSEEN_SHARE * 100)}% of mutations from an unseen signature "
                 f"({', '.join(f'{c}: {s}' for c, s in E3_UNSEEN.items())}); inverted prior = common and rare swapped")
    save_figure(fig, "eval_fig4_robustness_real", paths)
    plt.close(fig)

    # Figure 5: dominant-signature confusion (row-normalized), M0, M2, M3, M3 broad.
    # Row labels show how many test tumors have that true dominant signature:
    # rows with only a few tumors look like strong patterns but are not.
    methods = ["M0_nnls", "M2_multinomial_em", "M3_mlp", "M3_mlp_broad"]
    fig, axes = plt.subplots(3, 4, figsize=(17, 12.5))
    for i, cancer in enumerate(cancers):
        names = CANCER_SIGNATURES[cancer]
        for j, m in enumerate(methods):
            counts = confusion[(cancer, m)]
            n_row = counts.sum(1)
            cm = counts / np.maximum(n_row[:, None], 1)
            ax = axes[i, j]
            ax.imshow(cm, cmap="Blues", vmin=0, vmax=1)
            ax.set_xticks(range(len(names)), names, rotation=90, fontsize=7)
            ax.set_yticks(range(len(names)), [f"{s} (n={int(k)})" for s, k in zip(names, n_row)], fontsize=7)
            for a in range(len(names)):
                for b in range(len(names)):
                    if cm[a, b] >= 0.05:
                        ax.text(b, a, f"{cm[a, b]:.2f}", ha="center", va="center", fontsize=6,
                                color="white" if cm[a, b] > 0.6 else "black")
            ax.set_title(f"{cancer}: {LABELS[m]}", fontsize=10)
            if j == 0:
                ax.set_ylabel("true dominant signature")
            if i == 2:
                ax.set_xlabel("estimated dominant signature")
    fig.suptitle(y=1.03, t="Dominant-signature confusion (rows sum to 1; n = test tumors with that true dominant "
                 "signature; all mutation counts pooled)")
    fig.tight_layout()
    save_figure(fig, "eval_fig5_confusion_real", paths)
    plt.close(fig)

    # Figure 6: E5 on real LUAD tumors, every method.
    if clinical is not None:
        ms = ["M0_nnls", "M1_elastic_net", "M2_multinomial_em", "M3_mlp", "M3_mlp_broad"]
        fig, axes = plt.subplots(1, len(ms), figsize=(3.4 * len(ms), 3.8), sharey=True)
        from sigdet.realdata import per_patient
        for ax, m in zip(axes, ms):
            w = per_patient(real_fits[("LUAD", m)])["SBS4"]
            df = clinical.join(w.rename("w"), how="inner").dropna(subset=["smoker"])
            data = [df.loc[df.smoker == 0, "w"], df.loc[df.smoker == 1, "w"]]
            ax.boxplot(data, widths=0.5, showfliers=False)
            for k, d in enumerate(data):
                ax.scatter(np.random.default_rng(0).normal(k + 1, 0.06, len(d)), d, s=4, alpha=0.3, color=COLORS[m])
            row = e5[(e5.method == m) & (e5.patients == "all patients")].iloc[0]
            ax.set_xticks([1, 2], [f"never\n(n={len(data[0])})", f"smokers\n(n={len(data[1])})"])
            ax.set_title(f"{LABELS[m]}\np = {row.mannwhitney_p:.1e}", fontsize=10)
            ax.grid(axis="y", alpha=0.3)
        axes[0].set_ylabel("estimated SBS4 (tobacco) weight")
        fig.suptitle(y=1.10, t="E5, real LUAD tumors: SBS4 weight by smoking status (one-sided Mann-Whitney)")
        save_figure(fig, "eval_fig6_e5_smoking_real", paths)
        plt.close(fig)

    # Figure 7: M3 training curves.
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.6))
    for ax, cancer in zip(axes, cancers):
        for m, ls in [("M3_mlp", "-"), ("M3_mlp_broad", "--")]:
            h = history[(history.cancer == cancer) & (history.method == m)]
            ax.plot(h.epoch, h.train_loss, ls=ls, color="#555555", label=f"{LABELS[m]} train")
            ax.plot(h.epoch, h.val_loss, ls=ls, color=COLORS[m], label=f"{LABELS[m]} validation")
            best = h.loc[h.val_loss.idxmin()]
            ax.scatter([best.epoch], [best.val_loss], color=COLORS[m], zorder=3, s=18)
        ax.set_title(cancer)
        ax.set_xlabel("epoch")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("mean squared error of exposures")
    axes[0].legend(fontsize=7)
    fig.suptitle(y=1.03, t="M3 training (dots: epoch kept by early stopping). The broad curves sit higher because its "
                 "validation set also has harder flat-prior tumors")
    save_figure(fig, "eval_fig7_m3_training_real", paths)
    plt.close(fig)

    # Figure 8: real tumors, how well each method's estimate reconstructs the observed profile.
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.6), sharey=True)
    bins = [0, 50, 100, 250, 500, 1000, np.inf]
    labels = ["<50", "50-100", "100-250", "250-500", "500-1k", "1k+"]
    for ax, cancer in zip(axes, cancers):
        g = agreement[agreement.cancer == cancer].assign(bin=lambda d: pd.cut(d.n, bins, labels=labels, right=False))
        for m in _methods_in(g):
            med = g[g.method == m].groupby("bin", observed=False)["reconstruction_cosine"].median()
            ax.plot(range(len(med)), med.values, marker="o", ms=3, color=COLORS[m], label=LABELS[m])
        ax.set_xticks(range(len(labels)), labels, fontsize=8)
        ax.set_title(cancer)
        ax.set_xlabel("mutations per real tumor")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("median cosine, observed profile\nvs reconstructed profile")
    axes[0].legend(fontsize=7)
    fig.suptitle(y=1.03, t="Real tumors (no ground truth): how well each estimate explains the observed mutations")
    save_figure(fig, "eval_fig8_real_reconstruction_real", paths)
    plt.close(fig)


# ---------------------------------------------------------------- summary
def write_summary(per_tumor, paired, detection, robust, recall, e5, agreement, paths):
    def md(df):
        """A Markdown table without extra packages (requirements.txt has no tabulate)."""
        cols = list(df.columns)
        rows = ["| " + " | ".join(map(str, cols)) + " |", "|" + " --- |" * len(cols)]
        rows += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
        return "\n".join(rows)

    lines = ["# Evaluation summary (real-calibrated simulations and real TCGA tumors)", "",
             "Written by run_evaluation.py.", ""]

    lines += ["## Median cosine by mutation count (simulated test tumors)", ""]
    t = per_tumor.groupby(["cancer", "method", "n"])["cosine"].median().unstack("n").round(3).reset_index()
    lines += [md(t), ""]

    lines += ["## Paired Wilcoxon tests involving M3 (from results/method_tests_real.csv)", "",
              "All mutation counts pooled; Holm-corrected by run_checklist.py over every pair of methods.",
              "Median difference is method A minus method B: positive favours A for cosine, negative",
              "favours A for MAE. Per-count tests are in the same file.", ""]
    m3 = paired[(paired["n"].astype(str) == "all")
                & (paired["method_a"].str.startswith("M3") | paired["method_b"].str.startswith("M3"))].copy()
    m3["p_holm"] = m3["p_holm"].map(lambda p: f"{p:.1e}")
    lines += [md(m3[["cancer", "metric", "method_a", "method_b", "median_difference", "a_higher_fraction",
                     "p_holm", "significant"]].round(4)), ""]

    lines += ["## Detection at 25 and 1,000 mutations", ""]
    d = detection[detection.n.isin([25, 1000])][["cancer", "method", "n", "roc_auc_macro", "avg_precision_macro",
                                                  "hit_rate_at_fa5", "fdr"]].round(3)
    lines += [md(d), ""]

    lines += ["## Recall per signature (5% threshold), M2 vs M3", ""]
    r = recall[recall.method.isin(["M2_multinomial_em", "M3_mlp", "M3_mlp_broad"])]
    r = r.pivot_table(index=["cancer", "signature", "support"], columns="method", values="recall").round(2).reset_index()
    lines += [md(r.sort_values(["cancer", "support"])), ""]

    lines += ["## M3 robustness checks (median cosine; not the team's E3, which removes a listed signature)", ""]
    rb = robust.pivot_table(index=["cancer", "mutations", "method"], columns="test", values="median_cosine").round(3).reset_index()
    lines += [md(rb), ""]

    lines += ["## E5 smoking test, every method (real LUAD)", ""]
    cols = [c for c in ["method", "patients", "n_smokers", "n_never", "median_smokers", "median_never",
                        "mannwhitney_p", "spearman_rho_pack_years", "spearman_p"] if c in e5.columns]
    e5 = e5.copy()
    for c in ("mannwhitney_p", "spearman_p"):
        if c in e5.columns:
            e5[c] = e5[c].map(lambda p: f"{p:.1e}")
    lines += [md(e5[cols].round(4)), ""]

    lines += ["## Real tumors: reconstruction and agreement with M2 (median over tumors)", ""]
    a = agreement.groupby(["cancer", "method"])[["reconstruction_cosine", "cosine_vs_m2"]].median().round(3).reset_index()
    lines += [md(a), ""]

    lines += ["## Figures", ""] + [f"- figures/{p.name}" for p in sorted(paths["figures"].glob("eval_*.png"))]
    (paths["results"] / "evaluation_summary_real.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
