# M3 neural network, evaluation and figures (Bao Tran)

All numbers come from `run_checklist.py` and `run_evaluation.py` with
`--signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt`, on the saved patient split and simulated
sets. Tables: `results/run_summary_real.md`, `results/method_tests_real.csv` and
`results/evaluation_summary_real.md`. Figures: `figures/eval_fig*_real.png`.

## 1. What we built

**M3 (`sigdet/methods.py`, class `M3MLP`).** A small neural network built on scikit-learn's `MLPRegressor`
(already in `requirements.txt`), with the same `fit(train, val)` / `predict(counts)` interface as B0 to M2.
It is in the `methods` list of `run_checklist.py`, so it is scored on the same test tumors as every other
method and is included in `compare_methods`, the team's paired Wilcoxon tests with Holm's correction.

- Input: the square root of each context's share of the tumor's mutations (96 values) and
  log10(mutation count) / 4. The square root evens out count noise; the count tells the network how noisy
  the profile is.
- Two hidden layers (128 and 64 units, ReLU), one output per signature, squared-error loss, Adam,
  L2 penalty `alpha = 0.1`.
- **Valid output.** The output layer is linear, so predictions can be negative and need not sum to 1. We set
  negatives to 0 and divide each row by its sum, which gives a valid exposure vector (Q1).
- **Early stopping on our validation set (Q2).** We do not use `MLPRegressor(early_stopping=True)`: it holds
  out a random 10% of the training rows, and with thinning, copies of the same tumor would sit on both sides.
  We train one epoch at a time with `partial_fit` and keep the model with the lowest validation error.
- **Training data.** The 1,000 simulated training tumors per cancer type, augmented by **thinning** (a random
  subset of each tumor's mutations, 8 copies per tumor, counts spread on a log scale down to 10). A random
  subset of a tumor's mutations has the same signature mix, so the labels stay exact.
- **M3_mlp_broad** adds 4,000 tumors simulated with a flat prior over signatures (Dirichlet 0.5, same
  `simulate()` and noise settings), to test the risk that "M3 may learn the simulator rather than tumors".
- **Tuning.** Settings chosen on validation error only (`tune_m3.py`, `results/m3_tuning_real.csv`,
  20 settings per cancer type). Validation error is flat between `alpha` 0.03 and 0.1 (within 0.1%).
- Training takes 3 to 7 seconds per cancer type for M3 and 11 to 79 seconds for M3 broad, stopping after
  34 to 136 epochs (`figures/eval_fig7_m3_training_real.png`). Training is deterministic.

**Evaluation (`sigdet/evaluation.py`, `run_evaluation.py`).** On the same split and simulated sets, for
B0, M0, M1, M2, M3 and M3 broad: per-signature detection metrics, a limit of detection, two M3 robustness
checks, E5 for every method, and checks on real tumors. It writes only new files. The paired tests are not
repeated: the summary quotes the M3 rows of `results/method_tests_real.csv`.

## 2. Results

### 2.1 Accuracy and breakdown point (simulated test tumors)

| | M0 NNLS | M1 elastic-net | M2 EM | M3 | M3 broad |
| --- | --- | --- | --- | --- | --- |
| N*, LUAD / SKCM / COAD | 50 / 50 / 250 | 50 / 50 / 250 | 50 / 25 / 100 | 25 / 25 / 25 | 25 / 25 / 25 |
| MAE on present signatures, LUAD | 0.089 | 0.080 | 0.080 | 0.074 | 0.072 |
| MAE, SKCM | 0.090 | 0.088 | 0.073 | 0.066 | 0.063 |
| MAE, COAD | 0.101 | 0.092 | 0.082 | 0.072 | 0.065 |
| Macro F1, LUAD / SKCM / COAD | 0.870 / 0.762 / 0.749 | 0.880 / 0.756 / 0.758 | 0.873 / 0.800 / 0.779 | 0.878 / 0.778 / **0.645** | 0.879 / **0.814** / **0.786** |

Both networks have the lowest error and reach the 0.90 cosine threshold at the smallest tested count (25)
in every cancer type. The paired tests (`compare_methods`, Holm-corrected over all 15 method pairs) show
**where** the gain is:

- Pooled over mutation counts, M3 beats M2 on cosine only in COAD (median difference 0.005); in LUAD and
  SKCM the difference is not significant. M3 broad beats M2 in COAD (0.012) and SKCM (0.002).
- Per mutation count, M3 and M3 broad beat M2 significantly only at 25 mutations in COAD and SKCM, and at
  50 in COAD. From 100 mutations up, no difference is significant. In LUAD none is significant at any count.
- At 5,000 mutations both networks are slightly behind M2 (median cosine difference -0.005 to -0.022, not
  significant), probably because few training tumors have that many mutations; thinning only lowers counts.
- M3 broad beats plain M3 pooled in COAD (median difference 0.005, Holm p = 3e-4), not in LUAD or SKCM.

### 2.2 Where M3's gain comes from: it learned the usual signature mix

- **Recall per signature** (`figures/eval_fig3_recall_vs_prevalence_real.png`): plain M3 finds common
  signatures almost always (recall 0.99 to 1.00 for SBS1, SBS5 and SBS40a in COAD) but misses rare ones. In
  COAD it **never** detects SBS17b (recall 0.00; M1 0.92, M2 0.86) and finds SBS13 in 0.34 (M2 0.81); in
  SKCM it finds SBS38 and SBS7c in about 0.40 (M2 0.76). That gives plain M3 the lowest macro F1 of the
  non-baseline methods in COAD (0.645). M3 broad recovers most of it (SBS17b 0.78, SBS13 0.70) and has the
  best macro F1 of all methods in SKCM and COAD.
- **Detection per signature** (`figures/eval_fig2_detection_real.png`; B0 scores 0.5): at 25 mutations the
  networks are about level with M2 (ROC-AUC 0.75 to 0.85 for all three). At 1,000 mutations M2 is best in all
  three types (0.948 to 0.968), M3 broad within 0.007, plain M3 lower (0.914 to 0.952). So M3 does not detect
  signatures in individual tumors better; its lower error comes largely from predicting the typical mix.
- **Confusion** (`figures/eval_fig5_confusion_real.png`): for the 48 COAD tumors whose main signature is
  SBS40a, plain M3 names SBS40a in 10% of them (M0 52%, M2 54%, M3 broad 44%). Row labels show the number of
  tumors per row; rows with 1 or 2 tumors should not be read as patterns.

### 2.3 M3 robustness checks (`figures/eval_fig4_robustness_real.png`)

These are M3-specific checks, not the team's E3 (which removes a signature from the candidate list).
Median cosine at 250 to 5,000 mutations:

| Test | M0 | M1 | M2 | M3 | M3 broad |
| --- | --- | --- | --- | --- | --- |
| In-distribution, COAD | 0.953 | 0.938 | 0.971 | 0.962 | 0.974 |
| Inverted prior, COAD | 0.986 | 0.988 | 0.990 | **0.350** | 0.985 |
| Inverted prior, SKCM | 0.988 | 0.989 | 0.993 | **0.430** | 0.989 |
| Inverted prior, LUAD | 0.985 | 0.985 | 0.992 | 0.960 | 0.988 |
| Unseen signature (20% SBS10a), COAD | 0.852 | 0.916 | **0.834** | 0.957 | 0.944 |

- **Inverted prior** (common and rare signatures swapped): plain M3 collapses in COAD and SKCM; LUAD, whose
  mix is more even, suffers less. M0, M1, M2 and M3 broad are unaffected. This is the clearest evidence that
  plain M3 relies on the training mix.
- **Unseen signature** (20% of mutations from a signature outside the list): M2 is hit hardest because it
  explains every mutation with the listed signatures as faithfully as possible. M1's shrinkage halves that
  loss, and the networks lose least because they partly ignore mutations that do not fit their expected mix.
  B0 is also "robust" here (0.904 in COAD) for the same reason, so this is not evidence that M3 reads tumors
  better.

### 2.4 Real tumors (E5 and reconstruction)

| LUAD, all patients | M0 | M1 | M2 | M3 | M3 broad |
| --- | --- | --- | --- | --- | --- |
| Median SBS4, smokers (n = 424) | 0.503 | 0.499 | 0.514 | 0.458 | 0.450 |
| Median SBS4, never-smokers (n = 69) | 0.029 | 0.045 | 0.036 | **0.156** | 0.070 |
| One-sided Mann-Whitney p | 3.7e-21 | 4.2e-21 | 8.2e-22 | 1.7e-20 | 3.3e-21 |

All methods separate smokers from never-smokers strongly, also on test-partition patients only (p = 0.002 to
0.003, 83 vs 11 patients). But plain M3 gives **never-smokers a median tobacco weight of 0.156**, five times
M0's (`figures/eval_fig6_e5_smoking_real.png`): the learned "typical lung tumor" leaks into individual
estimates. M3 broad halves it. No method finds a significant pack-years correlation (Spearman p = 0.16 to 0.61).

On test-partition real tumors, the networks agree closely with M2 (median cosine 0.960 to 0.994) and
reconstruct the observed profiles slightly less well than M0 and M2 (0.889 vs 0.896 in LUAD), as expected,
since NNLS and EM directly maximize the fit (`figures/eval_fig8_real_reconstruction_real.png`).

## 3. What we conclude

1. A neural network trained on simulated tumors gives the lowest error and the lowest breakdown point, but
   it beats M2 significantly only at 25 to 50 mutations, and the gain comes from the prior it learns, not
   from better detection in individual tumors.
2. That prior is a liability when tumors differ from the training mix: rare signatures are missed (SBS17b
   never), an inverted mix breaks the network, and on real lung tumors it assigns tobacco damage to
   never-smokers.
3. Training on extra flat-prior tumors (M3 broad) removes most of these failures, keeps the low-count
   advantage, and gives the best macro F1 in SKCM and COAD. If we report one network, it should be M3 broad.
4. For conclusions about individual real tumors, M2 remains the safest method: never far behind the
   networks on accuracy, best at detection above a few hundred mutations, and no learned prior.

## 4. Limitations

- All accuracy numbers come from our simulator; the networks are trained and tested on the same simulator
  family, which favours them. The inverted-prior check measures part of that advantage, not all of it.
- One training run per network. M3's numbers can also differ in the third decimal between computers
  (floating-point differences in the math libraries); B0 to M2 reproduce exactly.
- Each mutation count has 50 test tumors per cancer type, so per-count detection metrics for rare
  signatures are noisy and have no error bars.
- The robustness checks use one setting each: 20% from one outside signature, and one inverted mix.
- We did not run E5 negative controls for M3: a network trained without UV signatures cannot assign weight
  to them, and one trained with them would learn they are almost always absent, so the control tests nothing.

## 5. Notes for the team

- **Holm family size.** `compare_methods` corrects over every method pair, so adding M3 (15 pairs instead
  of 6) makes all p-values more conservative. Pooled results are unchanged, and the progress report's
  claim (M2 beats M0 at 5 of 8 counts in SKCM, 3 in COAD, none in LUAD) still holds, but 16 borderline
  per-count comparisons among B0 to M2 (Holm p 0.026 to 0.048 before) are now above 0.05. Hugo may prefer
  a smaller family (for example, only comparisons against M2); that is a choice for the statistics owner.
- **Overlap with planned work.** The unseen-signature check overlaps Hugo's E3, and the inverted prior is
  close to his planned sparse-tumor (`n_active`) stress test. If his E3 lands, the report should use his
  version and keep the inverted prior as the M3-specific check. E5 for every method complements Anas's
  extended E5.
- **Fingerprints across machines.** With a different NumPy version, `splits/sim_fingerprints_*.json` change
  after a run although results agree to about 1e-16; restore them with `git checkout splits/`.

## 6. How to reproduce

```
pip install -r requirements.txt
pytest -q                                                                               # 37 tests
python run_checklist.py  --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt   # about 3 minutes
python run_evaluation.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt   # about 3 minutes
python tune_m3.py        --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt   # optional, about 4 minutes
```
