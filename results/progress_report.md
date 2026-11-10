# Progress report: Mutational Signature Detective (draft, 2026-11-10)

Team: Anas Farooqui (40347862, contact), Hugo Xu (40064100), Bao Tran Nguyen (40257379).
Repository: github.com/ralegionc/COMP6721-F26-PROJECT. Every number below comes from
`results/run_summary_real.md` and the CSV files next to it, produced with
`python run_checklist.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt`.

## 1. Summary

We estimate how much each known COSMIC mutational signature contributed to a tumor, using the tumor's
single-base substitutions counted over the 96 standard contexts (the substituted base plus one
neighbouring base on each side). We have built the full data pipeline for three TCGA cancer types, a
simulator that produces tumors with known answers, four of our five planned methods, and a statistical
comparison between them. On simulated tumors, the maximum-likelihood method (M2) is the most accurate
in all three cancer types, and the improvement over plain non-negative least squares (NNLS) is
statistically significant but small. On real lung tumors, the estimated tobacco signature (SBS4) is
much higher in smokers than in lifelong non-smokers. The main open problems are that our simulated
tumors are easier than real ones, and that signatures which should be absent still receive weight in
real tumors with few mutations.

## 2. Work completed

### 2.1 Data (Week 8)

| Step | Result |
| --- | --- |
| Download | 1,547 open-access TCGA "Masked Somatic Mutation" MAF files (LUAD 618, SKCM 472, COAD 457; 239 MB), each checked against its GDC MD5 checksum, and 1,516 clinical records. `data/manifest.csv` lists every file with its SHA-256 checksum. |
| Preprocessing | We keep single-base substitutions called by MuTect2 in tumor samples, remove duplicates, and fold each onto the pyrimidine strand to give 96 context counts. Result: 723,536 mutations in 1,536 samples from 1,451 patients. |
| Genome-build check | None of the 723,536 reference bases in the MAFs disagrees with the GRCh38 genome. |
| Signature lists | Non-artifact COSMIC v3.4 signatures that Alexandrov et al. (2020) found in at least 5% of each type's TCGA exomes: LUAD 6, SKCM 7, COAD 8. |
| Clinical labels | LUAD: 433 current or reformed smokers, 75 lifelong non-smokers, 77 with unknown status (excluded). |
| Split | 871 / 290 / 290 patients (train / validation / test), by 12-character patient ID, stratified by cancer type, seed 6721, saved once with a SHA-256 fingerprint. |

Median mutation counts per sample are 194 (LUAD), 390 (SKCM) and 121 (COAD). Between 18% (SKCM) and 34% (COAD)
of samples per type have fewer than 100 mutations, which is the range where the breakdown analysis matters.

### 2.2 Simulator and methods (Weeks 8 and 9)

The simulator draws each tumor's exposures from a Dirichlet distribution and its mutation count from a
log-normal distribution, both fitted to NNLS estimates on training patients only. It then adds two
sources of realism: extra count noise beyond the multinomial model (Dirichlet-multinomial sampling) and
a random distortion of each signature per tumor. For every cancer type we simulate 1,000 training
tumors, 200 validation tumors and 400 test tumors (50 at each of eight mutation counts from 25 to
5,000). A sparse option (`n_active`) gives tumors with only one to three active signatures; it came
from our earlier standalone simulator and is now part of the shared one.

| Method | Description |
| --- | --- |
| B0 cohort mean | Predicts the average training exposure for every tumor. It is the attribution equivalent of always guessing the majority class, and every method must beat it. |
| M0 NNLS | Non-negative least squares fit of the tumor's mutation profile. |
| M1 elastic-net NNLS | NNLS with an L1 penalty (which sets small weights to zero) and an L2 penalty (which shrinks weights). Penalty strength (1e-7 to 1e-3) and the L1 share (0.1, 0.5, 0.9) are chosen on validation tumors by mean absolute error on truly present signatures. |
| M2 multinomial EM | Maximum-likelihood fit of a multinomial count model with the expectation-maximization algorithm, with bootstrap intervals. |
| M3 MLP | Not finished; see Section 4. |

In Week 9 we merged a second, standalone simulator and classical-methods pipeline into the shared code
in `sigdet/`. The standalone version sampled from all 86 COSMIC signatures, including 19 sequencing
artifacts, ignored the patient split, called a signature present at a weight of 1e-6, and averaged
error over all 86 signatures. Those choices made precision look much worse (0.13 at the 1e-6 threshold
against 0.72 at 5% for the same estimates) and error look much better than under our agreed design.
Its statistical tests were correct and are now part of `sigdet/metrics.py`.

### 2.3 Results on simulated test tumors

All mutation counts are pooled. MAE is the mean absolute error on signatures that are truly present
(weight at least 5%); macro F1 and PR-AUC (area under the precision-recall curve) measure whether each
signature is correctly called present or absent.

| Cancer | Method | Median cosine | MAE (present) | Macro F1 | PR-AUC | N* |
| --- | --- | --- | --- | --- | --- | --- |
| LUAD | B0 | 0.859 | 0.131 | 0.647 | 0.614 | not reached |
| LUAD | M0 | 0.959 | 0.089 | 0.870 | 0.940 | 50 |
| LUAD | M1 | 0.961 | 0.080 | 0.880 | 0.942 | 50 |
| LUAD | M2 | 0.966 | 0.080 | 0.873 | 0.946 | 50 |
| SKCM | B0 | 0.931 | 0.111 | 0.409 | 0.473 | 25 |
| SKCM | M0 | 0.966 | 0.090 | 0.762 | 0.840 | 50 |
| SKCM | M1 | 0.968 | 0.088 | 0.756 | 0.838 | 50 |
| SKCM | M2 | 0.980 | 0.073 | 0.800 | 0.864 | 25 |
| COAD | B0 | 0.913 | 0.100 | 0.444 | 0.492 | 2,500 |
| COAD | M0 | 0.930 | 0.101 | 0.749 | 0.822 | 250 |
| COAD | M1 | 0.919 | 0.092 | 0.758 | 0.825 | 250 |
| COAD | M2 | 0.955 | 0.082 | 0.779 | 0.844 | 100 |

N* is the breakdown point: the smallest tested mutation count from which the median cosine similarity
between true and estimated exposures stays at or above 0.90. Figure:
`figures/cosine_vs_mutation_count_real.png`.

**M1 after widening the grid.** With the earlier grid (1e-7 to 1e-4, L1 share fixed at 0.9), validation
always picked the weakest penalty and M1 reproduced M0. With the wider grid, validation picks a penalty
of 3e-5 (LUAD, SKCM) or 1e-4 (COAD), always with the smallest L1 share (0.1), so the chosen penalty is
mostly L2 shrinkage. M1 now lowers MAE compared with M0 in all three types (most in LUAD and COAD), but
in COAD it lowers the median cosine (0.919 against 0.930). Because the L1 share sits at the edge of the
grid, a pure L2 penalty may do slightly better still.

**Statistical comparison.** For each cancer type we ran a paired Wilcoxon signed-rank test for every
pair of methods on the same 400 test tumors, for cosine and for MAE, with Holm's correction for
multiple comparisons (`results/method_tests_real.csv`). With all mutation counts pooled:

- M2 is better than M0 on both cosine and MAE in all three types (Holm-adjusted p from 6e-6 down to 9e-24).
  The median gain is small: 0.005 (LUAD) to 0.017 (COAD) in cosine.
- M2 is better than M1 in SKCM and COAD; in LUAD the two do not differ (p = 0.5 for cosine).
- M1 is better than M0 on MAE in all three types, by 0.002 to 0.007.
- In COAD, M0 is not significantly better than B0 on cosine (p = 1) or MAE (p = 0.87), although its
  F1 is much higher (0.749 against 0.444). In COAD, exposure cosine alone cannot separate a real method
  from the baseline.

When tested separately at each mutation count (50 tumors per test), M2's cosine advantage over M0 is
significant at 5 of 8 counts in SKCM, 3 of 8 in COAD and none in LUAD, so the pooled LUAD result comes
from many small, consistent differences.

### 2.4 Real tumors (E5)

Real tumors have no known answer, so we use checks whose direction is known.

- **Smoking.** NNLS gives smokers' lung tumors a median SBS4 weight of 0.50, against 0.03 for lifelong
  non-smokers (424 against 69 patients; one-sided Mann-Whitney p = 3.7e-21). Among smokers, SBS4 weight
  does not correlate with reported pack-years (Spearman rho = 0.08, p = 0.16, n = 348).
- **Negative controls.** When we add the UV signatures (SBS7a to SBS7d) to the LUAD fit, they receive
  at least 5% weight in 29% of lung tumors, which should carry no UV damage. Adding SBS4 to the melanoma
  fit gives the same result in 9% of SKCM tumors.

## 3. Limitations

- **Simulated accuracy is an upper bound.** Our simulated tumors follow the same signature model the
  methods assume, so breakdown points of 25 to 250 mutations are optimistic. A pilot in our project
  plan found that NNLS needed 500 to 1,000 mutations to reach a cosine of 0.90.
- **Cosine flatters the baseline.** B0 reaches a median cosine of 0.93 in SKCM, where SBS7a dominates
  most tumors, and is statistically indistinguishable from NNLS on cosine in COAD. We therefore treat
  presence F1 and PR-AUC as the main measures of whether a method identifies signatures.
- **Absent signatures absorb weight in real tumors.** The negative controls show that a 5% presence
  threshold is too low for real tumors with few mutations.
- **Data.** Open-access MAFs remove some low-frequency variants, which lowers counts in low-burden
  tumors. Signatures seen in fewer than 5% of a type's tumors (for example the POLE signatures SBS10a
  and SBS10b in about 2% of colorectal tumors) are not candidates, so tumors driven by them will be fitted
  poorly. The published attributions predate COSMIC v3.4, so newer signatures are not candidates.
- **M2 intervals not yet evaluated.** M2 produces bootstrap intervals, but we have not yet checked how
  often they contain the true value.

## 4. Future work (Weeks 10 to 13)

| Week | Work | Owner |
| --- | --- | --- |
| 10 | Finish M3, a small multilayer perceptron trained on the simulated training tumors and tuned on the validation tumors, using the same `fit` / `predict` interface so it is scored on the same test tumors. | Bao Tran |
| 10 | E3: refit with a signature removed from the candidate list (for example SBS10a/b or SBS17a) and measure the loss in accuracy; sparse-tumor stress test with `n_active`; extend M1's grid to a pure L2 penalty. | Hugo |
| 10 to 11 | E5 extended: more negative controls and a presence threshold that rises as the mutation count falls, chosen on simulated validation tumors. | Anas |
| 11 | Full comparison including M3, rerun of the statistical tests, and coverage of M2's bootstrap intervals (how often a 95% interval contains the true weight). | Hugo |
| 12 | Error analysis: which signatures are confused with each other (SBS1, SBS5 and SBS40a have similar flat profiles) and how errors depend on mutation count; final figures. | Bao Tran |
| 13 | Final report, reproducibility check from a clean clone of the repository, presentation. | All |

Optional, if time allows: E4 cohort shift by tissue source site, and a Bayesian fit with Markov chain
Monte Carlo sampling for comparison with M2's bootstrap intervals.

**Risks.** (1) If M3 only learns the simulator, it may score well on simulated tumors and badly on
real ones; we will compare it with M0 on the E5 real-data checks as well. (2) E3 may show much larger
errors than E1; that is a finding to report, not a failure. (3) Late changes to `sigdet/config.py` or
the saved split would change every number, so both stay fixed.

## References

Alexandrov, L. B., et al. (2020). The repertoire of mutational signatures in human cancer. *Nature*,
578, 94-101.

COSMIC Mutational Signatures v3.4. https://cancer.sanger.ac.uk/signatures/

Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of
Statistics*, 6(2), 65-70.
