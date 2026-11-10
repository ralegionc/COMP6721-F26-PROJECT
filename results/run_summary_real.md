# Run summary (real data, 2026-11-10)

## Checklist
- Q1: Exposure vector per tumor: non-negative shares over the cancer type's candidate COSMIC SBS signatures (artifacts excluded), summing to 1. Known only for simulated tumors; real tumors give directional checks only. Subsets: LUAD=6, SKCM=7, COAD=8
- Q2: Split by patient (train=871, val=290, test=290); no overlap; multi-sample patients kept together; simulator calibrated on training patients only.
- Q3: Presence: signature weight >= 5%. Smoking: current or reformed smoker (negative: lifelong non-smoker; unknown excluded). LUAD labels: 433 positive, 75 negative, 77 excluded.
- Q4: Max 'absent' rate 53% (SKCM); accuracy WOULD overstate success. Reported: per-signature F1/PR-AUC, cosine, MAE vs B0 (M0 macro F1 0.75-0.87).
- Q5: Split + 9 simulated sets saved once (seed 6721, SHA-256); all 4 methods scored on the same files.
- Q6: results/, figures/, splits/ committed; data/raw and data/sims gitignored (manifest + seeds).

## Test-set results (simulated tumors, all mutation counts pooled)

| cancer | method | median cosine | MAE (active) | macro F1 | PR-AUC |
| --- | --- | --- | --- | --- | --- |
| LUAD | B0_cohort_mean | 0.859 | 0.131 | 0.647 | 0.614 |
| LUAD | M0_nnls | 0.959 | 0.089 | 0.870 | 0.940 |
| LUAD | M1_elastic_net | 0.961 | 0.080 | 0.880 | 0.942 |
| LUAD | M2_multinomial_em | 0.966 | 0.080 | 0.873 | 0.946 |
| SKCM | B0_cohort_mean | 0.931 | 0.111 | 0.409 | 0.473 |
| SKCM | M0_nnls | 0.966 | 0.090 | 0.762 | 0.840 |
| SKCM | M1_elastic_net | 0.968 | 0.088 | 0.756 | 0.838 |
| SKCM | M2_multinomial_em | 0.980 | 0.073 | 0.800 | 0.864 |
| COAD | B0_cohort_mean | 0.913 | 0.100 | 0.444 | 0.492 |
| COAD | M0_nnls | 0.930 | 0.101 | 0.749 | 0.822 |
| COAD | M1_elastic_net | 0.919 | 0.092 | 0.758 | 0.825 |
| COAD | M2_multinomial_em | 0.955 | 0.082 | 0.779 | 0.844 |

## Breakdown point N* (smallest N with median cosine >= 0.90)

| method | COAD | LUAD | SKCM |
| --- | --- | --- | --- |
| B0_cohort_mean | 2500 | not reached | 25 |
| M0_nnls | 250 | 50 | 50 |
| M1_elastic_net | 250 | 50 | 50 |
| M2_multinomial_em | 100 | 50 | 25 |

M1 penalty chosen on validation data: {'LUAD': {'alpha': 3e-05, 'l1_ratio': 0.1}, 'SKCM': {'alpha': 3e-05, 'l1_ratio': 0.1}, 'COAD': {'alpha': 0.0001, 'l1_ratio': 0.1}}

## Paired Wilcoxon tests (all mutation counts pooled, Holm-corrected)

Median difference is method A minus method B; "A higher" is the share of tumors where A scores higher. Higher cosine and lower MAE are better.

| cancer | metric | A | B | median diff | A higher | p (Holm) |
| --- | --- | --- | --- | --- | --- | --- |
| COAD | cosine | B0_cohort_mean | M0_nnls | -0.0118 | 44% | 1 |
| COAD | mae_active | B0_cohort_mean | M0_nnls | +0.0054 | 52% | 0.87 |
| COAD | cosine | B0_cohort_mean | M1_elastic_net | -0.0008 | 50% | 1 |
| COAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0094 | 57% | 0.0031 |
| COAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.0276 | 35% | 1.8e-07 |
| COAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0211 | 64% | 3.8e-09 |
| COAD | cosine | M0_nnls | M1_elastic_net | +0.0068 | 53% | 1 |
| COAD | mae_active | M0_nnls | M1_elastic_net | +0.0066 | 55% | 0.0067 |
| COAD | cosine | M0_nnls | M2_multinomial_em | -0.0171 | 35% | 3.4e-11 |
| COAD | mae_active | M0_nnls | M2_multinomial_em | +0.0151 | 66% | 1.1e-13 |
| COAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0238 | 31% | 3.5e-12 |
| COAD | mae_active | M1_elastic_net | M2_multinomial_em | +0.0124 | 64% | 3.5e-06 |
| LUAD | cosine | B0_cohort_mean | M0_nnls | -0.0827 | 24% | 3.1e-29 |
| LUAD | mae_active | B0_cohort_mean | M0_nnls | +0.0445 | 72% | 2.3e-23 |
| LUAD | cosine | B0_cohort_mean | M1_elastic_net | -0.0885 | 17% | 4.1e-40 |
| LUAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0517 | 79% | 9.3e-36 |
| LUAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.0946 | 18% | 1.6e-39 |
| LUAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0562 | 77% | 1.2e-31 |
| LUAD | cosine | M0_nnls | M1_elastic_net | -0.0033 | 38% | 1.5e-08 |
| LUAD | mae_active | M0_nnls | M1_elastic_net | +0.0073 | 66% | 5e-09 |
| LUAD | cosine | M0_nnls | M2_multinomial_em | -0.0047 | 39% | 4.9e-06 |
| LUAD | mae_active | M0_nnls | M2_multinomial_em | +0.0080 | 61% | 6e-06 |
| LUAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0014 | 46% | 0.5 |
| LUAD | mae_active | M1_elastic_net | M2_multinomial_em | +0.0007 | 51% | 0.87 |
| SKCM | cosine | B0_cohort_mean | M0_nnls | -0.0263 | 34% | 4.8e-12 |
| SKCM | mae_active | B0_cohort_mean | M0_nnls | +0.0250 | 64% | 4.8e-08 |
| SKCM | cosine | B0_cohort_mean | M1_elastic_net | -0.0272 | 32% | 2.2e-14 |
| SKCM | mae_active | B0_cohort_mean | M1_elastic_net | +0.0266 | 66% | 6.2e-10 |
| SKCM | cosine | B0_cohort_mean | M2_multinomial_em | -0.0415 | 24% | 4.2e-25 |
| SKCM | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0424 | 73% | 4.1e-22 |
| SKCM | cosine | M0_nnls | M1_elastic_net | -0.0017 | 38% | 1.4e-06 |
| SKCM | mae_active | M0_nnls | M1_elastic_net | +0.0015 | 56% | 0.00029 |
| SKCM | cosine | M0_nnls | M2_multinomial_em | -0.0103 | 28% | 8.7e-24 |
| SKCM | mae_active | M0_nnls | M2_multinomial_em | +0.0154 | 72% | 7.6e-21 |
| SKCM | cosine | M1_elastic_net | M2_multinomial_em | -0.0090 | 33% | 4.3e-17 |
| SKCM | mae_active | M1_elastic_net | M2_multinomial_em | +0.0148 | 68% | 5.4e-15 |

## E5 smoking test (LUAD, real data)

{
 "n_smokers": 424,
 "n_never": 69,
 "median_smokers": 0.503240189366239,
 "median_never": 0.028861518922623154,
 "mannwhitney_U": 24904.5,
 "mannwhitney_p": 3.678330070371453e-21,
 "spearman_rho_pack_years": 0.07544613192581809,
 "spearman_p": 0.16020935141980677,
 "n_pack_years": 348
}
