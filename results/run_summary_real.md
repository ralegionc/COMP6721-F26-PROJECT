# Run summary (real data, 2026-11-16)

## Checklist
- Q1: Exposure vector per tumor: non-negative shares over the cancer type's candidate COSMIC SBS signatures (artifacts excluded), summing to 1. Known only for simulated tumors; real tumors give directional checks only. Subsets: LUAD=6, SKCM=7, COAD=8
- Q2: Split by patient (train=871, val=290, test=290); no overlap; multi-sample patients kept together; simulator calibrated on training patients only.
- Q3: Presence: signature weight >= 5%. Smoking: current or reformed smoker (negative: lifelong non-smoker; unknown excluded). LUAD labels: 433 positive, 75 negative, 77 excluded.
- Q4: Max 'absent' rate 53% (SKCM); accuracy WOULD overstate success. Reported: per-signature F1/PR-AUC, cosine, MAE vs B0 (M0 macro F1 0.75-0.87).
- Q5: Split + 9 simulated sets saved once (seed 6721, SHA-256); all 6 methods scored on the same files.
- Q6: results/, figures/, splits/ committed; data/raw and data/sims gitignored (manifest + seeds).

## Test-set results (simulated tumors, all mutation counts pooled)

| cancer | method | median cosine | MAE (active) | macro F1 | PR-AUC |
| --- | --- | --- | --- | --- | --- |
| LUAD | B0_cohort_mean | 0.859 | 0.131 | 0.647 | 0.614 |
| LUAD | M0_nnls | 0.959 | 0.089 | 0.870 | 0.940 |
| LUAD | M1_elastic_net | 0.961 | 0.080 | 0.880 | 0.942 |
| LUAD | M2_multinomial_em | 0.966 | 0.080 | 0.873 | 0.946 |
| LUAD | M3_mlp | 0.964 | 0.074 | 0.878 | 0.941 |
| LUAD | M3_mlp_broad | 0.969 | 0.072 | 0.879 | 0.944 |
| SKCM | B0_cohort_mean | 0.931 | 0.111 | 0.409 | 0.473 |
| SKCM | M0_nnls | 0.966 | 0.090 | 0.762 | 0.840 |
| SKCM | M1_elastic_net | 0.968 | 0.088 | 0.756 | 0.838 |
| SKCM | M2_multinomial_em | 0.980 | 0.073 | 0.800 | 0.864 |
| SKCM | M3_mlp | 0.979 | 0.066 | 0.778 | 0.872 |
| SKCM | M3_mlp_broad | 0.982 | 0.063 | 0.814 | 0.877 |
| COAD | B0_cohort_mean | 0.913 | 0.100 | 0.444 | 0.492 |
| COAD | M0_nnls | 0.930 | 0.101 | 0.749 | 0.822 |
| COAD | M1_elastic_net | 0.919 | 0.092 | 0.758 | 0.825 |
| COAD | M2_multinomial_em | 0.955 | 0.082 | 0.779 | 0.844 |
| COAD | M3_mlp | 0.958 | 0.072 | 0.645 | 0.824 |
| COAD | M3_mlp_broad | 0.965 | 0.065 | 0.786 | 0.839 |

## Breakdown point N* (smallest N with median cosine >= 0.90)

| method | COAD | LUAD | SKCM |
| --- | --- | --- | --- |
| B0_cohort_mean | 2500 | not reached | 25 |
| M0_nnls | 250 | 50 | 50 |
| M1_elastic_net | 250 | 50 | 50 |
| M2_multinomial_em | 100 | 50 | 25 |
| M3_mlp | 25 | 25 | 25 |
| M3_mlp_broad | 25 | 25 | 25 |

M1 penalty chosen on validation data: {'LUAD': {'alpha': 3e-05, 'l1_ratio': 0.1}, 'SKCM': {'alpha': 3e-05, 'l1_ratio': 0.1}, 'COAD': {'alpha': 0.0001, 'l1_ratio': 0.1}}

## Paired Wilcoxon tests (all mutation counts pooled, Holm-corrected)

Median difference is method A minus method B; "A higher" is the share of tumors where A scores higher. Higher cosine and lower MAE are better.

| cancer | metric | A | B | median diff | A higher | p (Holm) |
| --- | --- | --- | --- | --- | --- | --- |
| COAD | cosine | B0_cohort_mean | M0_nnls | -0.0118 | 44% | 1 |
| COAD | mae_active | B0_cohort_mean | M0_nnls | +0.0054 | 52% | 0.87 |
| COAD | cosine | B0_cohort_mean | M1_elastic_net | -0.0008 | 50% | 1 |
| COAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0094 | 57% | 0.0098 |
| COAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.0276 | 35% | 4.7e-07 |
| COAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0211 | 64% | 8.3e-09 |
| COAD | cosine | B0_cohort_mean | M3_mlp | -0.0373 | 18% | 1.2e-43 |
| COAD | mae_active | B0_cohort_mean | M3_mlp | +0.0254 | 81% | 8.9e-42 |
| COAD | cosine | B0_cohort_mean | M3_mlp_broad | -0.0453 | 19% | 2.5e-37 |
| COAD | mae_active | B0_cohort_mean | M3_mlp_broad | +0.0333 | 81% | 2e-41 |
| COAD | cosine | M0_nnls | M1_elastic_net | +0.0068 | 53% | 1 |
| COAD | mae_active | M0_nnls | M1_elastic_net | +0.0066 | 55% | 0.022 |
| COAD | cosine | M0_nnls | M2_multinomial_em | -0.0171 | 35% | 8.3e-11 |
| COAD | mae_active | M0_nnls | M2_multinomial_em | +0.0151 | 66% | 2.3e-13 |
| COAD | cosine | M0_nnls | M3_mlp | -0.0281 | 35% | 1.7e-15 |
| COAD | mae_active | M0_nnls | M3_mlp | +0.0263 | 66% | 7.1e-16 |
| COAD | cosine | M0_nnls | M3_mlp_broad | -0.0298 | 28% | 6.4e-28 |
| COAD | mae_active | M0_nnls | M3_mlp_broad | +0.0297 | 74% | 1.9e-29 |
| COAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0238 | 31% | 7.6e-12 |
| COAD | mae_active | M1_elastic_net | M2_multinomial_em | +0.0124 | 64% | 8.1e-06 |
| COAD | cosine | M1_elastic_net | M3_mlp | -0.0375 | 26% | 1.3e-25 |
| COAD | mae_active | M1_elastic_net | M3_mlp | +0.0192 | 69% | 8.3e-17 |
| COAD | cosine | M1_elastic_net | M3_mlp_broad | -0.0419 | 18% | 1.5e-43 |
| COAD | mae_active | M1_elastic_net | M3_mlp_broad | +0.0253 | 78% | 1e-34 |
| COAD | cosine | M2_multinomial_em | M3_mlp | -0.0054 | 44% | 0.003 |
| COAD | mae_active | M2_multinomial_em | M3_mlp | +0.0036 | 53% | 0.094 |
| COAD | cosine | M2_multinomial_em | M3_mlp_broad | -0.0124 | 36% | 1.7e-12 |
| COAD | mae_active | M2_multinomial_em | M3_mlp_broad | +0.0136 | 63% | 1.6e-09 |
| COAD | cosine | M3_mlp | M3_mlp_broad | -0.0053 | 41% | 0.00027 |
| COAD | mae_active | M3_mlp | M3_mlp_broad | +0.0068 | 63% | 2.8e-07 |
| LUAD | cosine | B0_cohort_mean | M0_nnls | -0.0827 | 24% | 7e-29 |
| LUAD | mae_active | B0_cohort_mean | M0_nnls | +0.0445 | 72% | 5e-23 |
| LUAD | cosine | B0_cohort_mean | M1_elastic_net | -0.0885 | 17% | 8.9e-40 |
| LUAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0517 | 79% | 2e-35 |
| LUAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.0946 | 18% | 3.6e-39 |
| LUAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0562 | 77% | 2.5e-31 |
| LUAD | cosine | B0_cohort_mean | M3_mlp | -0.0928 | 9% | 1.2e-55 |
| LUAD | mae_active | B0_cohort_mean | M3_mlp | +0.0552 | 87% | 6.3e-52 |
| LUAD | cosine | B0_cohort_mean | M3_mlp_broad | -0.0931 | 12% | 1.3e-51 |
| LUAD | mae_active | B0_cohort_mean | M3_mlp_broad | +0.0577 | 85% | 3.3e-48 |
| LUAD | cosine | M0_nnls | M1_elastic_net | -0.0033 | 38% | 3.7e-08 |
| LUAD | mae_active | M0_nnls | M1_elastic_net | +0.0073 | 66% | 1.2e-08 |
| LUAD | cosine | M0_nnls | M2_multinomial_em | -0.0047 | 39% | 1.5e-05 |
| LUAD | mae_active | M0_nnls | M2_multinomial_em | +0.0080 | 61% | 1.5e-05 |
| LUAD | cosine | M0_nnls | M3_mlp | -0.0080 | 41% | 9.7e-06 |
| LUAD | mae_active | M0_nnls | M3_mlp | +0.0138 | 60% | 8e-06 |
| LUAD | cosine | M0_nnls | M3_mlp_broad | -0.0075 | 40% | 6.5e-08 |
| LUAD | mae_active | M0_nnls | M3_mlp_broad | +0.0133 | 62% | 9e-08 |
| LUAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0014 | 46% | 0.88 |
| LUAD | mae_active | M1_elastic_net | M2_multinomial_em | +0.0007 | 51% | 0.87 |
| LUAD | cosine | M1_elastic_net | M3_mlp | -0.0026 | 45% | 0.016 |
| LUAD | mae_active | M1_elastic_net | M3_mlp | +0.0039 | 54% | 0.094 |
| LUAD | cosine | M1_elastic_net | M3_mlp_broad | -0.0027 | 46% | 0.0011 |
| LUAD | mae_active | M1_elastic_net | M3_mlp_broad | +0.0039 | 54% | 0.0098 |
| LUAD | cosine | M2_multinomial_em | M3_mlp | +0.0005 | 52% | 1 |
| LUAD | mae_active | M2_multinomial_em | M3_mlp | +0.0018 | 51% | 0.55 |
| LUAD | cosine | M2_multinomial_em | M3_mlp_broad | -0.0007 | 49% | 0.44 |
| LUAD | mae_active | M2_multinomial_em | M3_mlp_broad | +0.0039 | 54% | 0.057 |
| LUAD | cosine | M3_mlp | M3_mlp_broad | -0.0008 | 47% | 0.88 |
| LUAD | mae_active | M3_mlp | M3_mlp_broad | +0.0018 | 53% | 0.48 |
| SKCM | cosine | B0_cohort_mean | M0_nnls | -0.0263 | 34% | 1.1e-11 |
| SKCM | mae_active | B0_cohort_mean | M0_nnls | +0.0250 | 64% | 1.2e-07 |
| SKCM | cosine | B0_cohort_mean | M1_elastic_net | -0.0272 | 32% | 5e-14 |
| SKCM | mae_active | B0_cohort_mean | M1_elastic_net | +0.0266 | 66% | 1.4e-09 |
| SKCM | cosine | B0_cohort_mean | M2_multinomial_em | -0.0415 | 24% | 9.2e-25 |
| SKCM | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0424 | 73% | 9.3e-22 |
| SKCM | cosine | B0_cohort_mean | M3_mlp | -0.0397 | 12% | 2.1e-51 |
| SKCM | mae_active | B0_cohort_mean | M3_mlp | +0.0428 | 88% | 1.2e-50 |
| SKCM | cosine | B0_cohort_mean | M3_mlp_broad | -0.0413 | 18% | 3.9e-44 |
| SKCM | mae_active | B0_cohort_mean | M3_mlp_broad | +0.0487 | 82% | 1e-45 |
| SKCM | cosine | M0_nnls | M1_elastic_net | -0.0017 | 38% | 3.9e-06 |
| SKCM | mae_active | M0_nnls | M1_elastic_net | +0.0015 | 56% | 0.00083 |
| SKCM | cosine | M0_nnls | M2_multinomial_em | -0.0103 | 28% | 2e-23 |
| SKCM | mae_active | M0_nnls | M2_multinomial_em | +0.0154 | 72% | 1.8e-20 |
| SKCM | cosine | M0_nnls | M3_mlp | -0.0117 | 34% | 8.4e-13 |
| SKCM | mae_active | M0_nnls | M3_mlp | +0.0203 | 66% | 8.4e-14 |
| SKCM | cosine | M0_nnls | M3_mlp_broad | -0.0143 | 29% | 7.6e-23 |
| SKCM | mae_active | M0_nnls | M3_mlp_broad | +0.0236 | 72% | 4.5e-20 |
| SKCM | cosine | M1_elastic_net | M2_multinomial_em | -0.0090 | 33% | 9.6e-17 |
| SKCM | mae_active | M1_elastic_net | M2_multinomial_em | +0.0148 | 68% | 1.2e-14 |
| SKCM | cosine | M1_elastic_net | M3_mlp | -0.0094 | 35% | 9.9e-11 |
| SKCM | mae_active | M1_elastic_net | M3_mlp | +0.0165 | 66% | 6.4e-12 |
| SKCM | cosine | M1_elastic_net | M3_mlp_broad | -0.0121 | 29% | 4e-21 |
| SKCM | mae_active | M1_elastic_net | M3_mlp_broad | +0.0229 | 70% | 3.7e-18 |
| SKCM | cosine | M2_multinomial_em | M3_mlp | -0.0004 | 49% | 1 |
| SKCM | mae_active | M2_multinomial_em | M3_mlp | +0.0020 | 52% | 0.55 |
| SKCM | cosine | M2_multinomial_em | M3_mlp_broad | -0.0024 | 44% | 0.012 |
| SKCM | mae_active | M2_multinomial_em | M3_mlp_broad | +0.0058 | 57% | 0.022 |
| SKCM | cosine | M3_mlp | M3_mlp_broad | -0.0006 | 47% | 0.42 |
| SKCM | mae_active | M3_mlp | M3_mlp_broad | +0.0018 | 55% | 0.2 |

## E5 smoking test (LUAD, real data)

{
 "n_smokers": 424,
 "n_never": 69,
 "median_smokers": 0.5032401893662389,
 "median_never": 0.028861518922623168,
 "mannwhitney_U": 24904.5,
 "mannwhitney_p": 3.678330070371453e-21,
 "spearman_rho_pack_years": 0.07544613192581809,
 "spearman_p": 0.1602093514198068,
 "n_pack_years": 348
}
