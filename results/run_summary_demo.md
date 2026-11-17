# Run summary (demo data, 2026-11-16)

## Checklist
- Q1: Exposure vector per tumor: non-negative shares over the cancer type's candidate COSMIC SBS signatures (artifacts excluded), summing to 1. Known only for simulated tumors; real tumors give directional checks only. Subsets: LUAD=6, SKCM=7, COAD=8 (demo signatures)
- Q2: Split by patient (train=108, val=36, test=36); no overlap; multi-sample patients kept together; simulator calibrated on training patients only.
- Q3: Presence: signature weight >= 5%. Smoking: current or reformed smoker (negative: lifelong non-smoker; unknown excluded). LUAD labels: 36 positive, 12 negative, 12 excluded.
- Q4: Max 'absent' rate 40% (SKCM); absent labels are not the majority in this run. Reported: per-signature F1/PR-AUC, cosine, MAE vs B0 (M0 macro F1 0.78-0.82).
- Q5: Split + 9 simulated sets saved once (seed 6721, SHA-256); all 6 methods scored on the same files.
- Q6: results/, figures/, splits/ committed; data/raw and data/sims gitignored (manifest + seeds).

## Test-set results (simulated tumors, all mutation counts pooled)

| cancer | method | median cosine | MAE (active) | macro F1 | PR-AUC |
| --- | --- | --- | --- | --- | --- |
| LUAD | B0_cohort_mean | 0.687 | 0.152 | 0.750 | 0.603 |
| LUAD | M0_nnls | 0.943 | 0.098 | 0.800 | 0.888 |
| LUAD | M1_elastic_net | 0.946 | 0.091 | 0.811 | 0.891 |
| LUAD | M2_multinomial_em | 0.954 | 0.093 | 0.824 | 0.908 |
| LUAD | M3_mlp | 0.945 | 0.083 | 0.816 | 0.906 |
| LUAD | M3_mlp_broad | 0.946 | 0.083 | 0.818 | 0.909 |
| SKCM | B0_cohort_mean | 0.694 | 0.126 | 0.745 | 0.598 |
| SKCM | M0_nnls | 0.944 | 0.080 | 0.816 | 0.901 |
| SKCM | M1_elastic_net | 0.947 | 0.076 | 0.819 | 0.904 |
| SKCM | M2_multinomial_em | 0.961 | 0.072 | 0.829 | 0.921 |
| SKCM | M3_mlp | 0.953 | 0.065 | 0.837 | 0.920 |
| SKCM | M3_mlp_broad | 0.959 | 0.063 | 0.840 | 0.927 |
| COAD | B0_cohort_mean | 0.708 | 0.102 | 0.751 | 0.604 |
| COAD | M0_nnls | 0.916 | 0.082 | 0.778 | 0.859 |
| COAD | M1_elastic_net | 0.907 | 0.074 | 0.786 | 0.868 |
| COAD | M2_multinomial_em | 0.926 | 0.077 | 0.784 | 0.879 |
| COAD | M3_mlp | 0.916 | 0.066 | 0.797 | 0.877 |
| COAD | M3_mlp_broad | 0.916 | 0.066 | 0.802 | 0.881 |

## Breakdown point N* (smallest N with median cosine >= 0.90)

| method | COAD | LUAD | SKCM |
| --- | --- | --- | --- |
| B0_cohort_mean | not reached | not reached | not reached |
| M0_nnls | 250 | 50 | 100 |
| M1_elastic_net | 1000 | 100 | 50 |
| M2_multinomial_em | 100 | 50 | 50 |
| M3_mlp | 100 | 100 | 50 |
| M3_mlp_broad | 100 | 100 | 50 |

M1 penalty chosen on validation data: {'LUAD': {'alpha': 1e-05, 'l1_ratio': 0.1}, 'SKCM': {'alpha': 3e-05, 'l1_ratio': 0.1}, 'COAD': {'alpha': 0.0001, 'l1_ratio': 0.1}}

## Paired Wilcoxon tests (all mutation counts pooled, Holm-corrected)

Median difference is method A minus method B; "A higher" is the share of tumors where A scores higher. Higher cosine and lower MAE are better.

| cancer | metric | A | B | median diff | A higher | p (Holm) |
| --- | --- | --- | --- | --- | --- | --- |
| COAD | cosine | B0_cohort_mean | M0_nnls | -0.1805 | 18% | 1.1e-41 |
| COAD | mae_active | B0_cohort_mean | M0_nnls | +0.0244 | 68% | 3.5e-16 |
| COAD | cosine | B0_cohort_mean | M1_elastic_net | -0.1752 | 10% | 2.6e-57 |
| COAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0287 | 78% | 1.2e-34 |
| COAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.1917 | 15% | 2.3e-43 |
| COAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0319 | 71% | 2.5e-20 |
| COAD | cosine | B0_cohort_mean | M3_mlp | -0.1844 | 6% | 8.4e-63 |
| COAD | mae_active | B0_cohort_mean | M3_mlp | +0.0341 | 88% | 6.4e-54 |
| COAD | cosine | B0_cohort_mean | M3_mlp_broad | -0.1804 | 6% | 1.5e-61 |
| COAD | mae_active | B0_cohort_mean | M3_mlp_broad | +0.0344 | 87% | 6.7e-52 |
| COAD | cosine | M0_nnls | M1_elastic_net | -0.0042 | 42% | 0.023 |
| COAD | mae_active | M0_nnls | M1_elastic_net | +0.0065 | 62% | 8.8e-07 |
| COAD | cosine | M0_nnls | M2_multinomial_em | -0.0069 | 39% | 5.4e-05 |
| COAD | mae_active | M0_nnls | M2_multinomial_em | +0.0048 | 57% | 0.00043 |
| COAD | cosine | M0_nnls | M3_mlp | -0.0103 | 39% | 1.1e-05 |
| COAD | mae_active | M0_nnls | M3_mlp | +0.0138 | 69% | 9.4e-14 |
| COAD | cosine | M0_nnls | M3_mlp_broad | -0.0111 | 40% | 1.9e-06 |
| COAD | mae_active | M0_nnls | M3_mlp_broad | +0.0127 | 69% | 1.8e-15 |
| COAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0020 | 49% | 1 |
| COAD | mae_active | M1_elastic_net | M2_multinomial_em | -0.0029 | 48% | 1 |
| COAD | cosine | M1_elastic_net | M3_mlp | -0.0070 | 40% | 2.8e-06 |
| COAD | mae_active | M1_elastic_net | M3_mlp | +0.0062 | 64% | 8.3e-11 |
| COAD | cosine | M1_elastic_net | M3_mlp_broad | -0.0070 | 41% | 4.1e-05 |
| COAD | mae_active | M1_elastic_net | M3_mlp_broad | +0.0064 | 64% | 1.1e-10 |
| COAD | cosine | M2_multinomial_em | M3_mlp | -0.0028 | 48% | 0.23 |
| COAD | mae_active | M2_multinomial_em | M3_mlp | +0.0086 | 60% | 4.8e-06 |
| COAD | cosine | M2_multinomial_em | M3_mlp_broad | -0.0023 | 48% | 0.23 |
| COAD | mae_active | M2_multinomial_em | M3_mlp_broad | +0.0083 | 64% | 3.6e-07 |
| COAD | cosine | M3_mlp | M3_mlp_broad | +0.0004 | 51% | 1 |
| COAD | mae_active | M3_mlp | M3_mlp_broad | +0.0004 | 51% | 1 |
| LUAD | cosine | B0_cohort_mean | M0_nnls | -0.2222 | 10% | 1.6e-49 |
| LUAD | mae_active | B0_cohort_mean | M0_nnls | +0.0543 | 80% | 1.2e-33 |
| LUAD | cosine | B0_cohort_mean | M1_elastic_net | -0.2246 | 9% | 2.4e-53 |
| LUAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0594 | 82% | 5.2e-41 |
| LUAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.2344 | 10% | 2e-50 |
| LUAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0630 | 78% | 1e-35 |
| LUAD | cosine | B0_cohort_mean | M3_mlp | -0.2290 | 6% | 1.3e-59 |
| LUAD | mae_active | B0_cohort_mean | M3_mlp | +0.0666 | 88% | 1.4e-53 |
| LUAD | cosine | B0_cohort_mean | M3_mlp_broad | -0.2333 | 6% | 2e-61 |
| LUAD | mae_active | B0_cohort_mean | M3_mlp_broad | +0.0715 | 90% | 6.9e-56 |
| LUAD | cosine | M0_nnls | M1_elastic_net | -0.0014 | 42% | 0.01 |
| LUAD | mae_active | M0_nnls | M1_elastic_net | +0.0038 | 61% | 1.1e-05 |
| LUAD | cosine | M0_nnls | M2_multinomial_em | -0.0037 | 39% | 0.0005 |
| LUAD | mae_active | M0_nnls | M2_multinomial_em | +0.0059 | 61% | 0.0009 |
| LUAD | cosine | M0_nnls | M3_mlp | -0.0060 | 44% | 0.027 |
| LUAD | mae_active | M0_nnls | M3_mlp | +0.0091 | 58% | 0.00011 |
| LUAD | cosine | M0_nnls | M3_mlp_broad | -0.0069 | 42% | 0.00099 |
| LUAD | mae_active | M0_nnls | M3_mlp_broad | +0.0098 | 60% | 8.7e-05 |
| LUAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0029 | 42% | 0.23 |
| LUAD | mae_active | M1_elastic_net | M2_multinomial_em | +0.0001 | 50% | 1 |
| LUAD | cosine | M1_elastic_net | M3_mlp | -0.0019 | 47% | 0.34 |
| LUAD | mae_active | M1_elastic_net | M3_mlp | +0.0060 | 56% | 0.02 |
| LUAD | cosine | M1_elastic_net | M3_mlp_broad | -0.0036 | 45% | 0.014 |
| LUAD | mae_active | M1_elastic_net | M3_mlp_broad | +0.0061 | 57% | 0.019 |
| LUAD | cosine | M2_multinomial_em | M3_mlp | +0.0015 | 52% | 1 |
| LUAD | mae_active | M2_multinomial_em | M3_mlp | +0.0055 | 57% | 0.071 |
| LUAD | cosine | M2_multinomial_em | M3_mlp_broad | -0.0019 | 48% | 1 |
| LUAD | mae_active | M2_multinomial_em | M3_mlp_broad | +0.0044 | 57% | 0.075 |
| LUAD | cosine | M3_mlp | M3_mlp_broad | -0.0020 | 45% | 0.41 |
| LUAD | mae_active | M3_mlp | M3_mlp_broad | -0.0008 | 49% | 1 |
| SKCM | cosine | B0_cohort_mean | M0_nnls | -0.2336 | 9% | 1.5e-53 |
| SKCM | mae_active | B0_cohort_mean | M0_nnls | +0.0528 | 80% | 4.2e-37 |
| SKCM | cosine | B0_cohort_mean | M1_elastic_net | -0.2338 | 7% | 5e-58 |
| SKCM | mae_active | B0_cohort_mean | M1_elastic_net | +0.0551 | 85% | 1.9e-45 |
| SKCM | cosine | B0_cohort_mean | M2_multinomial_em | -0.2460 | 7% | 1e-57 |
| SKCM | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0624 | 82% | 1.9e-43 |
| SKCM | cosine | B0_cohort_mean | M3_mlp | -0.2453 | 4% | 1.5e-63 |
| SKCM | mae_active | B0_cohort_mean | M3_mlp | +0.0662 | 92% | 1.6e-59 |
| SKCM | cosine | B0_cohort_mean | M3_mlp_broad | -0.2468 | 3% | 1.7e-64 |
| SKCM | mae_active | B0_cohort_mean | M3_mlp_broad | +0.0656 | 94% | 5.7e-61 |
| SKCM | cosine | M0_nnls | M1_elastic_net | -0.0019 | 40% | 0.0012 |
| SKCM | mae_active | M0_nnls | M1_elastic_net | +0.0023 | 62% | 0.0011 |
| SKCM | cosine | M0_nnls | M2_multinomial_em | -0.0085 | 32% | 1.5e-13 |
| SKCM | mae_active | M0_nnls | M2_multinomial_em | +0.0079 | 64% | 2.6e-09 |
| SKCM | cosine | M0_nnls | M3_mlp | -0.0120 | 34% | 4.1e-10 |
| SKCM | mae_active | M0_nnls | M3_mlp | +0.0117 | 64% | 4.9e-11 |
| SKCM | cosine | M0_nnls | M3_mlp_broad | -0.0145 | 32% | 6.3e-15 |
| SKCM | mae_active | M0_nnls | M3_mlp_broad | +0.0130 | 68% | 5.6e-15 |
| SKCM | cosine | M1_elastic_net | M2_multinomial_em | -0.0079 | 36% | 9.7e-07 |
| SKCM | mae_active | M1_elastic_net | M2_multinomial_em | +0.0066 | 59% | 0.00027 |
| SKCM | cosine | M1_elastic_net | M3_mlp | -0.0095 | 36% | 1.5e-08 |
| SKCM | mae_active | M1_elastic_net | M3_mlp | +0.0087 | 64% | 8.5e-10 |
| SKCM | cosine | M1_elastic_net | M3_mlp_broad | -0.0111 | 34% | 1.5e-14 |
| SKCM | mae_active | M1_elastic_net | M3_mlp_broad | +0.0105 | 67% | 7e-15 |
| SKCM | cosine | M2_multinomial_em | M3_mlp | +0.0026 | 55% | 1 |
| SKCM | mae_active | M2_multinomial_em | M3_mlp | +0.0026 | 53% | 0.58 |
| SKCM | cosine | M2_multinomial_em | M3_mlp_broad | +0.0010 | 52% | 1 |
| SKCM | mae_active | M2_multinomial_em | M3_mlp_broad | +0.0041 | 56% | 0.041 |
| SKCM | cosine | M3_mlp | M3_mlp_broad | -0.0014 | 46% | 0.043 |
| SKCM | mae_active | M3_mlp | M3_mlp_broad | +0.0020 | 57% | 0.071 |

## E5 smoking test (LUAD, real data)

{
 "n_smokers": 36,
 "n_never": 12,
 "median_smokers": 0.08610760046461384,
 "median_never": 0.1075158708784014,
 "mannwhitney_U": 211.0,
 "mannwhitney_p": 0.5520991349176474
}
