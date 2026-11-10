# Run summary (demo data, 2026-11-10)

## Checklist
- Q1: Exposure vector per tumor: non-negative shares over the cancer type's candidate COSMIC SBS signatures (artifacts excluded), summing to 1. Known only for simulated tumors; real tumors give directional checks only. Subsets: LUAD=6, SKCM=7, COAD=8 (demo signatures)
- Q2: Split by patient (train=108, val=36, test=36); no overlap; multi-sample patients kept together; simulator calibrated on training patients only.
- Q3: Presence: signature weight >= 5%. Smoking: current or reformed smoker (negative: lifelong non-smoker; unknown excluded). LUAD labels: 36 positive, 12 negative, 12 excluded.
- Q4: Max 'absent' rate 40% (SKCM); absent labels are not the majority in this run. Reported: per-signature F1/PR-AUC, cosine, MAE vs B0 (M0 macro F1 0.78-0.82).
- Q5: Split + 9 simulated sets saved once (seed 6721, SHA-256); all 4 methods scored on the same files.
- Q6: results/, figures/, splits/ committed; data/raw and data/sims gitignored (manifest + seeds).

## Test-set results (simulated tumors, all mutation counts pooled)

| cancer | method | median cosine | MAE (active) | macro F1 | PR-AUC |
| --- | --- | --- | --- | --- | --- |
| LUAD | B0_cohort_mean | 0.687 | 0.152 | 0.750 | 0.603 |
| LUAD | M0_nnls | 0.943 | 0.098 | 0.800 | 0.888 |
| LUAD | M1_elastic_net | 0.946 | 0.091 | 0.811 | 0.891 |
| LUAD | M2_multinomial_em | 0.954 | 0.093 | 0.824 | 0.908 |
| SKCM | B0_cohort_mean | 0.694 | 0.126 | 0.745 | 0.598 |
| SKCM | M0_nnls | 0.944 | 0.080 | 0.816 | 0.901 |
| SKCM | M1_elastic_net | 0.947 | 0.076 | 0.819 | 0.904 |
| SKCM | M2_multinomial_em | 0.961 | 0.072 | 0.829 | 0.921 |
| COAD | B0_cohort_mean | 0.708 | 0.102 | 0.751 | 0.604 |
| COAD | M0_nnls | 0.916 | 0.082 | 0.778 | 0.859 |
| COAD | M1_elastic_net | 0.907 | 0.074 | 0.786 | 0.868 |
| COAD | M2_multinomial_em | 0.926 | 0.077 | 0.784 | 0.879 |

## Breakdown point N* (smallest N with median cosine >= 0.90)

| method | COAD | LUAD | SKCM |
| --- | --- | --- | --- |
| B0_cohort_mean | not reached | not reached | not reached |
| M0_nnls | 250 | 50 | 100 |
| M1_elastic_net | 1000 | 100 | 50 |
| M2_multinomial_em | 100 | 50 | 50 |

M1 penalty chosen on validation data: {'LUAD': {'alpha': 1e-05, 'l1_ratio': 0.1}, 'SKCM': {'alpha': 3e-05, 'l1_ratio': 0.1}, 'COAD': {'alpha': 0.0001, 'l1_ratio': 0.1}}

## Paired Wilcoxon tests (all mutation counts pooled, Holm-corrected)

Median difference is method A minus method B; "A higher" is the share of tumors where A scores higher. Higher cosine and lower MAE are better.

| cancer | metric | A | B | median diff | A higher | p (Holm) |
| --- | --- | --- | --- | --- | --- | --- |
| COAD | cosine | B0_cohort_mean | M0_nnls | -0.1805 | 18% | 3.6e-42 |
| COAD | mae_active | B0_cohort_mean | M0_nnls | +0.0244 | 68% | 1.1e-16 |
| COAD | cosine | B0_cohort_mean | M1_elastic_net | -0.1752 | 10% | 1.1e-57 |
| COAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0287 | 78% | 4.6e-35 |
| COAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.1917 | 15% | 7.9e-44 |
| COAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0319 | 71% | 8.8e-21 |
| COAD | cosine | M0_nnls | M1_elastic_net | -0.0042 | 42% | 0.0049 |
| COAD | mae_active | M0_nnls | M1_elastic_net | +0.0065 | 62% | 3.5e-07 |
| COAD | cosine | M0_nnls | M2_multinomial_em | -0.0069 | 39% | 1.9e-05 |
| COAD | mae_active | M0_nnls | M2_multinomial_em | +0.0048 | 57% | 0.00015 |
| COAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0020 | 49% | 0.6 |
| COAD | mae_active | M1_elastic_net | M2_multinomial_em | -0.0029 | 48% | 0.61 |
| LUAD | cosine | B0_cohort_mean | M0_nnls | -0.2222 | 10% | 6e-50 |
| LUAD | mae_active | B0_cohort_mean | M0_nnls | +0.0543 | 80% | 4.5e-34 |
| LUAD | cosine | B0_cohort_mean | M1_elastic_net | -0.2246 | 9% | 9.7e-54 |
| LUAD | mae_active | B0_cohort_mean | M1_elastic_net | +0.0594 | 82% | 2.2e-41 |
| LUAD | cosine | B0_cohort_mean | M2_multinomial_em | -0.2344 | 10% | 7.7e-51 |
| LUAD | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0630 | 78% | 4.2e-36 |
| LUAD | cosine | M0_nnls | M1_elastic_net | -0.0014 | 42% | 0.0026 |
| LUAD | mae_active | M0_nnls | M1_elastic_net | +0.0038 | 61% | 4.4e-06 |
| LUAD | cosine | M0_nnls | M2_multinomial_em | -0.0037 | 39% | 0.00016 |
| LUAD | mae_active | M0_nnls | M2_multinomial_em | +0.0059 | 61% | 0.00028 |
| LUAD | cosine | M1_elastic_net | M2_multinomial_em | -0.0029 | 42% | 0.043 |
| LUAD | mae_active | M1_elastic_net | M2_multinomial_em | +0.0001 | 50% | 0.85 |
| SKCM | cosine | B0_cohort_mean | M0_nnls | -0.2336 | 9% | 6.1e-54 |
| SKCM | mae_active | B0_cohort_mean | M0_nnls | +0.0528 | 80% | 1.7e-37 |
| SKCM | cosine | B0_cohort_mean | M1_elastic_net | -0.2338 | 7% | 2.3e-58 |
| SKCM | mae_active | B0_cohort_mean | M1_elastic_net | +0.0551 | 85% | 8.7e-46 |
| SKCM | cosine | B0_cohort_mean | M2_multinomial_em | -0.2460 | 7% | 4.5e-58 |
| SKCM | mae_active | B0_cohort_mean | M2_multinomial_em | +0.0624 | 82% | 8.3e-44 |
| SKCM | cosine | M0_nnls | M1_elastic_net | -0.0019 | 40% | 0.00036 |
| SKCM | mae_active | M0_nnls | M1_elastic_net | +0.0023 | 62% | 0.00028 |
| SKCM | cosine | M0_nnls | M2_multinomial_em | -0.0085 | 32% | 5e-14 |
| SKCM | mae_active | M0_nnls | M2_multinomial_em | +0.0079 | 64% | 1.1e-09 |
| SKCM | cosine | M1_elastic_net | M2_multinomial_em | -0.0079 | 36% | 3.1e-07 |
| SKCM | mae_active | M1_elastic_net | M2_multinomial_em | +0.0066 | 59% | 0.00011 |

## E5 smoking test (LUAD, real data)

{
 "n_smokers": 36,
 "n_never": 12,
 "median_smokers": 0.08610760046461383,
 "median_never": 0.10751587087840145,
 "mannwhitney_U": 211.0,
 "mannwhitney_p": 0.5520991349176474
}
