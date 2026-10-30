# Run summary (real data, 2026-10-22)

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
| LUAD | M1_elastic_net | 0.959 | 0.089 | 0.870 | 0.940 |
| LUAD | M2_multinomial_em | 0.966 | 0.080 | 0.873 | 0.946 |
| SKCM | B0_cohort_mean | 0.931 | 0.111 | 0.409 | 0.473 |
| SKCM | M0_nnls | 0.966 | 0.090 | 0.762 | 0.840 |
| SKCM | M1_elastic_net | 0.966 | 0.090 | 0.762 | 0.840 |
| SKCM | M2_multinomial_em | 0.980 | 0.073 | 0.800 | 0.864 |
| COAD | B0_cohort_mean | 0.913 | 0.100 | 0.444 | 0.492 |
| COAD | M0_nnls | 0.930 | 0.101 | 0.749 | 0.822 |
| COAD | M1_elastic_net | 0.930 | 0.101 | 0.749 | 0.822 |
| COAD | M2_multinomial_em | 0.955 | 0.082 | 0.779 | 0.844 |

## Breakdown point N* (smallest N with median cosine >= 0.90)

| method | COAD | LUAD | SKCM |
| --- | --- | --- | --- |
| B0_cohort_mean | 2500 | not reached | 25 |
| M0_nnls | 250 | 50 | 50 |
| M1_elastic_net | 250 | 50 | 50 |
| M2_multinomial_em | 100 | 50 | 25 |

M1 penalty chosen on validation data: {'LUAD': 1e-07, 'SKCM': 1e-07, 'COAD': 1e-07}

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
