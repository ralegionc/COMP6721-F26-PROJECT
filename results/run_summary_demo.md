# Run summary (demo data, 2026-10-30)

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
| LUAD | M1_elastic_net | 0.943 | 0.098 | 0.800 | 0.888 |
| LUAD | M2_multinomial_em | 0.954 | 0.093 | 0.824 | 0.908 |
| SKCM | B0_cohort_mean | 0.694 | 0.126 | 0.745 | 0.598 |
| SKCM | M0_nnls | 0.944 | 0.080 | 0.816 | 0.901 |
| SKCM | M1_elastic_net | 0.944 | 0.080 | 0.816 | 0.901 |
| SKCM | M2_multinomial_em | 0.961 | 0.072 | 0.829 | 0.921 |
| COAD | B0_cohort_mean | 0.708 | 0.102 | 0.751 | 0.604 |
| COAD | M0_nnls | 0.916 | 0.082 | 0.778 | 0.859 |
| COAD | M1_elastic_net | 0.916 | 0.082 | 0.777 | 0.859 |
| COAD | M2_multinomial_em | 0.926 | 0.077 | 0.784 | 0.879 |

## Breakdown point N* (smallest N with median cosine >= 0.90)

| method | COAD | LUAD | SKCM |
| --- | --- | --- | --- |
| B0_cohort_mean | not reached | not reached | not reached |
| M0_nnls | 250 | 50 | 100 |
| M1_elastic_net | 250 | 50 | 100 |
| M2_multinomial_em | 100 | 50 | 50 |

M1 penalty chosen on validation data: {'LUAD': 1e-07, 'SKCM': 1e-07, 'COAD': 1e-06}

## E5 smoking test (LUAD, real data)

{
 "n_smokers": 36,
 "n_never": 12,
 "median_smokers": 0.08610760046461383,
 "median_never": 0.10751587087840145,
 "mannwhitney_U": 211.0,
 "mannwhitney_p": 0.5520991349176474
}
