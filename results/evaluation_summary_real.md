# Evaluation summary (real-calibrated simulations and real TCGA tumors)

Written by run_evaluation.py.

## Median cosine by mutation count (simulated test tumors)

| cancer | method | 25 | 50 | 100 | 250 | 500 | 1000 | 2500 | 5000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COAD | B0_cohort_mean | 0.919 | 0.917 | 0.897 | 0.911 | 0.905 | 0.899 | 0.92 | 0.921 |
| COAD | M0_nnls | 0.814 | 0.859 | 0.895 | 0.944 | 0.941 | 0.962 | 0.961 | 0.945 |
| COAD | M1_elastic_net | 0.832 | 0.892 | 0.885 | 0.942 | 0.932 | 0.931 | 0.937 | 0.943 |
| COAD | M2_multinomial_em | 0.809 | 0.898 | 0.944 | 0.953 | 0.961 | 0.972 | 0.976 | 0.972 |
| COAD | M3_mlp | 0.938 | 0.953 | 0.957 | 0.966 | 0.96 | 0.968 | 0.966 | 0.944 |
| COAD | M3_mlp_broad | 0.913 | 0.949 | 0.962 | 0.973 | 0.966 | 0.978 | 0.979 | 0.971 |
| LUAD | B0_cohort_mean | 0.866 | 0.851 | 0.871 | 0.83 | 0.867 | 0.877 | 0.833 | 0.862 |
| LUAD | M0_nnls | 0.878 | 0.944 | 0.954 | 0.948 | 0.958 | 0.975 | 0.978 | 0.977 |
| LUAD | M1_elastic_net | 0.888 | 0.95 | 0.959 | 0.958 | 0.964 | 0.979 | 0.986 | 0.971 |
| LUAD | M2_multinomial_em | 0.877 | 0.939 | 0.952 | 0.962 | 0.973 | 0.983 | 0.983 | 0.982 |
| LUAD | M3_mlp | 0.933 | 0.952 | 0.945 | 0.959 | 0.966 | 0.986 | 0.981 | 0.972 |
| LUAD | M3_mlp_broad | 0.926 | 0.949 | 0.961 | 0.97 | 0.975 | 0.984 | 0.982 | 0.968 |
| SKCM | B0_cohort_mean | 0.945 | 0.925 | 0.939 | 0.913 | 0.92 | 0.935 | 0.945 | 0.919 |
| SKCM | M0_nnls | 0.86 | 0.942 | 0.958 | 0.972 | 0.982 | 0.978 | 0.983 | 0.973 |
| SKCM | M1_elastic_net | 0.887 | 0.952 | 0.966 | 0.971 | 0.982 | 0.98 | 0.982 | 0.972 |
| SKCM | M2_multinomial_em | 0.92 | 0.968 | 0.971 | 0.98 | 0.989 | 0.989 | 0.986 | 0.988 |
| SKCM | M3_mlp | 0.969 | 0.963 | 0.981 | 0.975 | 0.986 | 0.985 | 0.985 | 0.978 |
| SKCM | M3_mlp_broad | 0.958 | 0.975 | 0.98 | 0.981 | 0.99 | 0.99 | 0.988 | 0.985 |

## Paired Wilcoxon tests involving M3 (from results/method_tests_real.csv)

All mutation counts pooled; Holm-corrected by run_checklist.py over every pair of methods.
Median difference is method A minus method B: positive favours A for cosine, negative
favours A for MAE. Per-count tests are in the same file.

| cancer | metric | method_a | method_b | median_difference | a_higher_fraction | p_holm | significant |
| --- | --- | --- | --- | --- | --- | --- | --- |
| COAD | cosine | B0_cohort_mean | M3_mlp | -0.0373 | 0.1825 | 1.2e-43 | True |
| COAD | mae_active | B0_cohort_mean | M3_mlp | 0.0254 | 0.8075 | 8.9e-42 | True |
| COAD | cosine | B0_cohort_mean | M3_mlp_broad | -0.0453 | 0.1875 | 2.5e-37 | True |
| COAD | mae_active | B0_cohort_mean | M3_mlp_broad | 0.0333 | 0.81 | 2.0e-41 | True |
| COAD | cosine | M0_nnls | M3_mlp | -0.0281 | 0.3525 | 1.7e-15 | True |
| COAD | mae_active | M0_nnls | M3_mlp | 0.0263 | 0.665 | 7.1e-16 | True |
| COAD | cosine | M0_nnls | M3_mlp_broad | -0.0298 | 0.28 | 6.4e-28 | True |
| COAD | mae_active | M0_nnls | M3_mlp_broad | 0.0297 | 0.74 | 1.9e-29 | True |
| COAD | cosine | M1_elastic_net | M3_mlp | -0.0375 | 0.255 | 1.3e-25 | True |
| COAD | mae_active | M1_elastic_net | M3_mlp | 0.0192 | 0.6875 | 8.3e-17 | True |
| COAD | cosine | M1_elastic_net | M3_mlp_broad | -0.0419 | 0.175 | 1.5e-43 | True |
| COAD | mae_active | M1_elastic_net | M3_mlp_broad | 0.0253 | 0.78 | 1.0e-34 | True |
| COAD | cosine | M2_multinomial_em | M3_mlp | -0.0054 | 0.445 | 3.0e-03 | True |
| COAD | mae_active | M2_multinomial_em | M3_mlp | 0.0036 | 0.53 | 9.4e-02 | False |
| COAD | cosine | M2_multinomial_em | M3_mlp_broad | -0.0124 | 0.365 | 1.7e-12 | True |
| COAD | mae_active | M2_multinomial_em | M3_mlp_broad | 0.0136 | 0.6325 | 1.6e-09 | True |
| COAD | cosine | M3_mlp | M3_mlp_broad | -0.0053 | 0.4125 | 2.7e-04 | True |
| COAD | mae_active | M3_mlp | M3_mlp_broad | 0.0068 | 0.63 | 2.8e-07 | True |
| LUAD | cosine | B0_cohort_mean | M3_mlp | -0.0928 | 0.09 | 1.2e-55 | True |
| LUAD | mae_active | B0_cohort_mean | M3_mlp | 0.0552 | 0.8725 | 6.3e-52 | True |
| LUAD | cosine | B0_cohort_mean | M3_mlp_broad | -0.0931 | 0.1225 | 1.3e-51 | True |
| LUAD | mae_active | B0_cohort_mean | M3_mlp_broad | 0.0577 | 0.8475 | 3.3e-48 | True |
| LUAD | cosine | M0_nnls | M3_mlp | -0.008 | 0.41 | 9.7e-06 | True |
| LUAD | mae_active | M0_nnls | M3_mlp | 0.0138 | 0.6 | 8.0e-06 | True |
| LUAD | cosine | M0_nnls | M3_mlp_broad | -0.0075 | 0.395 | 6.5e-08 | True |
| LUAD | mae_active | M0_nnls | M3_mlp_broad | 0.0133 | 0.6225 | 9.0e-08 | True |
| LUAD | cosine | M1_elastic_net | M3_mlp | -0.0026 | 0.4525 | 1.6e-02 | True |
| LUAD | mae_active | M1_elastic_net | M3_mlp | 0.0039 | 0.5425 | 9.4e-02 | False |
| LUAD | cosine | M1_elastic_net | M3_mlp_broad | -0.0027 | 0.4625 | 1.1e-03 | True |
| LUAD | mae_active | M1_elastic_net | M3_mlp_broad | 0.0039 | 0.54 | 9.8e-03 | True |
| LUAD | cosine | M2_multinomial_em | M3_mlp | 0.0005 | 0.52 | 1.0e+00 | False |
| LUAD | mae_active | M2_multinomial_em | M3_mlp | 0.0018 | 0.51 | 5.5e-01 | False |
| LUAD | cosine | M2_multinomial_em | M3_mlp_broad | -0.0007 | 0.4875 | 4.4e-01 | False |
| LUAD | mae_active | M2_multinomial_em | M3_mlp_broad | 0.0039 | 0.535 | 5.7e-02 | False |
| LUAD | cosine | M3_mlp | M3_mlp_broad | -0.0008 | 0.47 | 8.8e-01 | False |
| LUAD | mae_active | M3_mlp | M3_mlp_broad | 0.0018 | 0.5275 | 4.8e-01 | False |
| SKCM | cosine | B0_cohort_mean | M3_mlp | -0.0397 | 0.1175 | 2.1e-51 | True |
| SKCM | mae_active | B0_cohort_mean | M3_mlp | 0.0428 | 0.875 | 1.2e-50 | True |
| SKCM | cosine | B0_cohort_mean | M3_mlp_broad | -0.0413 | 0.1775 | 3.9e-44 | True |
| SKCM | mae_active | B0_cohort_mean | M3_mlp_broad | 0.0487 | 0.8225 | 1.0e-45 | True |
| SKCM | cosine | M0_nnls | M3_mlp | -0.0117 | 0.345 | 8.4e-13 | True |
| SKCM | mae_active | M0_nnls | M3_mlp | 0.0203 | 0.665 | 8.4e-14 | True |
| SKCM | cosine | M0_nnls | M3_mlp_broad | -0.0143 | 0.2925 | 7.6e-23 | True |
| SKCM | mae_active | M0_nnls | M3_mlp_broad | 0.0236 | 0.715 | 4.5e-20 | True |
| SKCM | cosine | M1_elastic_net | M3_mlp | -0.0094 | 0.35 | 9.9e-11 | True |
| SKCM | mae_active | M1_elastic_net | M3_mlp | 0.0165 | 0.655 | 6.4e-12 | True |
| SKCM | cosine | M1_elastic_net | M3_mlp_broad | -0.0121 | 0.29 | 4.0e-21 | True |
| SKCM | mae_active | M1_elastic_net | M3_mlp_broad | 0.0229 | 0.7025 | 3.7e-18 | True |
| SKCM | cosine | M2_multinomial_em | M3_mlp | -0.0004 | 0.49 | 1.0e+00 | False |
| SKCM | mae_active | M2_multinomial_em | M3_mlp | 0.002 | 0.5175 | 5.5e-01 | False |
| SKCM | cosine | M2_multinomial_em | M3_mlp_broad | -0.0024 | 0.44 | 1.2e-02 | True |
| SKCM | mae_active | M2_multinomial_em | M3_mlp_broad | 0.0058 | 0.5675 | 2.2e-02 | True |
| SKCM | cosine | M3_mlp | M3_mlp_broad | -0.0006 | 0.4675 | 4.2e-01 | False |
| SKCM | mae_active | M3_mlp | M3_mlp_broad | 0.0018 | 0.545 | 2.0e-01 | False |

## Detection at 25 and 1,000 mutations

| cancer | method | n | roc_auc_macro | avg_precision_macro | hit_rate_at_fa5 | fdr |
| --- | --- | --- | --- | --- | --- | --- |
| LUAD | B0_cohort_mean | 25 | 0.5 | 0.673 | 0.243 | 0.02 |
| LUAD | B0_cohort_mean | 1000 | 0.5 | 0.61 | 0.268 | 0.02 |
| LUAD | M0_nnls | 25 | 0.837 | 0.912 | 0.698 | 0.066 |
| LUAD | M0_nnls | 1000 | 0.927 | 0.971 | 0.76 | 0.021 |
| LUAD | M1_elastic_net | 25 | 0.83 | 0.912 | 0.683 | 0.061 |
| LUAD | M1_elastic_net | 1000 | 0.945 | 0.977 | 0.749 | 0.021 |
| LUAD | M2_multinomial_em | 25 | 0.829 | 0.915 | 0.693 | 0.067 |
| LUAD | M2_multinomial_em | 1000 | 0.951 | 0.976 | 0.792 | 0.014 |
| LUAD | M3_mlp | 25 | 0.851 | 0.912 | 0.723 | 0.046 |
| LUAD | M3_mlp | 1000 | 0.939 | 0.971 | 0.738 | 0.022 |
| LUAD | M3_mlp_broad | 25 | 0.85 | 0.912 | 0.713 | 0.053 |
| LUAD | M3_mlp_broad | 1000 | 0.948 | 0.976 | 0.792 | 0.007 |
| SKCM | B0_cohort_mean | 25 | 0.5 | 0.288 | 0.581 | 0.0 |
| SKCM | B0_cohort_mean | 1000 | 0.5 | 0.37 | 0.615 | 0.01 |
| SKCM | M0_nnls | 25 | 0.667 | 0.569 | 0.779 | 0.107 |
| SKCM | M0_nnls | 1000 | 0.943 | 0.929 | 0.919 | 0.069 |
| SKCM | M1_elastic_net | 25 | 0.665 | 0.571 | 0.779 | 0.124 |
| SKCM | M1_elastic_net | 1000 | 0.948 | 0.927 | 0.913 | 0.075 |
| SKCM | M2_multinomial_em | 25 | 0.75 | 0.602 | 0.785 | 0.146 |
| SKCM | M2_multinomial_em | 1000 | 0.968 | 0.947 | 0.944 | 0.067 |
| SKCM | M3_mlp | 25 | 0.753 | 0.657 | 0.901 | 0.055 |
| SKCM | M3_mlp | 1000 | 0.952 | 0.926 | 0.907 | 0.076 |
| SKCM | M3_mlp_broad | 25 | 0.778 | 0.623 | 0.89 | 0.11 |
| SKCM | M3_mlp_broad | 1000 | 0.961 | 0.937 | 0.95 | 0.056 |
| COAD | B0_cohort_mean | 25 | 0.5 | 0.333 | 0.745 | 0.007 |
| COAD | B0_cohort_mean | 1000 | 0.5 | 0.403 | 0.754 | 0.04 |
| COAD | M0_nnls | 25 | 0.805 | 0.724 | 0.59 | 0.085 |
| COAD | M0_nnls | 1000 | 0.904 | 0.86 | 0.665 | 0.045 |
| COAD | M1_elastic_net | 25 | 0.826 | 0.745 | 0.575 | 0.087 |
| COAD | M1_elastic_net | 1000 | 0.911 | 0.858 | 0.644 | 0.039 |
| COAD | M2_multinomial_em | 25 | 0.829 | 0.767 | 0.65 | 0.116 |
| COAD | M2_multinomial_em | 1000 | 0.948 | 0.907 | 0.759 | 0.027 |
| COAD | M3_mlp | 25 | 0.823 | 0.726 | 0.8 | 0.048 |
| COAD | M3_mlp | 1000 | 0.914 | 0.866 | 0.827 | 0.06 |
| COAD | M3_mlp_broad | 25 | 0.816 | 0.752 | 0.785 | 0.071 |
| COAD | M3_mlp_broad | 1000 | 0.942 | 0.903 | 0.812 | 0.031 |

## Recall per signature (5% threshold), M2 vs M3

| cancer | signature | support | M2_multinomial_em | M3_mlp | M3_mlp_broad |
| --- | --- | --- | --- | --- | --- |
| COAD | SBS17b | 37 | 0.86 | 0.0 | 0.78 |
| COAD | SBS13 | 47 | 0.81 | 0.34 | 0.7 |
| COAD | SBS2 | 56 | 0.79 | 0.54 | 0.73 |
| COAD | SBS15 | 101 | 0.65 | 0.4 | 0.73 |
| COAD | SBS44 | 170 | 0.74 | 0.91 | 0.69 |
| COAD | SBS40a | 371 | 0.74 | 0.99 | 0.99 |
| COAD | SBS1 | 396 | 0.97 | 1.0 | 1.0 |
| COAD | SBS5 | 396 | 0.92 | 1.0 | 1.0 |
| LUAD | SBS1 | 111 | 0.83 | 0.81 | 0.85 |
| LUAD | SBS2 | 135 | 0.84 | 0.87 | 0.9 |
| LUAD | SBS13 | 146 | 0.9 | 0.92 | 0.9 |
| LUAD | SBS40a | 343 | 0.79 | 0.99 | 0.98 |
| LUAD | SBS5 | 345 | 0.81 | 0.97 | 0.97 |
| LUAD | SBS4 | 393 | 0.97 | 1.0 | 0.98 |
| SKCM | SBS38 | 34 | 0.76 | 0.41 | 0.71 |
| SKCM | SBS7c | 42 | 0.76 | 0.4 | 0.69 |
| SKCM | SBS7d | 69 | 0.71 | 0.71 | 0.67 |
| SKCM | SBS1 | 77 | 0.86 | 0.82 | 0.86 |
| SKCM | SBS5 | 310 | 0.85 | 0.96 | 0.91 |
| SKCM | SBS7b | 394 | 0.93 | 1.0 | 1.0 |
| SKCM | SBS7a | 399 | 0.97 | 1.0 | 1.0 |

## M3 robustness checks (median cosine; not the team's E3, which removes a listed signature)

| cancer | mutations | method | in-distribution | inverted prior | unseen signature |
| --- | --- | --- | --- | --- | --- |
| COAD | 25-100 | B0_cohort_mean | 0.913 | 0.106 | 0.903 |
| COAD | 25-100 | M0_nnls | 0.864 | 0.951 | 0.813 |
| COAD | 25-100 | M1_elastic_net | 0.874 | 0.961 | 0.846 |
| COAD | 25-100 | M2_multinomial_em | 0.899 | 0.969 | 0.787 |
| COAD | 25-100 | M3_mlp | 0.946 | 0.315 | 0.938 |
| COAD | 25-100 | M3_mlp_broad | 0.946 | 0.967 | 0.904 |
| COAD | 250-5000 | B0_cohort_mean | 0.912 | 0.109 | 0.904 |
| COAD | 250-5000 | M0_nnls | 0.953 | 0.986 | 0.852 |
| COAD | 250-5000 | M1_elastic_net | 0.938 | 0.988 | 0.916 |
| COAD | 250-5000 | M2_multinomial_em | 0.971 | 0.99 | 0.834 |
| COAD | 250-5000 | M3_mlp | 0.962 | 0.35 | 0.957 |
| COAD | 250-5000 | M3_mlp_broad | 0.974 | 0.985 | 0.944 |
| LUAD | 25-100 | B0_cohort_mean | 0.866 | 0.304 | 0.863 |
| LUAD | 25-100 | M0_nnls | 0.934 | 0.948 | 0.853 |
| LUAD | 25-100 | M1_elastic_net | 0.931 | 0.957 | 0.881 |
| LUAD | 25-100 | M2_multinomial_em | 0.935 | 0.976 | 0.801 |
| LUAD | 25-100 | M3_mlp | 0.943 | 0.902 | 0.92 |
| LUAD | 25-100 | M3_mlp_broad | 0.942 | 0.974 | 0.902 |
| LUAD | 250-5000 | B0_cohort_mean | 0.857 | 0.27 | 0.85 |
| LUAD | 250-5000 | M0_nnls | 0.974 | 0.985 | 0.903 |
| LUAD | 250-5000 | M1_elastic_net | 0.97 | 0.985 | 0.936 |
| LUAD | 250-5000 | M2_multinomial_em | 0.979 | 0.992 | 0.839 |
| LUAD | 250-5000 | M3_mlp | 0.972 | 0.96 | 0.945 |
| LUAD | 250-5000 | M3_mlp_broad | 0.976 | 0.988 | 0.921 |
| SKCM | 25-100 | B0_cohort_mean | 0.937 | 0.079 | 0.935 |
| SKCM | 25-100 | M0_nnls | 0.935 | 0.966 | 0.865 |
| SKCM | 25-100 | M1_elastic_net | 0.938 | 0.968 | 0.884 |
| SKCM | 25-100 | M2_multinomial_em | 0.96 | 0.975 | 0.866 |
| SKCM | 25-100 | M3_mlp | 0.972 | 0.415 | 0.95 |
| SKCM | 25-100 | M3_mlp_broad | 0.974 | 0.977 | 0.928 |
| SKCM | 250-5000 | B0_cohort_mean | 0.925 | 0.086 | 0.94 |
| SKCM | 250-5000 | M0_nnls | 0.977 | 0.988 | 0.931 |
| SKCM | 250-5000 | M1_elastic_net | 0.976 | 0.989 | 0.951 |
| SKCM | 250-5000 | M2_multinomial_em | 0.986 | 0.993 | 0.929 |
| SKCM | 250-5000 | M3_mlp | 0.984 | 0.43 | 0.965 |
| SKCM | 250-5000 | M3_mlp_broad | 0.986 | 0.989 | 0.963 |

## E5 smoking test, every method (real LUAD)

| method | patients | n_smokers | n_never | median_smokers | median_never | mannwhitney_p | spearman_rho_pack_years | spearman_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| M0_nnls | all patients | 424 | 69 | 0.5032 | 0.0289 | 3.7e-21 | 0.0754 | 1.6e-01 |
| M0_nnls | test patients | 83 | 11 | 0.4378 | 0.0482 | 3.1e-03 | 0.0915 | 4.5e-01 |
| M1_elastic_net | all patients | 424 | 69 | 0.4986 | 0.0452 | 4.2e-21 | 0.0706 | 1.9e-01 |
| M1_elastic_net | test patients | 83 | 11 | 0.4348 | 0.0692 | 3.2e-03 | 0.0862 | 4.8e-01 |
| M2_multinomial_em | all patients | 424 | 69 | 0.5143 | 0.0356 | 8.2e-22 | 0.0667 | 2.1e-01 |
| M2_multinomial_em | test patients | 83 | 11 | 0.4484 | 0.0225 | 2.1e-03 | 0.1046 | 3.9e-01 |
| M3_mlp | all patients | 424 | 69 | 0.4579 | 0.1559 | 1.7e-20 | 0.0356 | 5.1e-01 |
| M3_mlp | test patients | 83 | 11 | 0.4079 | 0.1894 | 3.0e-03 | 0.0668 | 5.9e-01 |
| M3_mlp_broad | all patients | 424 | 69 | 0.4495 | 0.0698 | 3.3e-21 | 0.0437 | 4.2e-01 |
| M3_mlp_broad | test patients | 83 | 11 | 0.3875 | 0.1027 | 2.1e-03 | 0.0629 | 6.1e-01 |

## Real tumors: reconstruction and agreement with M2 (median over tumors)

| cancer | method | reconstruction_cosine | cosine_vs_m2 |
| --- | --- | --- | --- |
| COAD | B0_cohort_mean | 0.913 | 0.872 |
| COAD | M0_nnls | 0.942 | 0.959 |
| COAD | M1_elastic_net | 0.938 | 0.954 |
| COAD | M2_multinomial_em | 0.938 | 1.0 |
| COAD | M3_mlp | 0.927 | 0.934 |
| COAD | M3_mlp_broad | 0.931 | 0.972 |
| LUAD | B0_cohort_mean | 0.811 | 0.887 |
| LUAD | M0_nnls | 0.903 | 0.982 |
| LUAD | M1_elastic_net | 0.903 | 0.979 |
| LUAD | M2_multinomial_em | 0.901 | 1.0 |
| LUAD | M3_mlp | 0.898 | 0.984 |
| LUAD | M3_mlp_broad | 0.896 | 0.983 |
| SKCM | B0_cohort_mean | 0.982 | 0.966 |
| SKCM | M0_nnls | 0.989 | 0.995 |
| SKCM | M1_elastic_net | 0.989 | 0.996 |
| SKCM | M2_multinomial_em | 0.988 | 1.0 |
| SKCM | M3_mlp | 0.986 | 0.99 |
| SKCM | M3_mlp_broad | 0.987 | 0.994 |

## Figures

- figures/eval_fig1_accuracy_vs_mutations_real.png
- figures/eval_fig2_detection_real.png
- figures/eval_fig3_recall_vs_prevalence_real.png
- figures/eval_fig4_robustness_real.png
- figures/eval_fig5_confusion_real.png
- figures/eval_fig6_e5_smoking_real.png
- figures/eval_fig7_m3_training_real.png
- figures/eval_fig8_real_reconstruction_real.png
