import os

import pandas as pd

from statistical_tests.statistical_analysis import analyze_by_mutation_count


# ---------------------------------------------------
# 1. Load per-tumor metrics and simulation metadata
# ---------------------------------------------------

classical_dir = "data/processed/classical"

nnls_metrics = pd.read_csv(
    f"{classical_dir}/nnls_per_tumor.csv",
    index_col=0
)

regularized_metrics = pd.read_csv(
    f"{classical_dir}/regularized_nnls_per_tumor.csv",
    index_col=0
)

metadata = pd.read_csv(
    "data/processed/simulated/simulation_metadata.csv"
)


# ---------------------------------------------------
# 2. Compare NNLS and regularized NNLS
# ---------------------------------------------------

summary, tests = analyze_by_mutation_count(
    nnls_metrics,
    regularized_metrics,
    metadata,
    metrics=("MAE", "RMSE", "Support F1")
)


# ---------------------------------------------------
# 3. Print results
# ---------------------------------------------------

print("\nDescriptive statistics:")
print(summary.round(4).to_string(index=False))

print("\nPaired Wilcoxon tests (Holm-adjusted per metric):")
print(tests.round(4).to_string(index=False))


# ---------------------------------------------------
# 4. Save results
# ---------------------------------------------------

output_dir = "data/processed/statistics"
os.makedirs(output_dir, exist_ok=True)

summary.to_csv(f"{output_dir}/summary_by_mutation_count.csv", index=False)
tests.to_csv(f"{output_dir}/wilcoxon_tests.csv", index=False)

print(f"\nResults saved to {output_dir}/")
