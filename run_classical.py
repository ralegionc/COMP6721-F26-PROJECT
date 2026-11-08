
import pandas as pd
import os

from classical_methods.nnls import estimate_exposures_nnls
from classical_methods.regularized_nnls import (
    estimate_exposures_regularized_nnls
)
from classical_methods.evaluation import evaluate_exposures


# ---------------------------------------------------
# 1. Load simulated data
# ---------------------------------------------------

signature_matrix = pd.read_csv(
    "data/processed/reference_signatures.csv",
    index_col=0
)

profiles = pd.read_csv(
    "data/processed/simulated/simulated_sbs96.csv",
    index_col=0
)

true_exposures = pd.read_csv(
    "data/processed/simulated/simulated_exposures.csv",
    index_col=0
)


# ---------------------------------------------------
# 2. Run NNLS
# ---------------------------------------------------

nnls_estimates = estimate_exposures_nnls(
    profiles,
    signature_matrix
)


# ---------------------------------------------------
# 3. Run regularized NNLS
# ---------------------------------------------------

regularized_estimates = estimate_exposures_regularized_nnls(
    profiles,
    signature_matrix,
    alpha=0.01
)


# ---------------------------------------------------
# 4. Evaluate both methods
# ---------------------------------------------------

nnls_summary, nnls_per_tumor = evaluate_exposures(
    true_exposures,
    nnls_estimates,
    profiles=profiles,
    signature_matrix=signature_matrix
)

regularized_summary, regularized_per_tumor = evaluate_exposures(
    true_exposures,
    regularized_estimates,
    profiles=profiles,
    signature_matrix=signature_matrix
)


# ---------------------------------------------------
# 5. Compare results
# ---------------------------------------------------

comparison = pd.DataFrame({
    "NNLS": nnls_summary,
    "Regularized NNLS": regularized_summary
})

print("\nClassical Method Comparison:")
print(comparison.round(4))


# ---------------------------------------------------
# 6. Save results
# ---------------------------------------------------

output_dir = "data/processed/classical"


os.makedirs(output_dir, exist_ok=True)

nnls_estimates.to_csv(
    f"{output_dir}/nnls_exposures.csv"
)

regularized_estimates.to_csv(
    f"{output_dir}/regularized_nnls_exposures.csv"
)

comparison.to_csv(
    f"{output_dir}/method_comparison.csv"
)

nnls_per_tumor.to_csv(
    f"{output_dir}/nnls_per_tumor.csv"
)

regularized_per_tumor.to_csv(
    f"{output_dir}/regularized_nnls_per_tumor.csv"
)

print("\nResults saved successfully!")
