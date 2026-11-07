
import pandas as pd

from statistical_tests.hypothesis_tests import (
    paired_wilcoxon,
    holm_correction
)


def analyze_by_mutation_count(
    nnls_metrics,
    regularized_metrics,
    metadata,
    metrics=("MAE", "RMSE", "Support F1")
):
    """
    Compare NNLS and regularized NNLS across mutation burdens.

    Parameters
    ----------
    nnls_metrics : pandas.DataFrame
        Per-tumor metrics produced by evaluation.py.

    regularized_metrics : pandas.DataFrame
        Per-tumor metrics produced by evaluation.py.

    metadata : pandas.DataFrame
        Simulation metadata containing tumor_id and mutation_count.

    metrics : tuple[str]
        Metrics to summarize and compare.

    Returns
    -------
    summary : pandas.DataFrame
        Descriptive statistics for each method and mutation burden.

    tests : pandas.DataFrame
        Paired Wilcoxon results and Holm-adjusted p-values.
    """

    if not nnls_metrics.index.is_unique:
        raise ValueError("NNLS tumor IDs must be unique.")

    if not regularized_metrics.index.is_unique:
        raise ValueError("Regularized NNLS tumor IDs must be unique.")

    if set(nnls_metrics.index) != set(regularized_metrics.index):
        raise ValueError("Both methods must evaluate the same tumors.")

    required_metadata = {"tumor_id", "mutation_count"}

    if not required_metadata.issubset(metadata.columns):
        raise ValueError(
            "Metadata must contain tumor_id and mutation_count."
        )

    if metadata["tumor_id"].duplicated().any():
        raise ValueError("Metadata tumor IDs must be unique.")

    if set(metadata["tumor_id"]) != set(nnls_metrics.index):
        raise ValueError(
            "Metadata and evaluation results must contain the same tumors."
        )

    for metric in metrics:
        if metric not in nnls_metrics.columns:
            raise ValueError(f"Missing NNLS metric: {metric}")

        if metric not in regularized_metrics.columns:
            raise ValueError(
                f"Missing regularized NNLS metric: {metric}"
            )

    # Align data using tumor IDs.
    tumor_ids = nnls_metrics.index

    metadata_by_tumor = metadata.set_index("tumor_id").reindex(tumor_ids)

    nnls_aligned = nnls_metrics.reindex(tumor_ids)
    regularized_aligned = regularized_metrics.reindex(tumor_ids)

    summary_rows = []
    test_rows = []

    mutation_burdens = sorted(
        metadata_by_tumor["mutation_count"].unique()
    )

    for mutation_count in mutation_burdens:

        # Identify tumors belonging to this mutation burden.
        mask = (
            metadata_by_tumor["mutation_count"] == mutation_count
        )

        nnls_group = nnls_aligned.loc[mask]
        regularized_group = regularized_aligned.loc[mask]

        for method_name, group in [
            ("NNLS", nnls_group),
            ("Regularized NNLS", regularized_group)
        ]:

            for metric in metrics:

                values = group[metric]

                if not values.map(
                    lambda value: pd.notna(value) and abs(value) != float("inf")
                ).all():
                    raise ValueError(
                        f"Invalid values found in {method_name}: {metric}"
                    )

                summary_rows.append({
                    "mutation_count": mutation_count,
                    "method": method_name,
                    "metric": metric,
                    "n": len(values),
                    "mean": values.mean(),
                    "std": values.std()
                })

        # Compare the two methods on the same tumors.
        for metric in metrics:

            result = paired_wilcoxon(
                nnls_group[metric],
                regularized_group[metric],
                alternative="two-sided"
            )

            test_rows.append({
                "mutation_count": mutation_count,
                "metric": metric,
                "test": result["test"],
                "statistic": result["statistic"],
                "p_value": result["p_value"],
                "median_difference_NNLS_minus_regularized":
                    result["median_difference"],
                "n": result["n"]
            })

    summary = pd.DataFrame(summary_rows)
    tests = pd.DataFrame(test_rows)

    # Apply Holm correction separately for each metric.
    tests["adjusted_p_value"] = (
        tests.groupby("metric")["p_value"]
        .transform(lambda values: holm_correction(values.to_numpy()))
    )

    tests["significant"] = tests["adjusted_p_value"] < 0.05

    return summary, tests
