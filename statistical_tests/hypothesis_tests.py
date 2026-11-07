
import numpy as np
import pandas as pd

from scipy.stats import (
    wilcoxon,
    mannwhitneyu,
    kruskal
)


def _validate_samples(*samples):
    """Validate numerical samples used by statistical tests."""

    for sample in samples:
        values = np.asarray(sample, dtype=float)

        if values.ndim != 1:
            raise ValueError("Samples must be one-dimensional.")

        if len(values) == 0:
            raise ValueError("Samples cannot be empty.")

        if not np.isfinite(values).all():
            raise ValueError("Samples must contain finite values.")


def paired_wilcoxon(
    sample_a,
    sample_b,
    alternative="two-sided"
):
    """
    Perform a Wilcoxon signed-rank test on paired observations.

    A positive median difference means sample_a tends to be larger.
    """

    _validate_samples(sample_a, sample_b)

    a = np.asarray(sample_a, dtype=float)
    b = np.asarray(sample_b, dtype=float)

    if len(a) != len(b):
        raise ValueError("Paired samples must have equal lengths.")

    if alternative not in {"two-sided", "less", "greater"}:
        raise ValueError("Invalid alternative hypothesis.")

    differences = a - b

    # Handle identical samples explicitly.
    if np.all(differences == 0):
        statistic = 0.0
        p_value = 1.0
    else:
        statistic, p_value = wilcoxon(
            a,
            b,
            alternative=alternative
        )

    return {
        "test": "Wilcoxon signed-rank",
        "statistic": float(statistic),
        "p_value": float(p_value),
        "median_difference": float(np.median(differences)),
        "n": len(a)
    }


def independent_mannwhitney(
    sample_a,
    sample_b,
    alternative="two-sided"
):
    """
    Perform a Mann-Whitney U test on independent samples.

    Effect size is rank-biserial correlation.
    """

    _validate_samples(sample_a, sample_b)

    a = np.asarray(sample_a, dtype=float)
    b = np.asarray(sample_b, dtype=float)

    if alternative not in {"two-sided", "less", "greater"}:
        raise ValueError("Invalid alternative hypothesis.")

    statistic, p_value = mannwhitneyu(
        a,
        b,
        alternative=alternative
    )

    # Rank-biserial effect size.
    effect_size = (
        2 * statistic / (len(a) * len(b))
    ) - 1

    return {
        "test": "Mann-Whitney U",
        "statistic": float(statistic),
        "p_value": float(p_value),
        "rank_biserial": float(effect_size),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "n_a": len(a),
        "n_b": len(b)
    }


def independent_kruskal(*groups):
    """
    Perform a Kruskal-Wallis test on two or more independent groups.

    Reports epsilon-squared as an effect size.
    """

    if len(groups) < 2:
        raise ValueError("At least two groups are required.")

    _validate_samples(*groups)

    arrays = [
        np.asarray(group, dtype=float)
        for group in groups
    ]

    statistic, p_value = kruskal(*arrays)

    total_n = sum(len(group) for group in arrays)
    number_of_groups = len(arrays)

    # Kruskal-Wallis epsilon-squared effect size: H / (n - 1).
    if total_n > 1:
        effect_size = statistic / (total_n - 1)
    else:
        effect_size = np.nan

    return {
        "test": "Kruskal-Wallis",
        "statistic": float(statistic),
        "p_value": float(p_value),
        "epsilon_squared": float(effect_size),
        "n_groups": number_of_groups,
        "total_n": total_n
    }


def holm_correction(p_values):
    """
    Apply Holm-Bonferroni correction to multiple p-values.

    Returns adjusted p-values in their original order.
    """

    p_values = np.asarray(p_values, dtype=float)

    if p_values.ndim != 1:
        raise ValueError("p_values must be one-dimensional.")

    if not np.isfinite(p_values).all():
        raise ValueError("P-values must be finite.")

    if np.any((p_values < 0) | (p_values > 1)):
        raise ValueError("P-values must be between 0 and 1.")

    n = len(p_values)

    if n == 0:
        return np.array([], dtype=float)

    order = np.argsort(p_values)
    sorted_p = p_values[order]

    # Holm adjusted p-values.
    adjusted_sorted = np.maximum.accumulate(
        (n - np.arange(n)) * sorted_p
    )

    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)

    # Restore the original ordering.
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = adjusted_sorted

    return adjusted
