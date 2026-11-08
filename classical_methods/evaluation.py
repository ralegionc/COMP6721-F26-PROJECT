
import numpy as np
import pandas as pd

from classical_methods.nnls import prepare_signature_matrix

def evaluate_exposures(
    true_exposures,
    estimated_exposures,
    support_threshold=1e-6,
    profiles=None,
    signature_matrix=None
):
    """
    Evaluate estimated exposures against ground-truth exposures.

    Parameters
    ----------
    true_exposures : pandas.DataFrame
        Ground-truth tumor x signature exposure proportions.

    estimated_exposures : pandas.DataFrame
        Estimated tumor x signature exposure proportions.

    support_threshold : float, default=1e-6
        Minimum exposure proportion considered active.

    profiles : pandas.DataFrame, optional
        Observed tumor mutation profiles. Required for reconstruction error.

    signature_matrix : pandas.DataFrame, optional
        COSMIC SBS96 signature matrix. Required for reconstruction error.

    Returns
    -------
    summary : dict
        Overall average evaluation metrics.

    per_tumor : pandas.DataFrame
        Evaluation metrics for individual tumors.
    """

    if not isinstance(true_exposures, pd.DataFrame):
        raise TypeError("true_exposures must be a pandas DataFrame.")

    if not isinstance(estimated_exposures, pd.DataFrame):
        raise TypeError("estimated_exposures must be a pandas DataFrame.")

    if not np.isfinite(support_threshold) or support_threshold < 0:
        raise ValueError("support_threshold must be finite and non-negative.")

    if not true_exposures.index.is_unique:
        raise ValueError("Ground-truth tumor IDs must be unique.")

    if not estimated_exposures.index.is_unique:
        raise ValueError("Estimated tumor IDs must be unique.")

    if not true_exposures.columns.is_unique:
        raise ValueError("Ground-truth signature names must be unique.")

    if not estimated_exposures.columns.is_unique:
        raise ValueError("Estimated signature names must be unique.")

    if set(true_exposures.index) != set(estimated_exposures.index):
        raise ValueError("Ground-truth and estimated tumor IDs must match.")

    if set(true_exposures.columns) != set(estimated_exposures.columns):
        raise ValueError(
            "Ground-truth and estimated signature names must match."
        )

    # Align both DataFrames before comparing them.
    estimated_exposures = estimated_exposures.reindex(
        index=true_exposures.index,
        columns=true_exposures.columns
    )

    true = true_exposures.to_numpy(dtype=float)
    estimated = estimated_exposures.to_numpy(dtype=float)

    if (
        not np.isfinite(true).all()
        or not np.isfinite(estimated).all()
        or np.any(true < 0)
        or np.any(estimated < 0)
    ):
        raise ValueError("Exposure values must be finite and non-negative.")

    if not np.allclose(true.sum(axis=1), 1.0):
        raise ValueError("Ground-truth exposures must sum to 1.")

    if not np.allclose(estimated.sum(axis=1), 1.0):
        raise ValueError("Estimated exposures must sum to 1.")

    # Calculate exposure estimation errors.
    absolute_errors = np.abs(true - estimated)
    squared_errors = (true - estimated) ** 2

    mae = absolute_errors.mean(axis=1)

    rmse = np.sqrt(
        squared_errors.mean(axis=1)
    )

    # Calculate cosine similarity for each tumor.
    dot_products = np.sum(true * estimated, axis=1)

    norms = (
        np.linalg.norm(true, axis=1)
        * np.linalg.norm(estimated, axis=1)
    )

    cosine_similarity = np.divide(
        dot_products,
        norms,
        out=np.zeros_like(dot_products),
        where=norms > 0
    )

    # Identify active signatures using the threshold.
    true_active = true > support_threshold
    estimated_active = estimated > support_threshold

    true_positive = np.sum(
        true_active & estimated_active,
        axis=1
    )

    false_positive = np.sum(
        ~true_active & estimated_active,
        axis=1
    )

    false_negative = np.sum(
        true_active & ~estimated_active,
        axis=1
    )

    precision = np.divide(
        true_positive,
        true_positive + false_positive,
        out=np.zeros_like(true_positive, dtype=float),
        where=(true_positive + false_positive) > 0
    )

    recall = np.divide(
        true_positive,
        true_positive + false_negative,
        out=np.zeros_like(true_positive, dtype=float),
        where=(true_positive + false_negative) > 0
    )

    f1 = np.divide(
        2 * precision * recall,
        precision + recall,
        out=np.zeros_like(precision),
        where=(precision + recall) > 0
    )

    per_tumor = pd.DataFrame({
        "MAE": mae,
        "RMSE": rmse,
        "Cosine Similarity": cosine_similarity,
        "Support Precision": precision,
        "Support Recall": recall,
        "Support F1": f1
    }, index=true_exposures.index)

    # Optionally evaluate mutation-profile reconstruction.
    if (profiles is None) != (signature_matrix is None):
        raise ValueError(
            "Provide both profiles and signature_matrix "
            "to calculate reconstruction error."
        )

    if profiles is not None:

        if not profiles.index.equals(true_exposures.index):
            raise ValueError(
                "Profile tumor IDs must match exposure tumor IDs "
                "in the same order."
            )

        if not profiles.columns.equals(signature_matrix.index):
            raise ValueError(
                "Profile contexts and signature matrix contexts "
                "must match in the same order."
            )

        A = prepare_signature_matrix(signature_matrix)

        Y = profiles.to_numpy(dtype=float)

        if not np.isfinite(Y).all() or np.any(Y < 0):
            raise ValueError(
                "Mutation profiles must be finite and non-negative."
            )

        totals = Y.sum(axis=1)

        if np.any(totals <= 0):
            raise ValueError("Every tumor must contain mutations.")

        Y_normalized = Y / totals[:, np.newaxis]

        # Reconstruct each tumor's SBS96 probability profile.
        reconstructed = estimated @ A.T

        errors = np.linalg.norm(
            Y_normalized - reconstructed,
            axis=1
        )

        denominators = np.linalg.norm(
            Y_normalized,
            axis=1
        )

        reconstruction_error = np.divide(
            errors,
            denominators,
            out=np.zeros_like(errors),
            where=denominators > 0
        )

        per_tumor["Reconstruction Error"] = reconstruction_error

    # Calculate average metrics across all tumors.
    summary = per_tumor.mean().to_dict()

    return summary, per_tumor
