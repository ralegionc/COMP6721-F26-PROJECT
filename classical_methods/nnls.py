
import numpy as np
import pandas as pd
from scipy.optimize import nnls


def prepare_signature_matrix(signature_matrix):
    """
    Validate and normalize the COSMIC SBS96 signature matrix.

    Returns
    -------
    numpy.ndarray
        Normalized matrix with shape (96, K).
    """

    if not isinstance(signature_matrix, pd.DataFrame):
        raise TypeError("signature_matrix must be a pandas DataFrame.")

    if signature_matrix.shape[0] != 96:
        raise ValueError("Expected exactly 96 SBS contexts.")

    if signature_matrix.empty:
        raise ValueError("Signature matrix cannot be empty.")

    if not signature_matrix.index.is_unique:
        raise ValueError("Mutation context names must be unique.")

    if not signature_matrix.columns.is_unique:
        raise ValueError("Signature names must be unique.")

    matrix = signature_matrix.to_numpy(dtype=float)

    if not np.isfinite(matrix).all() or np.any(matrix < 0):
        raise ValueError("Signature matrix must be finite and non-negative.")

    column_sums = matrix.sum(axis=0)

    if np.any(column_sums <= 0):
        raise ValueError("Every signature must have a positive sum.")

    return matrix / column_sums


def estimate_exposures_nnls(
    profiles,
    signature_matrix
):
    """
    Estimate signature exposure proportions using NNLS.

    Parameters
    ----------
    profiles : pandas.DataFrame
        Tumors x 96 mutation contexts.

    signature_matrix : pandas.DataFrame
        96 mutation contexts x K COSMIC signatures.

    Returns
    -------
    pandas.DataFrame
        Estimated exposure proportions.
        Rows represent tumors and columns represent signatures.
    """

    if not isinstance(profiles, pd.DataFrame):
        raise TypeError("profiles must be a pandas DataFrame.")

    if not profiles.index.is_unique:
        raise ValueError("Tumor IDs must be unique.")

    if not profiles.columns.is_unique:
        raise ValueError("Mutation context names must be unique.")

    if not profiles.columns.equals(signature_matrix.index):
        raise ValueError(
            "Profile context names and signature matrix contexts "
            "must match in the same order."
        )

    A = prepare_signature_matrix(signature_matrix)

    Y = profiles.to_numpy(dtype=float)

    if not np.isfinite(Y).all() or np.any(Y < 0):
        raise ValueError("Mutation profiles must be finite and non-negative.")

    totals = Y.sum(axis=1)

    if np.any(totals <= 0):
        raise ValueError("Every tumor must contain at least one mutation.")

    # Convert mutation counts into probability distributions.
    Y_normalized = Y / totals[:, np.newaxis]

    estimates = []

    for tumor_profile in Y_normalized:

        # Solve the non-negative least squares problem.
        coefficients, _ = nnls(A, tumor_profile)

        coefficient_sum = coefficients.sum()

        if coefficient_sum <= 0:
            raise ValueError("NNLS returned an empty exposure estimate.")

        # Convert coefficients into exposure proportions.
        exposures = coefficients / coefficient_sum

        estimates.append(exposures)

    return pd.DataFrame(
        estimates,
        index=profiles.index,
        columns=signature_matrix.columns
    )
