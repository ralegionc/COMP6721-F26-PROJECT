
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from classical_methods.nnls import prepare_signature_matrix

def estimate_exposures_regularized_nnls(
    profiles,
    signature_matrix,
    alpha=0.01
):
    """
    Estimate signature exposures using ridge-regularized NNLS.

    Parameters
    ----------
    profiles : pandas.DataFrame
        Tumors x 96 mutation contexts.

    signature_matrix : pandas.DataFrame
        96 mutation contexts x K COSMIC signatures.

    alpha : float, default=0.01
        Regularization strength. Must be non-negative.

    Returns
    -------
    pandas.DataFrame
        Estimated exposure proportions.
    """

    if not isinstance(profiles, pd.DataFrame):
        raise TypeError("profiles must be a pandas DataFrame.")

    if not profiles.index.is_unique:
        raise ValueError("Tumor IDs must be unique.")

    if not profiles.columns.is_unique:
        raise ValueError("Mutation context names must be unique.")

    if not profiles.columns.equals(signature_matrix.index):
        raise ValueError(
            "Profile contexts and signature matrix contexts "
            "must match in the same order."
        )

    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and non-negative.")

    A = prepare_signature_matrix(signature_matrix)

    Y = profiles.to_numpy(dtype=float)

    if not np.isfinite(Y).all() or np.any(Y < 0):
        raise ValueError("Mutation profiles must be finite and non-negative.")

    totals = Y.sum(axis=1)

    if np.any(totals <= 0):
        raise ValueError("Every tumor must contain at least one mutation.")

    Y_normalized = Y / totals[:, np.newaxis]

    # Augment the least-squares problem to include the L2 penalty.
    identity = np.eye(A.shape[1])

    A_augmented = np.vstack([
        A,
        np.sqrt(alpha) * identity
    ])

    estimates = []

    for tumor_profile in Y_normalized:

        y_augmented = np.concatenate([
            tumor_profile,
            np.zeros(A.shape[1])
        ])

        # Solve the regularized non-negative least squares problem.
        coefficients, _ = nnls(
            A_augmented,
            y_augmented
        )

        coefficient_sum = coefficients.sum()

        if coefficient_sum <= 0:
            raise ValueError(
                "Regularized NNLS returned an empty exposure estimate."
            )

        exposures = coefficients / coefficient_sum

        estimates.append(exposures)

    return pd.DataFrame(
        estimates,
        index=profiles.index,
        columns=signature_matrix.columns
    )
