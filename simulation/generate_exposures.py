
import numpy as np
import pandas as pd


def generate_exposures(
    signature_names,
    rng,
    min_active=1,
    max_active=3,
    alpha=1.0
):
    """
    Generate ground-truth signature exposure proportions.

    Parameters
    ----------
    signature_names : list[str]
        Available COSMIC signature names.

    rng : numpy.random.Generator
        Random number generator for reproducibility.

    min_active : int, default=1
        Minimum number of active signatures per tumor.

    max_active : int, default=3
        Maximum number of active signatures per tumor.

    alpha : float, default=1.0
        Dirichlet concentration parameter.

    Returns
    -------
    pandas.Series
        Exposure proportions indexed by signature name.
        Inactive signatures have zero exposure.
    """

    n_signatures = len(signature_names)

    if n_signatures == 0:
        raise ValueError("At least one signature is required.")

    if not 1 <= min_active <= max_active <= n_signatures:
        raise ValueError(
            "Require 1 <= min_active <= max_active "
            "<= number of available signatures."
        )

    if alpha <= 0:
        raise ValueError("alpha must be positive.")

    # Randomly determine the number of active signatures.
    n_active = rng.integers(
        min_active,
        max_active + 1
    )

    # Select distinct signatures.
    active_indices = rng.choice(
        n_signatures,
        size=n_active,
        replace=False
    )

    # Generate exposure proportions.
    active_weights = rng.dirichlet(
        np.full(n_active, alpha)
    )

    # Initialize all exposures to zero.
    exposures = np.zeros(n_signatures)

    exposures[active_indices] = active_weights

    return pd.Series(
        exposures,
        index=signature_names,
        dtype=float
    )
