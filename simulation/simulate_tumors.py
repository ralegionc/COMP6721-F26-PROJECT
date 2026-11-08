
from pathlib import Path

import numpy as np
import pandas as pd
from simulation.generate_exposures import generate_exposures


def generate_mutation_profile(
    signature_matrix,
    exposures,
    mutation_count,
    rng
):
    """
    Generate one synthetic SBS96 mutation profile.

    Parameters
    ----------
    signature_matrix : pandas.DataFrame
        COSMIC SBS96 matrix (96 contexts x K signatures).

    exposures : pandas.Series
        Ground-truth signature exposure proportions.

    mutation_count : int
        Total number of mutations to simulate.

    rng : numpy.random.Generator
        Random number generator.

    Returns
    -------
    pandas.Series
        Simulated mutation counts across 96 contexts.
    """

    if signature_matrix.shape[0] != 96:
        raise ValueError("Expected exactly 96 SBS contexts.")

    if mutation_count <= 0:
        raise ValueError("Mutation count must be positive.")

    if not signature_matrix.index.is_unique:
        raise ValueError("Mutation context names must be unique.")

    if not signature_matrix.columns.is_unique:
        raise ValueError("Signature names must be unique.")

    # Align exposures with COSMIC signature columns.
    if set(exposures.index) != set(signature_matrix.columns):
        raise ValueError(
            "Exposure names must match COSMIC signature columns."
        )

    weights = exposures.reindex(
        signature_matrix.columns
    ).to_numpy(dtype=float)

    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("Exposure weights must be finite and non-negative.")

    if not np.isclose(weights.sum(), 1.0):
        raise ValueError("Exposure weights must sum to 1.")

    # Convert COSMIC matrix to NumPy.
    signatures = signature_matrix.to_numpy(dtype=float)

    if not np.isfinite(signatures).all() or np.any(signatures < 0):
        raise ValueError("Signature matrix must be finite and non-negative.")

    # Normalize each signature independently.
    column_sums = signatures.sum(axis=0)

    if np.any(column_sums <= 0):
        raise ValueError("Every signature must have a positive sum.")

    normalized_signatures = signatures / column_sums

    # Combine signatures using their exposure weights.
    probabilities = normalized_signatures @ weights

    # Normalize the combined probability distribution.
    probabilities /= probabilities.sum()

    # Sample mutation counts.
    counts = rng.multinomial(
        mutation_count,
        probabilities
    )

    return pd.Series(
        counts,
        index=signature_matrix.index,
        dtype=int
    )


def simulate_dataset(
    signature_matrix,
    n_tumors=100,
    mutation_counts=(100, 500, 1000, 5000, 10000),
    min_active=1,
    max_active=3,
    alpha=1.0,
    random_state=42,
    output_dir="data/processed/simulated"
):
    """
    Generate multiple synthetic tumors and save their
    SBS96 profiles and ground-truth exposures.

    Returns
    -------
    profiles : pandas.DataFrame
        Tumors x 96 mutation contexts.

    exposures_df : pandas.DataFrame
        Tumors x K COSMIC signatures.

    metadata : pandas.DataFrame
        Simulation settings for each tumor.
    """

    if n_tumors <= 0:
        raise ValueError("n_tumors must be positive.")

    if not mutation_counts or any(n <= 0 for n in mutation_counts):
        raise ValueError("Mutation counts must be positive.")

    if signature_matrix.empty:
        raise ValueError("Signature matrix cannot be empty.")

    rng = np.random.default_rng(random_state)

    signature_names = signature_matrix.columns.tolist()

    profiles = []
    exposures_list = []
    metadata = []

    # Generate tumors at each mutation burden.
    for mutation_count in mutation_counts:

        for tumor_number in range(n_tumors):

            tumor_id = (
                f"TUMOR_{mutation_count}_{tumor_number + 1:04d}"
            )

            # Step 1: Generate ground-truth exposures.
            exposure = generate_exposures(
                signature_names=signature_names,
                rng=rng,
                min_active=min_active,
                max_active=max_active,
                alpha=alpha
            )

            # Step 2: Generate SBS96 mutation counts.
            profile = generate_mutation_profile(
                signature_matrix=signature_matrix,
                exposures=exposure,
                mutation_count=mutation_count,
                rng=rng
            )

            profiles.append(profile.rename(tumor_id))
            exposures_list.append(exposure.rename(tumor_id))

            metadata.append({
                "tumor_id": tumor_id,
                "mutation_count": mutation_count,
                "n_active_signatures": int(
                    (exposure > 0).sum()
                )
            })

    # Convert results to DataFrames.
    profiles_df = pd.DataFrame(profiles)
    exposures_df = pd.DataFrame(exposures_list)
    metadata_df = pd.DataFrame(metadata)

    # Save results.
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    profiles_df.to_csv(
        output_dir / "simulated_sbs96.csv",
        index_label="tumor_id"
    )

    exposures_df.to_csv(
        output_dir / "simulated_exposures.csv",
        index_label="tumor_id"
    )

    metadata_df.to_csv(
        output_dir / "simulation_metadata.csv",
        index=False
    )

    print("Simulation completed successfully.")
    print(f"Total tumors: {len(profiles_df)}")
    print(f"Mutation contexts: {profiles_df.shape[1]}")
    print(f"COSMIC signatures: {exposures_df.shape[1]}")

    return profiles_df, exposures_df, metadata_df
