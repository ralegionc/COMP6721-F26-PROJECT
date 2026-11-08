import pandas as pd

from simulation.simulate_tumors import simulate_dataset


# ---------------------------------------------------
# 1. Load the COSMIC reference signatures
# ---------------------------------------------------

signature_matrix = pd.read_csv(
    "data/processed/reference_signatures.csv",
    index_col=0
)


# ---------------------------------------------------
# 2. Simulate tumors and save them
# ---------------------------------------------------

output_dir = "data/processed/simulated"

profiles, exposures, metadata = simulate_dataset(
    signature_matrix,
    n_tumors=100,
    mutation_counts=(100, 500, 1000, 5000, 10000),
    min_active=1,
    max_active=3,
    alpha=1.0,
    random_state=42,
    output_dir=output_dir
)


# ---------------------------------------------------
# 3. Quick sanity check
# ---------------------------------------------------

print("\nTumors per mutation burden:")
print(metadata["mutation_count"].value_counts().sort_index())

print("\nActive signatures per tumor:")
print(metadata["n_active_signatures"].value_counts().sort_index())

print(f"\nSaved to {output_dir}/")