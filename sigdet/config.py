"""Our project-wide settings.

SEED and SPLIT_FRACTIONS define our train/validation/test split. We fixed
them once, before any results existed; changing them now would change every
number we report.
"""
from pathlib import Path

# Project root (the folder above sigdet/), so paths work from any directory.
ROOT = Path(__file__).resolve().parent.parent

# One seed for everything random: the split, the simulator and the bootstraps.
SEED = 6721
SPLIT_FRACTIONS = {"train": 0.6, "val": 0.2, "test": 0.2}

# Q3: we count a signature as "present" when it explains at least 5% of a
# tumor's mutations. Below that, our estimates are mostly noise.
PRESENCE_THRESHOLD = 0.05

# Breakdown point N*: the smallest mutation count at which the median cosine
# similarity between true and estimated exposures reaches this value.
BREAKDOWN_COSINE = 0.90

# Mutation counts we test in E2 (accuracy vs mutation count).
N_GRID = [25, 50, 100, 250, 500, 1000, 2500, 5000]

# The 96 SBS contexts: 6 substitution types x 4 bases before x 4 bases after,
# written in COSMIC style, for example "A[C>T]G" (a C>T change between A and G).
MUTATION_TYPES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = "ACGT"
CONTEXTS = [f"{a}[{m}]{b}" for m in MUTATION_TYPES for a in BASES for b in BASES]

# Signatures COSMIC labels as likely sequencing artifacts. We never fit these.
ARTIFACT_SIGNATURES = (
    {"SBS27", "SBS43", "SBS95"} | {f"SBS{i}" for i in range(45, 61)}
)

# Our candidate signatures per cancer type, from derive_signature_subsets.py:
# non-artifact signatures that Alexandrov et al. 2020 attributed in at least
# 5% of that type's TCGA exomes (details in results/signature_subsets.csv).
# We write the old SBS40 as its COSMIC v3.4 successor SBS40a.
CANCER_SIGNATURES = {
    "LUAD": ["SBS1", "SBS2", "SBS4", "SBS5", "SBS13", "SBS40a"],
    "SKCM": ["SBS1", "SBS5", "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"],
    "COAD": ["SBS1", "SBS2", "SBS5", "SBS13", "SBS15", "SBS17b", "SBS40a", "SBS44"],
}
