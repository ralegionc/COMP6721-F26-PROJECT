"""Build data/processed/reference_signatures.csv from the COSMIC v3.4 file,
reordered to match the context order of the TCGA SBS96 catalogs."""
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent  # the project folder, wherever you run from


def find_file(name, folders):
    """Return the first existing BASE/<folder>/<name>, or raise a clear error."""
    candidates = [BASE / folder / name for folder in folders]
    for path in candidates:
        if path.exists():
            return path
    searched = "\n  ".join(str(p) for p in candidates)
    raise FileNotFoundError(f"Could not find {name}. Looked in:\n  {searched}")


COSMIC = find_file(
    "COSMIC_v3.4_SBS_GRCh38_exome.txt",
    ["data/reference", "data/raw", "data", "."],
)
CATALOG = find_file(
    "COAD_sbs96.csv",  # all three catalogs share the same column order
    ["data/processed", "data", "."],
)

cosmic = pd.read_csv(COSMIC, sep="\t", index_col=0)
contexts = pd.read_csv(CATALOG, index_col=0, nrows=1).columns

if set(cosmic.index) != set(contexts):
    raise ValueError("COSMIC and catalog contexts differ.")

cosmic = cosmic.reindex(contexts)  # same 96 contexts, catalog order
cosmic.index.name = "Type"

output = BASE / "data" / "processed" / "reference_signatures.csv"
output.parent.mkdir(parents=True, exist_ok=True)
cosmic.to_csv(output)

print(cosmic.shape, "->", output)