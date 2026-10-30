"""We derive the candidate signature list for each cancer type (sigdet/config.py).

    python derive_signature_subsets.py

Our source: the per-sample SigProfiler attributions for TCGA exomes published
with Alexandrov et al. 2020, "The repertoire of mutational signatures in human
cancer" (Nature 578:94-101), file TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv.
Our copy in data/reference/ comes from the PCAWG7 R package
(github.com/steverozen/PCAWG7, branch v0.1.4-branch, data-raw/sig.profiler.exposures/).

Our rule: we keep a signature for a cancer type when the published fit
attributed at least one mutation to it in at least MIN_SAMPLE_FRACTION of
that type's TCGA exomes, and COSMIC does not label it a sequencing artifact.
We leave out rarer signatures (for example the POLE signatures SBS10a/b in
about 2% of colorectal tumors); they serve as "unseen" signatures in E3.

The 2020 attributions use COSMIC v3.0 names. COSMIC v3.4 split SBS40 into
SBS40a, SBS40b and SBS40c (defined from kidney tumors); SBS40a is the closest
single successor (COSMIC reports cosine 0.91 with the old SBS40) and the one
COSMIC describes as attributable across many cancer types, so we map SBS40 to it.

We write results/signature_subsets.csv and print a CANCER_SIGNATURES block
to paste into sigdet/config.py.
"""
import argparse

import pandas as pd

from sigdet.config import ARTIFACT_SIGNATURES
from sigdet.storage import ensure_dirs

MIN_SAMPLE_FRACTION = 0.05
# Our cancer codes and the names the attribution file uses for them.
TYPE_NAMES = {"LUAD": "Lung-AdenoCa", "SKCM": "Skin-Melanoma", "COAD": "ColoRect-AdenoCa"}
# Old COSMIC names that changed by v3.4.
V34_NAMES = {"SBS40": "SBS40a"}


def signature_activity(attributions, type_name):
    """Per signature, we compute the fraction of samples with any attributed mutation, and the mean share."""
    x = attributions[attributions["Cancer Types"] == type_name]
    if x.empty:
        raise ValueError(f"No samples labelled {type_name!r} in the attribution file")
    sigs = [c for c in x.columns if c.startswith("SBS")]
    counts = x[sigs]  # mutations attributed to each signature in each sample
    return pd.DataFrame({
        "n_samples": len(x),
        "fraction_of_samples": (counts > 0).mean(),
        "mean_share": counts.div(counts.sum(axis=1), axis=0).mean(),
        "artifact": [s in ARTIFACT_SIGNATURES for s in sigs],
    }, index=pd.Index(sigs, name="signature"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--attributions",
                    default="data/reference/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
    ap.add_argument("--min-fraction", type=float, default=MIN_SAMPLE_FRACTION)
    args = ap.parse_args()
    paths = ensure_dirs()

    attributions = pd.read_csv(args.attributions)
    tables, chosen = [], {}
    for cancer, type_name in TYPE_NAMES.items():
        t = signature_activity(attributions, type_name)
        # Our selection rule: common enough, and not a sequencing artifact.
        t["selected"] = (t["fraction_of_samples"] >= args.min_fraction) & ~t["artifact"]
        t["v34_name"] = [V34_NAMES.get(s, s) for s in t.index]
        chosen[cancer] = [t.at[s, "v34_name"] for s in t.index if t.at[s, "selected"]]  # COSMIC order
        # For the record we keep every signature seen at all, selected or not.
        tables.append(t[t["fraction_of_samples"] > 0].assign(cancer=cancer).reset_index())

    out = pd.concat(tables)[["cancer", "signature", "v34_name", "n_samples", "fraction_of_samples",
                             "mean_share", "artifact", "selected"]]
    out.to_csv(paths["results"] / "signature_subsets.csv", index=False, float_format="%.4f")

    print("CANCER_SIGNATURES = {")
    for cancer, sigs in chosen.items():
        print(f'    "{cancer}": {sigs},')
    print("}")
    print(f"Details: {paths['results'] / 'signature_subsets.csv'}")


if __name__ == "__main__":
    main()
