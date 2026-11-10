"""We rerun the paired method comparisons without redoing the whole experiment.

    python run_method_tests.py            # real data (results/per_tumor_cosine_real.csv)
    python run_method_tests.py --demo     # demo data

run_checklist.py already does this at the end of every run; this script is
for re-testing saved per-tumor scores. For every cancer type and every pair
of methods we run a paired Wilcoxon signed-rank test on per-tumor cosine and
MAE (active signatures), pooled over mutation counts and at each count, and
apply Holm's correction. We write results/method_tests_<source>.csv.
"""
import argparse

import pandas as pd

from sigdet.metrics import compare_methods
from sigdet.storage import ensure_dirs, save_results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true", help="use the demo run's per-tumor scores")
    args = ap.parse_args()
    source = "demo" if args.demo else "real"
    paths = ensure_dirs()  # absolute paths, so the script works from any folder

    per_tumor = pd.read_csv(paths["results"] / f"per_tumor_cosine_{source}.csv")
    if "tumor" not in per_tumor.columns:
        raise SystemExit("This per-tumor file predates the 'tumor' column; rerun run_checklist.py first.")
    tests = compare_methods(per_tumor)
    save_results(tests, f"method_tests_{source}", paths)

    pooled = tests[tests["n"] == "all"]
    cols = ["cancer", "metric", "method_a", "method_b", "median_difference", "a_higher_fraction", "p_holm"]
    print(pooled[cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"\n{int(tests['significant'].sum())} of {len(tests)} comparisons significant after Holm "
          f"correction; wrote results/method_tests_{source}.csv")


if __name__ == "__main__":
    main()
