"""We turn raw MAFs into 96-context catalogs plus a clinical table (data/processed/).

    python preprocess_data.py --reference GRCh38.d1.vd1.fa

Outputs per cancer type:
    data/processed/<TYPE>_sbs96.csv   samples x 96 contexts (counts)
    data/processed/qc_<TYPE>.json     rows kept after each filter step
and data/processed/clinical.csv (smoking status, pack-years, smoker label).
"""
import argparse
import json

from sigdet.clinical import load_clinical, unmapped_smoking_values
from sigdet.preprocess import preprocess_mafs
from sigdet.reference import FastaReference
from sigdet.storage import ensure_dirs, verify_manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, help="uncompressed FASTA matching the MAF build")
    ap.add_argument("--cancers", nargs="+", default=["LUAD", "SKCM", "COAD"])
    ap.add_argument("--no-caller-filter", action="store_true",
                    help="keep all callers (only if the MAF has no caller column; report as deviation)")
    args = ap.parse_args()
    paths = ensure_dirs()

    # Before using the raw files, we confirm they are exactly the ones we downloaded.
    if paths["manifest"].exists():
        problems = verify_manifest(paths["manifest"], paths["raw"])
        if problems:
            raise SystemExit(f"Raw data differs from manifest: {problems[:5]}")
        print("Raw files match data/manifest.csv")

    with FastaReference(args.reference) as ref:
        for cancer in args.cancers:
            mafs = sorted((paths["raw"] / f"TCGA-{cancer}").rglob("*.maf*"))
            mafs = [m for m in mafs if not m.name.endswith(".part")]  # skip unfinished downloads
            print(f"{cancer}: {len(mafs)} MAF files")
            catalog, qc = preprocess_mafs(mafs, ref, caller=None if args.no_caller_filter else "mutect")
            catalog.to_csv(paths["processed"] / f"{cancer}_sbs96.csv")
            (paths["processed"] / f"qc_{cancer}.json").write_text(json.dumps(qc, indent=2))
            print(f"  {qc['samples']} samples, median {qc.get('mutations_per_sample_median', 0):.0f} "
                  f"mutations, ref mismatch {qc['ref_mismatch_rate']:.2%}")

    clin_json = paths["raw"] / "clinical" / "cases.json"
    if clin_json.exists():
        clin = load_clinical(clin_json)
        clin.to_csv(paths["processed"] / "clinical.csv")
        print(f"Clinical: {len(clin)} patients, smoker={int((clin['smoker'] == 1).sum())}, "
              f"never={int((clin['smoker'] == 0).sum())}")
        # Any smoking value we could not map is printed, so we never drop labels silently.
        odd = unmapped_smoking_values(clin)
        if len(odd):
            print("Unrecognised smoking values (check labels.py):\n" + odd.to_string())


if __name__ == "__main__":
    main()
