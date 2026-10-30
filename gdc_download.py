"""We download open-access TCGA mutation MAFs and clinical data from the GDC.

    python gdc_download.py --dry-run          # list what would be downloaded
    python gdc_download.py                    # download LUAD, SKCM, COAD

Files land in data/raw/TCGA-<TYPE>/<file_id>/<file_name>; clinical data in
data/raw/clinical/cases.json. We check each MAF against the MD5 checksum
the GDC reports, and rewrite data/manifest.csv at the end (Q6).
No login is needed: we only request open-access files.
"""
import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

from sigdet.storage import ensure_dirs, write_manifest

API = "https://api.gdc.cancer.gov"
# The GDC workflow that merges several callers' results and masks germline/low-quality calls.
WORKFLOW = "Aliquot Ensemble Somatic Variant Merging and Masking"


def _post(endpoint, payload, retries=3):
    """We send one JSON query to the GDC API and return the parsed reply."""
    data = json.dumps(payload).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(f"{API}/{endpoint}", data=data,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except Exception as e:  # network hiccup: retry with backoff
            if attempt == retries - 1:
                raise
            print(f"  retrying {endpoint} after error: {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))


def maf_filters(projects):
    """Our GDC search filter: open-access masked somatic MAFs from the ensemble workflow."""
    return {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": projects}},
        {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}},
        {"op": "=", "content": {"field": "analysis.workflow_type", "value": WORKFLOW}},
        {"op": "=", "content": {"field": "access", "value": "open"}},
    ]}


def list_mafs(projects):
    """We list matching MAF files with the metadata we need to download and check them."""
    res = _post("files", {
        "filters": maf_filters(projects),
        "fields": "file_id,file_name,file_size,md5sum,cases.project.project_id,cases.submitter_id",
        "format": "JSON", "size": "20000"})  # size: more than enough to get every hit in one page
    return res["data"]["hits"]


def fetch_cases(projects):
    """We fetch each patient's smoking status and pack-years."""
    res = _post("cases", {
        "filters": {"op": "in", "content": {"field": "project.project_id", "value": projects}},
        "fields": "submitter_id,project.project_id,exposures.tobacco_smoking_status,exposures.pack_years_smoked",
        "format": "JSON", "size": "20000"})
    return res["data"]["hits"]


def md5_file(path):
    """We compute a file's MD5 checksum in 1 MB chunks (the GDC publishes MD5)."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        while block := f.read(1 << 20):
            h.update(block)
    return h.hexdigest()


def download(hit, raw_dir):
    """We download one MAF unless an intact copy is already on disk."""
    project = hit["cases"][0]["project"]["project_id"]
    out = raw_dir / project / hit["file_id"] / hit["file_name"]
    # Re-running is cheap: files that already match their checksum are kept.
    if out.exists() and md5_file(out) == hit["md5sum"]:
        return out, "kept"
    out.parent.mkdir(parents=True, exist_ok=True)
    # We download to a .part file first, so an interrupted download never looks complete.
    tmp = out.with_suffix(out.suffix + ".part")
    urllib.request.urlretrieve(f"{API}/data/{hit['file_id']}", tmp)
    if md5_file(tmp) != hit["md5sum"]:
        tmp.unlink()
        raise RuntimeError(f"MD5 mismatch for {hit['file_name']}")
    tmp.replace(out)
    return out, "downloaded"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancers", nargs="+", default=["LUAD", "SKCM", "COAD"])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    projects = [f"TCGA-{c}" for c in args.cancers]
    paths = ensure_dirs()

    # We list what is available first and print counts per project.
    hits = list_mafs(projects)
    by_project = {}
    for h in hits:
        by_project.setdefault(h["cases"][0]["project"]["project_id"], []).append(h)
    for p in projects:
        files = by_project.get(p, [])
        mb = sum(int(h.get("file_size", 0)) for h in files) / 1e6
        print(f"{p}: {len(files)} MAF files, {mb:.1f} MB")
    if not hits:
        sys.exit("No files matched. Check the GDC filters in maf_filters().")
    if args.dry_run:
        return

    for i, h in enumerate(hits, 1):
        _, status = download(h, paths["raw"])
        if i % 50 == 0 or i == len(hits):
            print(f"  {i}/{len(hits)} MAFs ({status})")

    # Clinical data comes from a separate endpoint and is small, so we save it as one JSON file.
    cases = fetch_cases(projects)
    clin = paths["raw"] / "clinical" / "cases.json"
    clin.parent.mkdir(parents=True, exist_ok=True)
    clin.write_text(json.dumps(cases, indent=1))
    print(f"Clinical records: {len(cases)} -> {clin}")

    # The manifest is what we commit instead of the raw files.
    write_manifest(paths["raw"], paths["manifest"])
    print(f"Manifest written: {paths['manifest']}")


if __name__ == "__main__":
    main()
