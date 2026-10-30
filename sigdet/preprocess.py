"""We turn MAF files into 96-context mutation catalogs (one row per tumor sample).

Our steps, each counted in a QC report so losses are visible:
1. We read the MAF (tab-separated, '#' header lines skipped).
2. We keep single-base substitutions (Variant_Type == SNP, A/C/G/T alleles).
3. We keep variants called by MuTect2 (caller column detected automatically:
   "callers" in GDC ensemble MAFs, "CENTERS" in the MC3 public MAF).
4. We keep tumor sample types 01 (primary) and 06 (metastatic) by default.
5. We remove duplicate calls (same sample, position and alternate base).
6. We look up the base before and after each mutation in the reference
   genome. If the MAF's reference base disagrees with the genome too often,
   the genome build is wrong (for example GRCh37 FASTA with GRCh38 MAFs) and
   we stop.
7. We fold onto the pyrimidine strand: a G>A on one strand is a C>T on the
   other, so we reverse-complement mutations of G or A. This gives the
   standard 96 contexts (6 substitution types x 16 flanking-base pairs).
"""
import gzip
from pathlib import Path

import pandas as pd

from .config import CONTEXTS

# MAF columns we cannot work without, and ones we read if present.
REQUIRED = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position",
            "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type"]
OPTIONAL = ["NCBI_Build"]
CALLER_COLUMNS = ["callers", "CALLERS", "CENTERS"]
# Base-pairing table for complementing DNA (N stays N).
_COMP = str.maketrans("ACGTN", "TGCAN")


class GenomeBuildError(RuntimeError):
    """MAF reference alleles do not match the reference genome."""


def revcomp(seq):
    """We return the reverse complement: the same DNA read from the other strand."""
    return seq.translate(_COMP)[::-1]


def _open(path):
    """We open a MAF as text, whether or not it is gzip-compressed."""
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def read_maf(path):
    """We read the columns our pipeline needs from one MAF file."""
    # First pass: count the '#' comment lines and read the header row.
    with _open(path) as f:
        skip = 0
        for line in f:
            if line.startswith("#"):
                skip += 1
            else:
                header = line.rstrip("\n").split("\t")
                break
        else:
            raise ValueError(f"{path}: no header line")
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    # Second pass: load only the columns we use, all as text, to save memory.
    wanted = set(REQUIRED + OPTIONAL + CALLER_COLUMNS)
    with _open(path) as f:
        df = pd.read_csv(f, sep="\t", skiprows=skip, dtype=str, usecols=lambda c: c in wanted,
                         low_memory=False, keep_default_na=False)
    df["source_file"] = Path(path).name
    return df


def _caller_column(df):
    """We return the name of whichever caller column this MAF has, or None."""
    for c in CALLER_COLUMNS:
        if c in df.columns:
            return c
    return None


def filter_mutations(df, caller="mutect", sample_types=("01", "06"), qc=None):
    """We apply steps 2-5. qc (dict) receives the row count after each step."""
    qc = {} if qc is None else qc
    qc["rows_read"] = len(df)

    # Step 2: one base changed into a different single base.
    alleles = set("ACGT")
    snv = (df["Variant_Type"] == "SNP") & df["Reference_Allele"].isin(alleles) \
        & df["Tumor_Seq_Allele2"].isin(alleles) & (df["Reference_Allele"] != df["Tumor_Seq_Allele2"])
    df = df[snv]
    qc["after_snv_filter"] = len(df)

    # Step 3: we stop rather than silently skip the caller filter.
    if caller:
        col = _caller_column(df)
        if col is None:
            raise ValueError(
                f"No caller column ({CALLER_COLUMNS}) found, so MuTect2 calls cannot be selected. "
                "Inspect the MAF header, or rerun with --no-caller-filter and report it as a deviation.")
        df = df[df[col].str.lower().str.contains(caller.lower(), regex=False)]
        qc["caller_column"] = col
    qc["after_caller_filter"] = len(df)

    # Step 4: characters 14-15 of the barcode give the sample type (01 primary, 06 metastatic).
    if sample_types:
        st = df["Tumor_Sample_Barcode"].str.slice(13, 15)
        df = df[st.isin(sample_types)]
    qc["after_sample_type_filter"] = len(df)

    # Step 5: the same mutation listed twice for one sample counts once.
    df = df.drop_duplicates(["Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Tumor_Seq_Allele2"])
    qc["after_dedup"] = len(df)
    return df


def assign_contexts(df, reference, max_ref_mismatch=0.01, qc=None):
    """We apply steps 6-7 and return df with a 'context' column, unusable rows removed."""
    qc = {} if qc is None else qc
    labels, mismatches, unresolved = [], 0, 0
    for chrom, pos, ref, alt in zip(df["Chromosome"], df["Start_Position"].astype(int),
                                    df["Reference_Allele"], df["Tumor_Seq_Allele2"]):
        # The three bases centred on the mutation.
        tri = reference.fetch(chrom, pos - 1, pos + 1)
        if len(tri) != 3 or "N" in tri:
            unresolved += 1
            labels.append(None)
            continue
        # The middle base must equal the MAF's reference base; if not, we count a mismatch.
        if tri[1] != ref:
            mismatches += 1
            labels.append(None)
            continue
        # Step 7: we express G and A mutations from the opposite strand.
        if ref in "GA":
            tri, ref, alt = revcomp(tri), ref.translate(_COMP), alt.translate(_COMP)
        labels.append(f"{tri[0]}[{ref}>{alt}]{tri[2]}")
    checked = len(df) - unresolved
    rate = mismatches / checked if checked else 0.0
    qc.update({"unresolved_positions": unresolved, "ref_mismatches": mismatches,
               "ref_mismatch_rate": round(rate, 5)})
    # A few mismatches are normal; many mean we are using the wrong genome build.
    if checked and rate > max_ref_mismatch:
        builds = sorted(set(df["NCBI_Build"])) if "NCBI_Build" in df.columns else ["unknown"]
        raise GenomeBuildError(
            f"{rate:.1%} of reference alleles disagree with the genome FASTA (MAF build: {builds}). "
            "Use the matching reference: GRCh38 for GDC MAFs, GRCh37 for MC3.")
    out = df.assign(context=labels)
    out = out[out["context"].notna()]
    qc["mutations_in_catalog"] = len(out)
    return out


def build_catalog(mutations):
    """We count mutations per sample and context (samples x 96, standard COSMIC order)."""
    tab = pd.crosstab(mutations["Tumor_Sample_Barcode"], mutations["context"])
    # Contexts that never occur still get a column of zeros.
    return tab.reindex(columns=CONTEXTS, fill_value=0).astype(int).rename_axis(index="sample", columns=None)


def preprocess_mafs(maf_paths, reference, caller="mutect", sample_types=("01", "06"),
                    max_ref_mismatch=0.01):
    """We run the full pipeline for one cancer type and return (catalog, qc dict)."""
    maf_paths = [Path(p) for p in maf_paths]
    if not maf_paths:
        raise FileNotFoundError("No MAF files given")
    qc = {"maf_files": len(maf_paths)}
    df = pd.concat([read_maf(p) for p in maf_paths], ignore_index=True)
    df = filter_mutations(df, caller=caller, sample_types=sample_types, qc=qc)
    df = assign_contexts(df, reference, max_ref_mismatch=max_ref_mismatch, qc=qc)
    catalog = build_catalog(df)
    qc["samples"] = int(len(catalog))
    if len(catalog):
        # Mutation burden per sample, which tells us how many tumors sit below N*.
        burden = catalog.sum(axis=1)
        qc["mutations_per_sample_median"] = float(burden.median())
        qc["samples_under_100_mutations"] = int((burden < 100).sum())
    return catalog, qc
