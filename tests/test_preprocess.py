"""We test preprocessing on a tiny fake genome and MAF."""
import json

import numpy as np
import pandas as pd
import pytest

from sigdet.clinical import flatten_gdc_cases, load_clinical
from sigdet.methods import M0NNLS, M1ElasticNetNNLS, M2MultinomialEM
from sigdet.preprocess import GenomeBuildError, preprocess_mafs, revcomp
from sigdet.reference import FastaReference
from sigdet.target import demo_signatures, signature_subset
from sigdet.simulate import simulate

#            1234567890 1234567890 1234567890
CHR1 = "ACGTTGCAAC" "GTACGATCGA" "TTGCCAGTAA"


@pytest.fixture
def genome(tmp_path):
    fa = tmp_path / "ref.fa"
    lines = [CHR1[i:i + 7] for i in range(0, len(CHR1), 7)]  # 7 bases per line
    fa.write_text(">chr1 test\n" + "\n".join(lines) + "\n>chr2\nAAAACCCCGGGG\n")
    return fa


def write_maf(path, rows, caller_col="callers"):
    header = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Reference_Allele",
              "Tumor_Seq_Allele2", "Variant_Type", "NCBI_Build", caller_col]
    body = "\n".join("\t".join(map(str, r)) for r in rows)
    path.write_text("#version gdc-1.0.0\n" + "\t".join(header) + "\n" + body + "\n")
    return path


T1 = "TCGA-05-0001-01A-01D-A000-08"
NORMAL = "TCGA-05-0001-11A-01D-A000-08"


def test_reference_fetch_across_lines_and_names(genome):
    with FastaReference(genome) as ref:
        assert ref.fetch("chr1", 6, 9) == CHR1[5:9]          # crosses a line break
        assert ref.fetch("1", 1, 3) == "ACG"                 # "1" resolves to "chr1"
        assert ref.fetch("chr1", 0, 2) == "NAC"              # padded before start
        assert ref.fetch("chr2", 11, 13) == "GGN"            # padded after end
        assert ref.fetch("chr9", 1, 3) == ""


def test_catalog_counts_and_strand_folding(genome, tmp_path):
    # pos 3: G (context C G T) G>A folds to C>T with context A[C>T]G
    # pos 4: T (context G T T) T>C stays T>C, context G[T>C]T
    rows = [
        (T1, "chr1", 3, "G", "A", "SNP", "GRCh38", "muse;mutect2"),
        (T1, "chr1", 4, "T", "C", "SNP", "GRCh38", "mutect2;varscan2"),
        (T1, "chr1", 4, "T", "C", "SNP", "GRCh38", "mutect2"),           # duplicate
        (T1, "chr1", 6, "G", "T", "SNP", "GRCh38", "muse"),              # not MuTect2
        (T1, "chr1", 7, "CA", "TT", "DNP", "GRCh38", "mutect2"),         # not an SNV
        (NORMAL, "chr1", 8, "A", "G", "SNP", "GRCh38", "mutect2"),       # normal sample
    ]
    maf = write_maf(tmp_path / "a.maf", rows)
    with FastaReference(genome) as ref:
        cat, qc = preprocess_mafs([maf], ref)
    assert cat.shape == (1, 96)
    assert cat.loc[T1, "A[C>T]G"] == 1 and cat.loc[T1, "G[T>C]T"] == 1
    assert cat.values.sum() == 2
    assert qc["after_snv_filter"] == 5 and qc["after_caller_filter"] == 4
    assert qc["after_sample_type_filter"] == 3 and qc["after_dedup"] == 2


def test_wrong_genome_build_is_caught(genome, tmp_path):
    rows = [(T1, "chr1", p, "C", "T", "SNP", "GRCh37", "mutect2") for p in (1, 3, 4, 5, 9)]
    maf = write_maf(tmp_path / "b.maf", rows)  # reference bases are mostly not C
    with FastaReference(genome) as ref, pytest.raises(GenomeBuildError):
        preprocess_mafs([maf], ref)


def test_missing_caller_column_is_explicit(genome, tmp_path):
    maf = write_maf(tmp_path / "c.maf", [(T1, "chr1", 3, "G", "A", "SNP", "GRCh38", "x")], caller_col="other")
    with FastaReference(genome) as ref:
        with pytest.raises(ValueError, match="caller"):
            preprocess_mafs([maf], ref)
        cat, _ = preprocess_mafs([maf], ref, caller=None)
        assert cat.values.sum() == 1


def test_mc3_style_caller_column(genome, tmp_path):
    maf = write_maf(tmp_path / "d.maf", [(T1, "1", 3, "G", "A", "SNP", "GRCh37", "MUTECT|MUSE")],
                    caller_col="CENTERS")
    with FastaReference(genome) as ref:
        cat, _ = preprocess_mafs([maf], ref)
    assert cat.values.sum() == 1


def test_revcomp():
    assert revcomp("ACG") == "CGT"


def test_clinical_from_gdc_json(tmp_path):
    hits = [
        {"submitter_id": "TCGA-05-0001", "project": {"project_id": "TCGA-LUAD"},
         "exposures": [{"tobacco_smoking_status": "Current Smoker", "pack_years_smoked": 40}]},
        {"submitter_id": "TCGA-05-0002", "project": {"project_id": "TCGA-LUAD"},
         "exposures": [{"tobacco_smoking_status": "Lifelong Non-Smoker", "pack_years_smoked": None}]},
        {"submitter_id": "TCGA-05-0003", "project": {"project_id": "TCGA-LUAD"},
         "exposures": [{"tobacco_smoking_status": "7"}]},
        {"submitter_id": "TCGA-05-0004", "project": {"project_id": "TCGA-LUAD"}, "exposures": []},
    ]
    p = tmp_path / "cases.json"
    p.write_text(json.dumps(hits))
    clin = load_clinical(p)
    assert clin["smoker"].tolist()[:2] == [1.0, 0.0] and clin["smoker"].iloc[2:].isna().all()
    assert clin.loc["TCGA-05-0001", "pack_years"] == 40


def test_stronger_methods_recover_exposures():
    sigs = demo_signatures()
    S = signature_subset(sigs, "LUAD")
    stats = {"signatures": list(S.columns), "alpha": [0.5] * S.shape[1], "log10_n_mean": 3, "log10_n_sd": 0.1}
    rng = np.random.default_rng(3)
    val = simulate(stats, sigs, 30, rng, n_mutations=500, overdispersion=0, perturb=0)
    test = simulate(stats, sigs, 30, rng, n_mutations=2000, overdispersion=0, perturb=0)
    for m in [M0NNLS(S.to_numpy()), M1ElasticNetNNLS(S.to_numpy()), M2MultinomialEM(S.to_numpy())]:
        pred = m.fit(None, val=val).predict(test["counts"])
        assert np.allclose(pred.sum(axis=1), 1)
        assert np.abs(pred - test["exposures"]).mean() < 0.05, m.name


def test_m2_bootstrap_interval_contains_estimate():
    sigs = demo_signatures()
    S = signature_subset(sigs, "SKCM").to_numpy()
    counts = np.random.default_rng(0).multinomial(800, S @ np.full(S.shape[1], 1 / S.shape[1]))[None, :]
    m = M2MultinomialEM(S)
    est = m.predict(counts)
    lo, hi = m.predict_interval(counts, n_boot=30)
    assert (lo <= est + 1e-6).all() and (est <= hi + 1e-6).all()
