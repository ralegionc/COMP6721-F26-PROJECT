"""We keep one test group per checklist question. Run with:  pytest -q"""
import numpy as np
import pandas as pd
import pytest

from sigdet.config import ARTIFACT_SIGNATURES, CANCER_SIGNATURES
from sigdet.labels import presence_labels, smoking_labels
from sigdet.methods import B0CohortMean, M0NNLS
from sigdet.metrics import breakdown_point, majority_class_check, presence_report
from sigdet.simulate import assert_calibrated_on_train, calibrate, simulate
from sigdet.splits import (LeakageError, assert_no_patient_leak, assign_partition, load_split,
                           make_patient_split, remove_overlapping_donors, save_split)
from sigdet.storage import ensure_dirs, load_sim_set, save_sim_set, verify_manifest, write_manifest
from sigdet.target import demo_signatures, signature_subset, validate_exposures


def _barcodes(n=40, tss="05"):
    out, cancers = [], []
    for i in range(n):
        for s in ["01A", "06A"] if i % 3 == 0 else ["01A"]:
            out.append(f"TCGA-{tss}-{i:04d}-{s}-01D-A000-08")
            cancers.append("LUAD")
    return out, cancers


# ---- Q1 -------------------------------------------------------------------
def test_q1_exposures_must_be_valid_shares():
    validate_exposures(np.array([[0.7, 0.3], [0.5, 0.5]]))
    with pytest.raises(ValueError):
        validate_exposures(np.array([[0.7, 0.4]]))      # sums to 1.1
    with pytest.raises(ValueError):
        validate_exposures(np.array([[1.2, -0.2]]))     # negative share


def test_q1_subsets_exclude_artifacts():
    sigs = demo_signatures()
    for cancer in CANCER_SIGNATURES:
        assert not set(signature_subset(sigs, cancer).columns) & ARTIFACT_SIGNATURES


# ---- Q2 -------------------------------------------------------------------
def test_q2_split_is_by_patient():
    bcs, cancers = _barcodes()
    split = make_patient_split(bcs, cancers)
    assert_no_patient_leak(split)
    part = pd.Series(assign_partition(bcs, split), index=[b[:12] for b in bcs])
    assert (part.groupby(level=0).nunique() == 1).all()


def test_q2_leak_is_detected():
    with pytest.raises(LeakageError):
        assert_no_patient_leak({"train": ["TCGA-05-0001"], "val": [], "test": ["TCGA-05-0001"]})


def test_q2_pcawg_overlap_removed():
    keep = remove_overlapping_donors(["TCGA-05-0001-01A", "DO12345"], ["TCGA-05-0001"])
    assert keep.tolist() == [False, True]


def test_q2_simulator_sees_only_training_patients():
    bcs, cancers = _barcodes()
    split = make_patient_split(bcs, cancers)
    E = pd.DataFrame(np.random.default_rng(0).dirichlet([1, 1, 1], len(bcs)), index=bcs, columns=list("abc"))
    stats = calibrate(E, pd.Series(500, index=bcs), split)
    assert set(stats["train_patients"]) <= set(split["train"])
    assert_calibrated_on_train(stats, split)


# ---- Q3 -------------------------------------------------------------------
def test_q3_presence_threshold_is_inclusive():
    assert presence_labels([[0.05, 0.0499, 0.9]]).tolist() == [[1, 0, 1]]


def test_q3_smoking_classes():
    s = smoking_labels(["Lifelong Non-Smoker", "Current Smoker", 4, "Not Reported", None])
    assert s.tolist()[:3] == [0.0, 1.0, 1.0] and s.iloc[3:].isna().all()
    long = smoking_labels(["Current Reformed Smoker for > 15 yrs"], exclude_long_reformed=True)
    assert long.isna().all()


# ---- Q4 -------------------------------------------------------------------
def test_q4_always_absent_looks_accurate_but_is_useless():
    E = np.zeros((100, 10))
    E[:, 0] = 1.0                          # one signature present, nine absent
    check = majority_class_check(E)
    assert check["always_absent_accuracy"] == pytest.approx(0.9)
    assert check["always_absent_f1"] == 0.0 and check["accuracy_would_overstate"]
    _, summ = presence_report(E, np.zeros_like(E), list("abcdefghij"))
    assert summ["accuracy"] == pytest.approx(0.9) and summ["micro_f1"] == 0.0


def test_q4_nnls_beats_b0_on_clean_data():
    sigs = demo_signatures()
    S = signature_subset(sigs, "LUAD")
    stats = {"signatures": list(S.columns), "alpha": [0.5] * S.shape[1], "log10_n_mean": 3, "log10_n_sd": 0.1}
    rng = np.random.default_rng(1)
    train = simulate(stats, sigs, 200, rng, n_mutations=2000)
    test = simulate(stats, sigs, 200, rng, n_mutations=2000)
    b0 = B0CohortMean().fit(train).predict(test["counts"])
    m0 = M0NNLS(S.to_numpy()).predict(test["counts"])
    err = lambda p: np.abs(p - test["exposures"]).mean()
    assert err(m0) < err(b0)


def test_q4_breakdown_point():
    df = pd.DataFrame({"method": "m", "n": [10, 10, 100, 100, 1000, 1000],
                       "cosine": [0.5, 0.6, 0.95, 0.91, 0.97, 0.99]})
    assert breakdown_point(df)["m"] == 100


# ---- Q5 -------------------------------------------------------------------
def test_q5_split_is_reproducible_saved_once_and_tamper_checked(tmp_path):
    bcs, cancers = _barcodes()
    a, b = make_patient_split(bcs, cancers), make_patient_split(bcs, cancers)
    assert a == b
    path = tmp_path / "split.json"
    save_split(a, path)
    with pytest.raises(FileExistsError):
        save_split(a, path)
    loaded, _ = load_split(path)
    assert loaded == a
    path.write_text(path.read_text().replace(a["test"][0], "TCGA-ZZ-9999"))
    with pytest.raises(ValueError):
        load_split(path)


# ---- Q6 -------------------------------------------------------------------
def test_q6_layout_manifest_and_sim_checksums(tmp_path):
    paths = ensure_dirs(tmp_path)
    for key in ["raw", "processed", "sims", "splits", "results", "figures"]:
        assert paths[key].is_dir()
    assert "data/raw/" in (tmp_path / ".gitignore").read_text()

    f = paths["raw"] / "0a1b2c3d-1111-2222-3333-444455556666" / "x.maf"
    f.parent.mkdir()
    f.write_text("mutation data")
    df = write_manifest(paths["raw"], paths["manifest"])
    assert df["gdc_file_uuid"].iloc[0] == "0a1b2c3d-1111-2222-3333-444455556666"
    assert verify_manifest(paths["manifest"], paths["raw"]) == []
    f.write_text("edited")
    assert verify_manifest(paths["manifest"], paths["raw"]) == ["changed: 0a1b2c3d-1111-2222-3333-444455556666/x.maf"]

    sim = paths["sims"] / "t.npz"
    save_sim_set(sim, {"x": 1}, counts=np.ones((2, 96)))
    data, meta = load_sim_set(sim)
    assert meta == {"x": 1} and data["counts"].shape == (2, 96)
