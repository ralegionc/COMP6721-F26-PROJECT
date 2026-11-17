"""Tests for M3 (neural network) and the evaluation helpers. Run with:  pytest -q"""
import numpy as np

from sigdet.config import CANCER_SIGNATURES
from sigdet.evaluation import (detection_report, detection_threshold, inverted_prior,
                               limit_of_detection, simulate_unseen)
from sigdet.methods import M3MLP, _thin
from sigdet.simulate import simulate
from sigdet.target import demo_signatures, signature_subset, validate_exposures


def _sets(cancer="LUAD", n_train=300, n_val=80, n_test=60):
    sigs = demo_signatures()
    S = signature_subset(sigs, cancer)
    K = S.shape[1]
    stats = {"signatures": list(S.columns), "alpha": list(np.linspace(0.5, 3, K)),
             "log10_n_mean": 2.3, "log10_n_sd": 0.4}
    rng = np.random.default_rng(0)
    parts = {p: simulate(stats, sigs, n, rng) for p, n in [("train", n_train), ("val", n_val), ("test", n_test)]}
    return sigs, S, stats, parts


def _small_m3(S, **kw):
    return M3MLP(S, hidden=(32, 16), n_augment=2, max_epochs=30, patience=5, **kw)


def test_m3_outputs_valid_exposures():
    _, S, _, parts = _sets()
    m = _small_m3(S).fit(parts["train"], parts["val"])
    pred = m.predict(parts["test"]["counts"])
    validate_exposures(pred, list(S.columns))  # non-negative, rows sum to 1 (Q1)
    assert pred.shape == (60, S.shape[1])
    # A tumor with no mutations still gets a valid answer.
    validate_exposures(m.predict(np.zeros((1, 96))))


def test_m3_learns_something():
    """M3 must beat a constant guess on held-out simulated tumors."""
    _, S, _, parts = _sets()
    m = _small_m3(S).fit(parts["train"], parts["val"])
    truth = parts["test"]["exposures"]
    err_m3 = np.abs(m.predict(parts["test"]["counts"]) - truth).mean()
    err_mean = np.abs(parts["train"]["exposures"].mean(0) - truth).mean()
    assert err_m3 < err_mean


def test_m3_is_reproducible_and_broad_differs():
    _, S, _, parts = _sets()
    a = _small_m3(S).fit(parts["train"], parts["val"]).predict(parts["test"]["counts"])
    b = _small_m3(S).fit(parts["train"], parts["val"]).predict(parts["test"]["counts"])
    np.testing.assert_allclose(a, b, atol=1e-12)  # same seed, same model (Q5)
    m = _small_m3(S, broad=True, n_broad=200)
    assert m.name == "M3_mlp_broad"
    c = m.fit(parts["train"], parts["val"]).predict(parts["test"]["counts"])
    assert not np.allclose(a, c)


def test_thinning_keeps_labels_and_lowers_counts():
    rng = np.random.default_rng(1)
    counts = rng.multinomial(500, np.full(96, 1 / 96), size=5)
    E = rng.dirichlet(np.ones(3), size=5)
    C2, E2 = _thin(counts, E, n_copies=3, min_count=10, rng=rng)
    assert C2.shape == (20, 96) and E2.shape == (20, 3)
    np.testing.assert_array_equal(E2[5:10], E)            # labels unchanged
    assert (C2[5:] <= np.tile(counts, (3, 1))).all()       # a subset of the original mutations
    assert (C2[5:].sum(1) >= 10).all() and (C2[5:].sum(1) <= 500).all()


def test_detection_constant_predictor_scores_half():
    """Per-signature ROC-AUC of a constant prediction is 0.5, whatever the prevalence."""
    rng = np.random.default_rng(2)
    E = rng.dirichlet([3.0, 1.0, 0.3], size=300)
    const = np.tile(E.mean(0), (300, 1))
    rep = detection_report(E, const, call_threshold=0.05)
    assert abs(rep["roc_auc_macro"] - 0.5) < 1e-12
    perfect = detection_report(E, E, call_threshold=0.05 - 1e-9)
    assert perfect["roc_auc_macro"] == 1.0 and perfect["hit_rate_at_fa5"] == 1.0


def test_detection_threshold_gives_requested_false_alarms():
    rng = np.random.default_rng(3)
    E = rng.dirichlet([1.0, 1.0, 0.2, 0.2], size=2000)
    noisy = np.clip(E + rng.normal(0, 0.05, E.shape), 0, None)
    thr = detection_threshold(E, noisy, false_alarm=0.05)
    absent = E < 0.05
    assert abs((noisy[absent] > thr).mean() - 0.05) < 0.01


def test_limit_of_detection_bins_cover_present_pairs():
    rng = np.random.default_rng(4)
    E = rng.dirichlet([1.0, 1.0, 1.0], size=100)
    n = np.full(100, 200)
    table = limit_of_detection(E, E, n, call_threshold=0.0)
    assert table["n_pairs"].sum() == (E >= 0.05).sum()
    assert (table["hit_rate"] == 1.0).all()


def test_inverted_prior_swaps_common_and_rare():
    stats = {"signatures": ["a", "b", "c"], "alpha": [5.0, 1.0, 0.2]}
    flipped = inverted_prior(stats)["alpha"]
    assert flipped == [0.2, 1.0, 5.0]
    assert stats["alpha"] == [5.0, 1.0, 0.2]  # the original is not changed


def test_unseen_signature_set_has_extra_mutations_and_valid_truth():
    sigs, S, stats, _ = _sets()
    unseen = [s for s in sigs.columns if s not in CANCER_SIGNATURES["LUAD"]][0]
    data = simulate_unseen(stats, sigs, unseen, 0.2, 5, np.random.default_rng(5), n_grid=[100, 1000])
    assert data["counts"].shape == (10, 96)
    np.testing.assert_array_equal(data["counts"].sum(1), data["n_mutations"])
    validate_exposures(data["exposures"], list(S.columns))
