"""We test the statistical comparisons, M1's tuning and the sparse simulator."""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import wilcoxon

from sigdet.config import PRESENCE_THRESHOLD
from sigdet.methods import M0NNLS, M1ElasticNetNNLS
from sigdet.metrics import compare_methods, holm_correction, mae_active_per_tumor, paired_wilcoxon
from sigdet.simulate import simulate
from sigdet.target import demo_signatures, signature_subset


def test_holm_matches_hand_computation():
    # Sorted: 0.01*3 = 0.03, 0.02*2 = 0.04, 0.04*1 = 0.04 (kept non-decreasing).
    adj = holm_correction([0.04, 0.01, 0.02])
    assert np.allclose(adj, [0.04, 0.03, 0.04])
    assert np.allclose(holm_correction([0.5, 0.9]), [1.0, 1.0])  # capped at 1
    with pytest.raises(ValueError):
        holm_correction([1.5])


def test_paired_wilcoxon_agrees_with_scipy_and_handles_ties():
    rng = np.random.default_rng(0)
    a, b = rng.random(40), rng.random(40)
    res = paired_wilcoxon(a, b)
    assert res["p_value"] == pytest.approx(wilcoxon(a, b).pvalue)
    assert paired_wilcoxon(a, a)["p_value"] == 1.0  # identical methods: nothing to test
    with pytest.raises(ValueError):
        paired_wilcoxon(a, b[:-1])


def test_compare_methods_pairs_tumors_and_corrects():
    rng = np.random.default_rng(1)
    base = rng.random(60)
    rows = []
    for method, shift in [("good", 0.2), ("bad", 0.0)]:
        rows.append(pd.DataFrame({"cancer": "LUAD", "method": method, "tumor": np.arange(60),
                                  "n": np.repeat([100, 1000], 30), "cosine": base + shift,
                                  "mae_active": 1 - base - shift}))
    tests = compare_methods(pd.concat(rows))
    pooled = tests[(tests["n"] == "all") & (tests["metric"] == "cosine")].iloc[0]
    # Methods are paired in alphabetical order, so "bad" is A and scores lower on every tumor.
    assert pooled["method_a"] == "bad" and pooled["median_difference"] == pytest.approx(-0.2)
    assert pooled["a_higher_fraction"] == 0 and pooled["significant"]
    assert set(tests["n"]) == {"all", "100", "1000"}
    assert (tests["p_holm"] >= tests["p_value"]).all()


def test_mae_active_per_tumor_ignores_absent_signatures():
    true = np.array([[0.6, 0.4, 0.0], [1.0, 0.0, 0.0]])
    pred = np.array([[0.5, 0.3, 0.2], [0.8, 0.2, 0.0]])
    assert np.allclose(mae_active_per_tumor(true, pred), [0.1, 0.2])


def test_m1_searches_full_grid_on_validation_only():
    sigs = demo_signatures()
    S = signature_subset(sigs, "LUAD")
    stats = {"signatures": list(S.columns), "alpha": [0.5] * S.shape[1], "log10_n_mean": 2.3, "log10_n_sd": 0.1}
    val = simulate(stats, sigs, 20, np.random.default_rng(2))
    m = M1ElasticNetNNLS(S.to_numpy(), alphas=(1e-7, 1e-5, 1e-3), l1_ratios=(0.1, 0.9)).fit(None, val=val)
    assert len(m.val_mae_) == 6 and (m.alpha_, m.l1_ratio_) in m.val_mae_
    assert m.val_mae_[(m.alpha_, m.l1_ratio_)] == min(m.val_mae_.values())


def test_sparse_simulation_uses_one_to_three_signatures():
    sigs = demo_signatures()
    S = signature_subset(sigs, "COAD")
    stats = {"signatures": list(S.columns), "alpha": [1.0] * S.shape[1], "log10_n_mean": 3, "log10_n_sd": 0.1}
    sim = simulate(stats, sigs, 200, np.random.default_rng(4), n_mutations=500, n_active=(1, 3))
    active = (sim["exposures"] > 0).sum(axis=1)
    assert active.min() == 1 and active.max() == 3
    assert np.allclose(sim["exposures"].sum(axis=1), 1) and (sim["counts"].sum(axis=1) == 500).all()
    # NNLS should still recover sparse tumors well at 500 mutations, judged at the 5% threshold.
    pred = M0NNLS(S.to_numpy()).predict(sim["counts"])
    assert np.nanmean(mae_active_per_tumor(sim["exposures"], pred, PRESENCE_THRESHOLD)) < 0.15
    with pytest.raises(ValueError):
        simulate(stats, sigs, 5, np.random.default_rng(0), n_active=(0, 2))
