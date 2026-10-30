"""Our attribution methods share one interface, so the same evaluation code
runs on all of them (Q5):

    method.fit(train, val=None)   train/val = dicts with "counts", "exposures"
    method.predict(counts)        -> exposures (tumors x K), rows sum to 1

B0  trivial baseline (Q4)            M0  NNLS baseline
M1  elastic-net NNLS                 M2  multinomial maximum likelihood (EM)
Our M3 (small MLP) will plug in with the same two methods.
None of these is Bayesian: each returns one best estimate. M2 adds
bootstrap intervals, which are frequentist.
"""
import warnings

import numpy as np
from scipy.optimize import nnls
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet


def _normalize(x):
    """We rescale weights to sum to 1; if all are zero we spread them evenly."""
    s = x.sum()
    return x / s if s > 0 else np.full_like(x, 1.0 / len(x))


def _unit_columns(S):
    """We make each signature column a probability profile over the 96 contexts.
    The tiny constant keeps every context possible, avoiding log(0) and 0/0."""
    S = np.asarray(S, dtype=float) + 1e-12
    return S / S.sum(axis=0, keepdims=True)


class B0CohortMean:
    """We ignore the tumor and predict the average training exposure profile.
    This is the attribution analogue of "always predict the majority class"."""
    name = "B0_cohort_mean"

    def fit(self, train, val=None):
        self.mean_ = _normalize(np.asarray(train["exposures"], dtype=float).mean(axis=0))
        return self

    def predict(self, counts):
        # Every tumor gets the same answer, whatever its mutations look like.
        return np.tile(self.mean_, (len(counts), 1))


class M0NNLS:
    """Non-negative least squares: we find the non-negative mix of signatures
    that best reproduces the tumor's mutation profile (squared error)."""
    name = "M0_nnls"

    def __init__(self, signature_matrix):
        self.S = _unit_columns(signature_matrix)

    def fit(self, train, val=None):
        return self  # nothing to learn

    def predict(self, counts):
        counts = np.asarray(counts, dtype=float)
        out = np.empty((len(counts), self.S.shape[1]))
        for i, c in enumerate(counts):
            # A tumor with no mutations carries no information: we return uniform weights.
            if c.sum() == 0:
                out[i] = 1.0 / self.S.shape[1]
                continue
            # We fit the profile (counts / total), not raw counts, so weights are shares.
            x, _ = nnls(self.S, c / c.sum())
            out[i] = _normalize(x)
        return out


class M1ElasticNetNNLS:
    """NNLS plus a penalty that pushes small weights to exactly zero (L1)
    and shrinks correlated weights together (L2). We choose the penalty
    strength on the VALIDATION set only, never on test data (Q2)."""
    name = "M1_elastic_net"

    def __init__(self, signature_matrix, alphas=(1e-7, 1e-6, 1e-5, 1e-4), l1_ratio=0.9):
        self.S = _unit_columns(signature_matrix)
        self.alphas, self.l1_ratio = alphas, l1_ratio
        # Until fit() sees validation data, we use the weakest penalty.
        self.alpha_ = alphas[0]

    def _fit_one(self, c, alpha):
        """We fit one tumor with penalty strength alpha."""
        if c.sum() == 0:
            return np.full(self.S.shape[1], 1.0 / self.S.shape[1])
        with warnings.catch_warnings():
            # Tiny penalties converge slowly; the result is still usable, so we silence the warning.
            warnings.simplefilter("ignore", ConvergenceWarning)
            m = ElasticNet(alpha=alpha, l1_ratio=self.l1_ratio, positive=True,
                           fit_intercept=False, max_iter=5000).fit(self.S, c / c.sum())
        return _normalize(np.clip(m.coef_, 0, None))

    def fit(self, train, val=None):
        if val is not None:
            counts = np.asarray(val["counts"], dtype=float)
            truth = np.asarray(val["exposures"], dtype=float)
            # We try every penalty on the validation tumors and keep the lowest error.
            errs = {a: np.abs(np.array([self._fit_one(c, a) for c in counts]) - truth).mean()
                    for a in self.alphas}
            self.alpha_ = min(errs, key=errs.get)
            self.val_mae_ = errs
        return self

    def predict(self, counts):
        return np.array([self._fit_one(c, self.alpha_) for c in np.asarray(counts, dtype=float)])


class M2MultinomialEM:
    """Our maximum-likelihood fit for count data. We assume each mutation comes
    from one signature; EM alternates between (E) estimating which signature
    produced each context's mutations and (M) updating the weights. It suits
    low mutation counts better than squared-error fitting."""
    name = "M2_multinomial_em"

    def __init__(self, signature_matrix, max_iter=1000, tol=1e-7):
        self.S = _unit_columns(signature_matrix)
        self.max_iter, self.tol = max_iter, tol

    def fit(self, train, val=None):
        return self  # nothing to learn; each tumor is fitted on its own

    def _em(self, c):
        """We run EM for one tumor's 96 counts and return its weights."""
        K = self.S.shape[1]
        if c.sum() == 0:
            return np.full(K, 1.0 / K)
        e = np.full(K, 1.0 / K)  # start from equal weights
        for _ in range(self.max_iter):
            mix = self.S * e                         # 96 x K
            # E step: the share of each context's mutations that each signature explains.
            resp = mix / mix.sum(axis=1, keepdims=True)
            # M step: new weight = expected number of mutations from each signature / total.
            new = (c[:, None] * resp).sum(axis=0) / c.sum()
            # We stop once no weight moves by more than tol.
            if np.abs(new - e).max() < self.tol:
                e = new
                break
            e = new
        return e

    def predict(self, counts):
        return np.array([self._em(c) for c in np.asarray(counts, dtype=float)])

    def predict_interval(self, counts, n_boot=100, level=0.95, seed=0):
        """Our bootstrap interval per weight: we resample each tumor's mutations,
        refit, and take percentiles. Returns (low, high) arrays."""
        rng = np.random.default_rng(seed)
        counts = np.asarray(counts, dtype=float)
        lo, hi = np.empty((len(counts), self.S.shape[1])), np.empty((len(counts), self.S.shape[1]))
        q = (100 * (1 - level) / 2, 100 * (1 + level) / 2)  # e.g. 2.5th and 97.5th percentiles
        for i, c in enumerate(counts):
            n = int(c.sum())
            p = c / n if n else np.full(len(c), 1 / len(c))
            # Each bootstrap tumor has the same total, drawn from the observed profile.
            boots = np.array([self._em(rng.multinomial(n, p).astype(float)) for _ in range(n_boot)])
            lo[i], hi[i] = np.percentile(boots, q, axis=0)
        return lo, hi
