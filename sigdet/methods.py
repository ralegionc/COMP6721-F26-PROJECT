"""Our attribution methods share one interface, so the same evaluation code
runs on all of them (Q5):

    method.fit(train, val=None)   train/val = dicts with "counts", "exposures"
    method.predict(counts)        -> exposures (tumors x K), rows sum to 1

B0  trivial baseline (Q4)            M0  NNLS baseline
M1  elastic-net NNLS                 M2  multinomial maximum likelihood (EM)
M3  small neural network (MLP) trained on simulated tumors; M3_mlp_broad
    is the same network with extra flat-prior training tumors.
None of these is Bayesian: each returns one best estimate. M2 adds
bootstrap intervals, which are frequentist.
"""
import warnings

import numpy as np
from scipy.optimize import nnls
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet

from .config import PRESENCE_THRESHOLD


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


# M1's search grid. The penalty strength spans four orders of magnitude, from
# effectively plain NNLS (1e-7) to strong enough to zero out most signatures
# (1e-3); the L1 share runs from mostly ridge (0.1) to mostly lasso (0.9).
M1_ALPHAS = (1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
M1_L1_RATIOS = (0.1, 0.5, 0.9)


class M1ElasticNetNNLS:
    """NNLS plus a penalty that pushes small weights to exactly zero (L1)
    and shrinks correlated weights together (L2). We choose the penalty
    strength and the L1 share on the VALIDATION set only, never on test
    data (Q2), by the lowest MAE on truly present signatures."""
    name = "M1_elastic_net"

    def __init__(self, signature_matrix, alphas=M1_ALPHAS, l1_ratios=M1_L1_RATIOS):
        self.S = _unit_columns(signature_matrix)
        self.alphas, self.l1_ratios = tuple(alphas), tuple(l1_ratios)
        # Until fit() sees validation data, we use the weakest penalty.
        self.alpha_, self.l1_ratio_ = self.alphas[0], self.l1_ratios[-1]

    def _fit_one(self, c, alpha, l1_ratio):
        """We fit one tumor with penalty strength alpha and L1 share l1_ratio."""
        if c.sum() == 0:
            return np.full(self.S.shape[1], 1.0 / self.S.shape[1])
        with warnings.catch_warnings():
            # Tiny penalties converge slowly; the result is still usable, so we silence the warning.
            warnings.simplefilter("ignore", ConvergenceWarning)
            m = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, positive=True,
                           fit_intercept=False, max_iter=5000).fit(self.S, c / c.sum())
        return _normalize(np.clip(m.coef_, 0, None))

    def fit(self, train, val=None):
        if val is not None:
            counts = np.asarray(val["counts"], dtype=float)
            truth = np.asarray(val["exposures"], dtype=float)
            present = truth >= PRESENCE_THRESHOLD
            # We try every (penalty, L1 share) pair on the validation tumors
            # and keep the lowest error on the signatures that are truly present.
            errs = {}
            for a in self.alphas:
                for r in self.l1_ratios:
                    pred = np.array([self._fit_one(c, a, r) for c in counts])
                    errs[(a, r)] = float(np.abs(pred - truth)[present].mean())
            self.alpha_, self.l1_ratio_ = min(errs, key=errs.get)
            self.val_mae_ = errs
            # A best penalty at the edge of the grid means the grid may be too narrow.
            self.at_grid_edge_ = self.alpha_ in (self.alphas[0], self.alphas[-1])
        return self

    def predict(self, counts):
        return np.array([self._fit_one(c, self.alpha_, self.l1_ratio_)
                         for c in np.asarray(counts, dtype=float)])


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


# ---------------------------------------------------------------------------
# M3: neural network
# ---------------------------------------------------------------------------

def _mlp_features(counts):
    """Our network input: the square root of each context's share of the tumor's
    mutations (96 values), plus log10(mutation count) / 4.

    The square root evens out the noise: a context with few expected mutations
    varies relatively more than a busy one, and sqrt is the standard fix for
    count data. The mutation-count feature tells the network how noisy the
    profile is, so it can learn to be more cautious with low-count tumors."""
    C = np.asarray(counts, dtype=float)
    n = C.sum(axis=1, keepdims=True)
    share = np.divide(C, n, out=np.full_like(C, 1.0 / C.shape[1]), where=n > 0)
    return np.hstack([np.sqrt(share), np.log10(np.maximum(n, 1)) / 4.0])


def _thin(counts, exposures, n_copies, min_count, rng):
    """Data augmentation by thinning: we keep a random subset of each tumor's
    mutations. A random subset of a tumor's mutations comes from the same mix
    of signatures, so the label (true exposures) stays exactly right. This
    turns 1,000 training tumors into many tumors at every mutation count,
    including the low counts where attribution is hard."""
    C, E = [counts], [exposures]
    for _ in range(n_copies):
        Ci = counts.copy()
        for i, c in enumerate(counts):
            n = int(c.sum())
            if n > min_count:
                # New count drawn evenly on a log scale between min_count and n.
                m = int(round(10 ** rng.uniform(np.log10(min_count), np.log10(n))))
                Ci[i] = rng.multivariate_hypergeometric(c.astype(np.int64), m)
        C.append(Ci)
        E.append(exposures)
    return np.vstack(C), np.vstack(E)


class M3MLP:
    """Our neural network (multi-layer perceptron), trained on simulated tumors,
    built on scikit-learn's MLPRegressor.

    Input: 97 features from _mlp_features. Two hidden layers with ReLU.
    Output: one value per signature of the cancer type. MLPRegressor's output
    layer is linear, so values can be negative and need not sum to 1; we set
    negatives to 0 and divide each row by its sum, which gives a valid
    exposure vector (Q1). Loss: squared error against the true exposures.
    Adam optimizer; alpha is the L2 penalty on the weights (weight decay).

    Early stopping uses OUR validation set (Q2). We do not use MLPRegressor's
    early_stopping=True: it holds out a random 10% of the TRAINING rows, and
    with thinning, copies of the same tumor would sit on both sides of that
    split. Instead we train one epoch at a time with partial_fit and keep the
    model with the lowest validation error.

    Training data: the simulated training tumors, thinned to many mutation
    counts (_thin). With broad=True we add tumors simulated with a flat,
    uninformative prior over signatures (same simulator and noise settings),
    so the network cannot rely only on the calibrated cancer-type mix.
    That is our check of the risk that "M3 may only learn the simulator".

    Defaults (hidden=(128, 64), alpha=0.1, n_augment=8) were chosen on
    validation error only, by tune_m3.py (results/m3_tuning_real.csv).
    Validation error is flat between alpha 0.03 and 0.1 (within 0.1%).
    """

    def __init__(self, signatures, hidden=(128, 64), n_augment=8, min_count=10,
                 broad=False, n_broad=4000, broad_alpha=0.5, lr=1e-3, alpha=0.1,
                 batch=128, max_epochs=300, patience=15, seed=None):
        from .config import SEED
        self.signatures = signatures  # DataFrame, 96 contexts x K signatures
        self.hidden, self.n_augment, self.min_count = tuple(hidden), n_augment, min_count
        self.broad, self.n_broad, self.broad_alpha = broad, n_broad, broad_alpha
        self.lr, self.alpha, self.batch = lr, alpha, batch
        self.max_epochs, self.patience = max_epochs, patience
        self.seed = SEED if seed is None else seed
        self.name = "M3_mlp_broad" if broad else "M3_mlp"

    # ---- building the training data ---------------------------------------
    def _broad_tumors(self, train, n, rng):
        """Extra training tumors with a flat Dirichlet prior over signatures,
        from our own simulator, with mutation counts like the training set."""
        from .simulate import simulate
        names = list(self.signatures.columns)
        logn = np.log10(np.maximum(np.asarray(train["counts"]).sum(axis=1), 1))
        stats = {"signatures": names, "alpha": [self.broad_alpha] * len(names),
                 "log10_n_mean": float(logn.mean()), "log10_n_sd": float(logn.std(ddof=1))}
        sim = simulate(stats, self.signatures, n, rng)
        return sim["counts"], sim["exposures"]

    def _data(self, part, rng, n_broad):
        C = np.asarray(part["counts"], dtype=np.int64)
        E = np.asarray(part["exposures"], dtype=float)
        if self.broad:
            extra = self._broad_tumors(part, n_broad, rng)
            C, E = np.vstack([C, extra[0]]), np.vstack([E, extra[1]])
        C, E = _thin(C, E, self.n_augment, self.min_count, rng)
        return _mlp_features(C), E

    @staticmethod
    def _to_exposures(raw):
        """Negatives to 0, then each row divided by its sum (uniform if all were negative)."""
        P = np.clip(raw, 0, None)
        s = P.sum(axis=1, keepdims=True)
        return np.divide(P, s, out=np.full_like(P, 1.0 / P.shape[1]), where=s > 0)

    # ---- training ------------------------------------------------------------
    def fit(self, train, val=None):
        import copy
        from sklearn.neural_network import MLPRegressor
        rng = np.random.default_rng([self.seed, 3, int(self.broad)])  # own stream per variant
        X, Y = self._data(train, rng, self.n_broad)
        if val is not None:
            # We thin the validation tumors too, so early stopping also checks low counts.
            Xv, Yv = self._data(val, rng, max(200, self.n_broad // 5))
        else:
            Xv, Yv = X, Y

        net = MLPRegressor(hidden_layer_sizes=self.hidden, activation="relu", solver="adam",
                           alpha=self.alpha, batch_size=self.batch, learning_rate_init=self.lr,
                           shuffle=True, random_state=self.seed + int(self.broad))
        best, bad = np.inf, 0
        self.history_ = []
        sub = rng.choice(len(X), size=min(len(X), len(Xv)), replace=False)  # fixed rows to track training error
        for epoch in range(self.max_epochs):
            net.partial_fit(X, Y)  # one pass over the training data
            val_err = float(((net.predict(Xv) - Yv) ** 2).mean())
            self.history_.append({"epoch": epoch + 1, "val_loss": val_err,
                                  "train_loss": float(((net.predict(X[sub]) - Y[sub]) ** 2).mean())})
            # Early stopping: we keep the model with the lowest validation error.
            if val_err < best - 1e-7:
                best, bad = val_err, 0
                self.net_ = copy.deepcopy(net)
            else:
                bad += 1
                if bad >= self.patience:
                    break
        self.best_val_loss_ = best
        return self

    def predict(self, counts):
        return self._to_exposures(self.net_.predict(_mlp_features(counts)))
