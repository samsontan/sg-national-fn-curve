"""Core frequency-severity model for national building-fire F-N curves.

Severity: discrete power law with upper limit (DPLDwUL; Kaneko et al. 2015,
J Mar Sci Technol 20:14-36), pmf p(n) = n^-b / sum_{i=1}^{Nmax} i^-b.
Frequency: Poisson event process. Estimation: MLE with interval likelihoods
for binned sources, profile-likelihood CIs, parametric-bootstrap G-test
goodness of fit (Clauset, Shalizi & Newman 2009 protocol).

Data are read from data/international_severity_bins.csv (repo-relative).
Every quantity in the paper's Sections 4.9-4.15 regenerates from here.
"""
import csv
import io
import os

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA_CSV = os.path.join(REPO, "data", "international_severity_bins.csv")
RESULTS = os.path.join(REPO, "results")


def load_data(path=DATA_CSV):
    """CSV -> {jurisdiction: {bins: [(lo, hi, count)], T, pop, scope, grade}}"""
    data = {}
    for r in csv.DictReader(io.open(path, encoding="utf-8")):
        j = r["jurisdiction"]
        d = data.setdefault(j, dict(bins=[], T=int(r["T_years"]),
                                    pop=float(r["pop_millions"]),
                                    scope=r["scope"], grade=r["provenance_grade"]))
        d["bins"].append((int(r["bin_lo"]), int(r["bin_hi"]), int(r["events"])))
    return data


def bin_prob(b, Nmax, lo, hi):
    hi = min(hi, Nmax)
    if lo > Nmax:
        return 0.0
    i = np.arange(lo, hi + 1, dtype=float)
    return np.sum(i ** (-b)) / np.sum(np.arange(1, Nmax + 1, dtype=float) ** (-b))


def survival(b, Nmax, n):
    """S(n) = P(severity >= n)."""
    if n > Nmax:
        return 0.0
    i = np.arange(1, Nmax + 1, dtype=float)
    w = i ** (-b)
    return w[int(n) - 1:].sum() / w.sum()


def zsum(b, Nmax, s=0):
    """sum_{i=1}^{Nmax} i^(s-b); s=1 gives the numerator of E[N]."""
    i = np.arange(1, Nmax + 1, dtype=float)
    return np.sum(i ** (s - b))


def nll(b, bins, Nmax):
    ll = 0.0
    for lo, hi, c in bins:
        p = bin_prob(b, Nmax, lo, hi)
        if p <= 0:
            return 1e12
        ll += c * np.log(p)
    return -ll


def fit_b(bins, Nmax):
    r = minimize_scalar(nll, bounds=(0.01, 8.0), args=(bins, Nmax), method="bounded")
    return r.x, -r.fun


def profile_ci(bins, Nmax, level=0.95, grid=None):
    if grid is None:
        grid = np.linspace(0.02, 8.0, 3200)
    lls = np.array([-nll(b, bins, Nmax) for b in grid])
    cut = lls.max() - chi2.ppf(level, 1) / 2
    keep = grid[lls >= cut]
    return float(keep.min()), float(keep.max())


def gof_pboot(bins, Nmax, bhat, nsim=2000, seed=1):
    """Parametric-bootstrap G-test on the observed binning."""
    rng = np.random.default_rng(seed)
    n_tot = sum(c for _, _, c in bins)
    probs = np.array([bin_prob(bhat, Nmax, lo, hi) for lo, hi, _ in bins])
    probs = np.append(probs, max(1e-12, 1 - probs.sum()))
    obs = np.array([c for _, _, c in bins] + [0], dtype=float)

    def g(o, e):
        m = (o > 0) & (e > 0)
        return 2 * np.sum(o[m] * np.log(o[m] / e[m]))

    e = probs * n_tot
    g_obs = g(obs, e)
    gs = np.array([g(rng.multinomial(n_tot, probs).astype(float), e) for _ in range(nsim)])
    return float(g_obs), float(np.mean(gs >= g_obs))


def within_regime_bins(data):
    """The common comparative basis: England excluding Grenfell (N<=6),
    China excluding its three events of N>=30; others unchanged."""
    wr = {j: list(d["bins"]) for j, d in data.items()}
    wr["England"] = [b for b in wr["England"] if b[1] <= 6]
    wr["China"] = [b for b in wr["China"] if b[0] < 30]
    return wr


def lrt_common_b(data, names, nmax_by=None):
    """Likelihood-ratio test of one common exponent across the named sets."""
    if nmax_by is None:
        nmax_by = {n: max(hi for _, hi, _ in data[n]["bins"]) for n in names}

    def negsum(b):
        return sum(nll(b, data[n]["bins"], nmax_by[n]) for n in names)

    r = minimize_scalar(negsum, bounds=(0.01, 8.0), method="bounded")
    ll_common = -r.fun
    ll_free = sum(fit_b(data[n]["bins"], nmax_by[n])[1] for n in names)
    stat = 2 * (ll_free - ll_common)
    df = len(names) - 1
    return dict(b_common=float(r.x), lrt=float(stat), df=df, p=float(chi2.sf(stat, df)))
