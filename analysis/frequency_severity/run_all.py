"""Reproduce every model-based quantity in the paper (Sections 4.9-4.15).

Outputs to results/: fits_per_jurisdiction.json, common_basis.json,
singapore_tail.json, grenfell_dragon_king.json, derived_indices.json.
The residential-basis comparison and the England bounds are in run_revision.py.

Run:  python analysis/frequency_severity/run_all.py
Runtime under a minute (bootstrap fit tests dominate).
"""
import io
import json
import os

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2

from fs_model import (RESULTS, fit_b, gof_pboot, load_data,
                      profile_ci, survival, within_regime_bins, zsum)

os.makedirs(RESULTS, exist_ok=True)
rng = np.random.default_rng(42)
data = load_data()
names = list(data)


def save(name, obj):
    io.open(os.path.join(RESULTS, name), "w", encoding="utf-8").write(json.dumps(obj, indent=1))
    print("wrote", name)


# ---------------- 1. per-jurisdiction fits (observed support) ----------------
fits = {}
for j, d in data.items():
    nmax = max(hi for _, hi, _ in d["bins"])
    b, ll = fit_b(d["bins"], nmax)
    ci = profile_ci(d["bins"], nmax)
    g, p = gof_pboot(d["bins"], nmax, b)
    n_ev = sum(c for _, _, c in d["bins"])
    fits[j] = dict(n_events=n_ev, T=d["T"], pop_millions=d["pop"],
                   lam_per_million_yr=round(n_ev / d["T"] / d["pop"], 4),
                   b_obs_support=round(b, 3), ci=[round(x, 3) for x in ci],
                   Nmax_obs=nmax, gof_p=round(p, 3),
                   mfi_observed=round(sum(c for lo, _, c in d["bins"] if lo >= 2) / n_ev, 4))
save("fits_per_jurisdiction.json", fits)

# ---------------- 2. common basis (within-regime, untruncated) ----------------
wr = within_regime_bins(data)
common = {}
for j, bins in wr.items():
    b, ll = fit_b(bins, 10 ** 4)
    ci = profile_ci(bins, 10 ** 4, grid=np.linspace(0.5, 8, 3000))
    n_ev = sum(c for _, _, c in bins)
    common[j] = dict(b=round(b, 3), ci=[round(x, 2) for x in ci], n_within_regime=n_ev,
                     mfi_within_regime=round(sum(c for lo, _, c in bins if lo >= 2) / n_ev, 4))
save("common_basis.json", common)

# ---------------- 3-4. (removed in revision) ----------------
# The common-exponent likelihood-ratio tests and the hierarchical pooling of the
# original analysis are not reported: the datasets differ in incident scope,
# treatment of deliberate fires and collection method (paper Sections 3.3, 4.11).
obs_nmax = {j: max(hi for _, hi, _ in data[j]["bins"]) for j in names}

# ---------------- 5. Singapore tail scenarios + monitoring ----------------
sg = data["Singapore"]["bins"]
lam_sg = sum(c for _, _, c in sg) / data["Singapore"]["T"]
scen = {}
for nmax, label in [(8, "Nmax8_large_household"), (17, "Nmax17_dormitory"),
                    (10 ** 4, "untruncated_bounding")]:
    b, _ = fit_b(sg, nmax)
    scen[label] = dict(b_refit=round(b, 3),
                       RP_5plus_yr=round(1 / (lam_sg * survival(b, nmax, 5)), 1),
                       RP_10plus_yr=(round(1 / (lam_sg * survival(b, nmax, 10)), 1)
                                     if survival(b, nmax, 10) > 0 else None),
                       RP_30plus_yr=(round(1 / (lam_sg * survival(b, nmax, 30)), 1)
                                     if survival(b, nmax, 30) > 0 else None))
# monitoring Bayes factors, marginal over posterior of untruncated b
bu, _ = fit_b(sg, 10 ** 4)
ci_u = profile_ci(sg, 10 ** 4, grid=np.linspace(1.2, 6, 2000))
sd = (np.log(ci_u[1]) - np.log(ci_u[0])) / (2 * 1.96)
# Deterministic quadrature over the lognormal posterior (revision: replaces a
# 4,000-draw Monte Carlo estimate whose last digit depended on the random stream).
zq = np.linspace(-6, 6, 2001)
wq = np.exp(-zq ** 2 / 2); wq /= wq.sum()
s5 = np.array([survival(b, 10 ** 4, 5) for b in np.exp(np.log(bu) + sd * zq)])
bf = {f"after_{y}yr": round(1 / float(np.sum(wq * (1 - s5) ** (lam_sg * y))), 1) for y in (10, 20, 30)}
scen["monitoring_bayes_factors_sharp_cap_vs_uncapped"] = bf
save("singapore_tail.json", scen)

# ---------------- 6. Grenfell dragon-king test ----------------
eng_wr = wr["England"]
b_e, _ = fit_b(eng_wr, 10 ** 4)
lam_e = sum(c for _, _, c in data["England"]["bins"]) / data["England"]["T"]
c70 = survival(b_e, 10 ** 4, 70)
n_ord = sum(c for _, _, c in eng_wr)
gren = dict(b_within_regime=round(b_e, 3),
            P_ge70_per_event=float(f"{c70:.2e}"),
            return_period_yr=round(1 / (lam_e * c70)),
            P_any_of_ordinary_events_ge70=round(1 - (1 - c70) ** (n_ord + 1), 5),
            note="Selection-free: the with-Grenfell refit (4.23) still leaves the event "
                 "at p about 6e-4 across the window (paper Section 4.12).")
save("grenfell_dragon_king.json", gren)

# ---------------- 7. (removed in revision) ----------------
# No pooled catastrophic-fire rate is estimated (paper Section 4.13).

# ---------------- 8. derived indices of the single law, and its exponent from the share ----------------
idx = {}
for j in names:
    b_c = common[j]["b"]
    nmax = obs_nmax[j] if j != "England" else 6
    if j == "China":
        nmax = 29
    en_t = zsum(b_c, nmax, 1) / zsum(b_c, nmax, 0)
    en_u = zsum(b_c, 10 ** 4, 1) / zsum(b_c, 10 ** 4, 0)
    idx[j] = dict(b_common=b_c,
                  lam_x_EN_obs_support=round(fits[j]["lam_per_million_yr"] * en_t, 2),
                  truncation_ratio=round(en_u / en_t, 2),
                  F10_upper_bound_per_M_yr=float(f"{fits[j]['lam_per_million_yr'] * survival(b_c, 10**4, 10):.1e}"))


def b_from_mfi(mfi, nmax=10 ** 4):
    return brentq(lambda b: survival(b, nmax, 2) - mfi, 0.05, 12)


idx["estimator_examples"] = {f"MFI_{m}": round(b_from_mfi(m), 2) for m in (0.06, 0.09, 0.13, 0.22, 0.26)}
idx["estimator_se_note"] = ("SE(b) = |db/dMFI| * sqrt(MFI(1-MFI)/n); about 0.37 at MFI=0.10 "
                            "with n=100, about 0.21 at MFI=0.25 with n=100.")
save("derived_indices.json", idx)

print("done.")
