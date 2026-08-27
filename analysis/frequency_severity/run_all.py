"""Reproduce every model-based quantity in the paper (Sections 4.9-4.15).

Outputs to results/: fits_per_jurisdiction.json, common_basis.json,
universality_lrt.json, hierarchical_posterior.json, singapore_tail.json,
grenfell_dragon_king.json, catastrophic_pooling.json, derived_indices.json.

Run:  python analysis/frequency_severity/run_all.py
Runtime a few minutes (hierarchical MCMC and bootstrap fit tests dominate).
"""
import io
import json
import os

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2

from fs_model import (DATA_CSV, RESULTS, fit_b, gof_pboot, load_data, lrt_common_b,
                      nll, profile_ci, survival, within_regime_bins, zsum)

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

# ---------------- 3. universality ----------------
lrt = {
    "all_seven": lrt_common_b(data, names),
    "five_western": lrt_common_b(data, ["USA", "Sweden", "Australia", "NewZealand", "England"]),
    "singapore_china": lrt_common_b(data, ["Singapore", "China"]),
    "note": "Reported in the supplement only; at these sample sizes with residual "
            "misspecification the p-values overstate crispness (paper Section 4.11).",
}
save("universality_lrt.json", lrt)

# ---------------- 4. hierarchical pooling ----------------
J = len(names)
obs_nmax = {j: max(hi for _, hi, _ in data[j]["bins"]) for j in names}


def loglik_j(j, b):
    return -nll(b, data[j]["bins"], obs_nmax[j])


def logpost(th):
    logb, mu, lsig = th[:J], th[J], th[J + 1]
    sig = np.exp(lsig)
    lp = -0.5 * ((mu - np.log(3)) / 1.0) ** 2
    lp += -0.5 * (sig / 0.5) ** 2 + lsig
    lp += np.sum(-0.5 * ((logb - mu) / sig) ** 2 - lsig)
    for k in range(J):
        lp += loglik_j(names[k], np.exp(logb[k]))
    return lp


th = np.concatenate([np.log([fits[j]["b_obs_support"] for j in names]), [np.log(3.0), np.log(0.3)]])
cur = logpost(th)
step = np.concatenate([np.full(J, 0.05), [0.08, 0.15]])
step[names.index("Singapore")] = 0.15
step[names.index("NewZealand")] = 0.08
chain = []
for it in range(40000):
    prop = th + rng.normal(0, step)
    lp = logpost(prop)
    if np.log(rng.random()) < lp - cur:
        th, cur = prop, lp
    if it >= 10000 and it % 10 == 0:
        chain.append(th.copy())
chain = np.array(chain)
b_post = np.exp(chain[:, :J])
hier = {names[k]: [round(x, 3) for x in np.percentile(b_post[:, k], [2.5, 50, 97.5])]
        for k in range(J)}
hier["sd_log_b"] = [round(x, 3) for x in np.percentile(np.exp(chain[:, J + 1]), [2.5, 50, 97.5])]
hier["note"] = "Used only as a stated sensitivity (paper Section 3.4)."
save("hierarchical_posterior.json", hier)

# ---------------- 5. Singapore tail scenarios + monitoring ----------------
sg = data["Singapore"]["bins"]
lam_sg = 47 / 14
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
bd = np.exp(rng.normal(np.log(bu), sd, 4000))
s5 = np.array([survival(b, 10 ** 4, 5) for b in bd])
bf = {f"after_{y}yr": round(1 / np.mean((1 - s5) ** (lam_sg * y)), 1) for y in (10, 20, 30)}
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

# ---------------- 7. catastrophic pooling ----------------
expo = {j: data[j]["pop"] * data[j]["T"] for j in names}
tot = sum(expo.values())
k = 4  # China 38, 74, 309; England 71
rate = k / tot
lo = chi2.ppf(0.025, 2 * k) / 2 / tot
hi = chi2.ppf(0.975, 2 * (k + 1)) / 2 / tot
cat = dict(events=k, exposure_M_person_yr=round(tot),
           pooled_rate_per_M_person_yr=float(f"{rate:.2e}"),
           ci=[float(f"{lo:.2e}"), float(f"{hi:.2e}")],
           singapore_RP_yr=round(1 / (rate * 5.65)),
           leave_one_out_RP_yr={j: round(1 / (((k - (3 if j == "China" else 1 if j == "England" else 0))
                                               / (tot - expo[j])) * 5.65))
                                for j in ("China", "England", "USA")
                                if (k - (3 if j == "China" else 1 if j == "England" else 0)) > 0},
           note="Any-cause reference class; size N>=30 as proxy for regime failure; "
                "assumptions per paper Section 4.13.")
save("catastrophic_pooling.json", cat)

# ---------------- 8. derived indices + two-number estimator ----------------
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
