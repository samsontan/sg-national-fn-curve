"""Revision analysis (R1): comparisons restricted to stated, consistent bases.

Reads data/revision_bases.csv. Writes results/revision_bases.json.
Run:  python analysis/frequency_severity/run_revision.py

Bases
  residential : fatal fires in residential buildings / dwellings
  own scope   : each source on the scope it publishes
England's source resolves deaths per incident only for all fires attended by
fire and rescue services. Its building-fire and dwelling-fire severity are
therefore BOUNDED, not estimated: the upper bound assigns every multi-fatality
fire in the source table to the subset.
"""
import csv
import io
import json
import os

import numpy as np
from scipy.stats import beta, fisher_exact, chi2_contingency

from fs_model import REPO, RESULTS, fit_b, gof_pboot, profile_ci, survival

CSV = os.path.join(REPO, "data", "revision_bases.csv")


def load():
    d = {}
    for r in csv.DictReader(io.open(CSV, encoding="utf-8")):
        k = r["dataset"]
        e = d.setdefault(k, dict(bins=[], T=float(r["T_years"]), pop=float(r["pop_millions"]),
                                 basis=r["basis"], scope=r["scope"]))
        e["bins"].append((int(r["bin_lo"]), int(r["bin_hi"]), int(r["events"])))
    return d


def cp(k, n, level=0.95):
    a = (1 - level) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


data = load()
out = {}
for k, d in data.items():
    bins = d["bins"]
    n = sum(c for _, _, c in bins)
    multi = sum(c for lo, _, c in bins if lo >= 2)
    b, _ = fit_b(bins, 10 ** 4)
    ci = profile_ci(bins, 10 ** 4, grid=np.linspace(0.5, 8, 3000))
    nmax = max(hi for _, hi, _ in bins)
    b_obs, _ = fit_b(bins, nmax)
    g, p = gof_pboot(bins, nmax, b_obs)
    lo, hi = cp(multi, n)
    out[k] = dict(basis=d["basis"], scope=d["scope"], events=n, T=d["T"], pop_millions=d["pop"],
                  rate_per_M_yr=round(n / d["T"] / d["pop"], 3),
                  multi_fatality_events=multi, mfi=round(multi / n, 4),
                  mfi_ci95=[round(lo, 3), round(hi, 3)],
                  b_untruncated=round(b, 2), b_ci95=[round(x, 2) for x in ci],
                  gof_p_obs_support=round(p, 3))

# ---- England bounds (source resolves severity for all fires only) ----
ENG_T, ENG_POP = 9.0, 54.35
eng = dict(all_fires=2356, dwellings=1808, other_buildings=144, road_vehicles=247, outdoor=157,
           multi_all=143, shoreham_outdoor=1, grenfell_dwelling=1)
multi_bldg_max = eng["multi_all"] - eng["shoreham_outdoor"]          # 142, incl. Grenfell
bounds = dict(
    dwelling_rate_per_M_yr=round(eng["dwellings"] / ENG_T / ENG_POP, 2),
    building_rate_per_M_yr=round((eng["dwellings"] + eng["other_buildings"]) / ENG_T / ENG_POP, 2),
    all_fires_rate_per_M_yr=round(eng["all_fires"] / ENG_T / ENG_POP, 2),
    mfi_all_fires=round(eng["multi_all"] / eng["all_fires"], 4),
    mfi_dwelling_upper_bound=round(multi_bldg_max / eng["dwellings"], 4),
    mfi_building_upper_bound=round(multi_bldg_max / (eng["dwellings"] + eng["other_buildings"]), 4),
)
sg_res = out["Singapore_residential"]
sg_all = out["Singapore_buildings"]
bounds["ratio_SGres_over_ENGdwelling_min"] = round(sg_res["mfi"] / bounds["mfi_dwelling_upper_bound"], 2)
bounds["ratio_SGbldg_over_ENGbuilding_min"] = round(sg_all["mfi"] / bounds["mfi_building_upper_bound"], 2)
bounds["rate_ratio_ENGdwelling_over_SGres"] = round(bounds["dwelling_rate_per_M_yr"] / sg_res["rate_per_M_yr"], 1)
bounds["rate_ratio_ENGbuilding_over_SGbldg"] = round(bounds["building_rate_per_M_yr"] / sg_all["rate_per_M_yr"], 1)
# exact tests of Singapore against the England upper-bound tables
t_dw = [[sg_res["multi_fatality_events"], sg_res["events"] - sg_res["multi_fatality_events"]],
        [multi_bldg_max, eng["dwellings"] - multi_bldg_max]]
t_bl = [[sg_all["multi_fatality_events"], sg_all["events"] - sg_all["multi_fatality_events"]],
        [multi_bldg_max, eng["dwellings"] + eng["other_buildings"] - multi_bldg_max]]
bounds["fisher_p_SGres_vs_ENGdwelling_bound"] = float(f"{fisher_exact(t_dw)[1]:.1e}")
bounds["fisher_p_SGbldg_vs_ENGbuilding_bound"] = float(f"{fisher_exact(t_bl)[1]:.1e}")
out["England_bounds"] = bounds

# ---- Grenfell: all-fires law, and building-only bounding sensitivity ----
def grenfell(bins, n_events_total):
    b, _ = fit_b(bins, 10 ** 4)
    ci = profile_ci(bins, 10 ** 4, grid=np.linspace(0.5, 8, 3000))
    lam = n_events_total / ENG_T
    res = dict(b=round(b, 2), ci=[round(x, 2) for x in ci], events_per_yr=round(lam, 1))
    for lab, bb in (("point", b), ("flat_edge", ci[0])):
        s = survival(bb, 10 ** 4, 70)
        res[f"P_ge70_{lab}"] = float(f"{s:.1e}")
        res[f"RP_yr_{lab}"] = int(round(1 / (lam * s), -2))
    return res

all_wr = [(1, 1, 2213), (2, 2, 108), (3, 3, 17), (4, 4, 8), (5, 5, 5), (6, 6, 3)]
n_b = eng["dwellings"] + eng["other_buildings"]                      # 1,952 incl. Grenfell
bld_wr = [(1, 1, n_b - 1 - 141), (2, 2, 108), (3, 3, 17), (4, 4, 8), (5, 5, 5), (6, 6, 3)]
n_d = eng["dwellings"]                                               # 1,808 incl. Grenfell
dwl_wr = [(1, 1, n_d - 1 - 141), (2, 2, 108), (3, 3, 17), (4, 4, 8), (5, 5, 5), (6, 6, 3)]
out["Grenfell"] = dict(all_fires_law=grenfell(all_wr, 2355),
                       dwelling_only_heaviest_bound=grenfell(dwl_wr, n_d),
                       building_only_heaviest_bound=grenfell(bld_wr, n_b),
                       note="Building-only bound assigns all 141 ordinary multi-fatality fires to "
                            "buildings: the heaviest building-fire severity the source admits.")

# ---- residential-basis contrasts ----
def tab(k):
    return [out[k]["multi_fatality_events"], out[k]["events"] - out[k]["multi_fatality_events"]]

res_keys = ["Singapore_residential", "USA_residential", "Australia_preventable_residential",
            "NewZealand_unintentional_domestic"]
out["residential_contrasts"] = dict(
    SGres_vs_USA_ratio=round(out["Singapore_residential"]["mfi"] / out["USA_residential"]["mfi"], 2),
    SGres_vs_USA_note="USA share is a published rounded percentage (13%, admissible 12.5-13.5%); "
                      "the comparison uses Singapore's exact binomial interval against that value.",
    SGres_mfi_ci95=out["Singapore_residential"]["mfi_ci95"])

# ---- sensitivities of the Singapore-England contrast (paper Section 4.8) ----
from fs_model import singapore_events
_res = singapore_events("residential"); _res18 = singapore_events("residential", (2012, 2018))
SG_N, SG_K = len(_res), sum(1 for r in _res if int(r["fatalities"]) >= 2)            # residential fatal fires, 2012-2025
SG18_N, SG18_K = len(_res18), sum(1 for r in _res18 if int(r["fatalities"]) >= 2)    # the same, 2012-2018
assert max(int(r["fatalities"]) for r in _res) == 4                                 # the largest event is the lodging-house fire
SG_DWELLINGS_MEAN = 1_434_347          # mean residential stock 2012-2025 (data/exposure_national.csv)
ENG_DWELLINGS_MEAN = 23_635_441        # mean of 31 March 2011-2019, MHCLG Live Table 104
bound_d = bounds["mfi_dwelling_upper_bound"]
out["contrast_sensitivities"] = dict(
    ratio_at_SG_lower_limit_residential=round(sg_res["mfi_ci95"][0] / bound_d, 2),
    ratio_at_SG_lower_limit_buildings=round(sg_all["mfi_ci95"][0] / bounds["mfi_building_upper_bound"], 2),
    without_largest_event=dict(k=SG_K - 1, n=SG_N - 1, mfi=round((SG_K - 1) / (SG_N - 1), 3), ci95=[round(x, 3) for x in cp(SG_K - 1, SG_N - 1)]),
    window_2012_2018=dict(k=SG18_K, n=SG18_N, mfi=round(SG18_K / SG18_N, 3), ci95=[round(x, 3) for x in cp(SG18_K, SG18_N)]),
    missed_single_death_fires_to_reach_bound=int(np.ceil(sg_res["multi_fatality_events"] / bound_d - sg_res["events"])),
    per_dwelling_rate_SG=float(f"{sg_res['events'] / sg_res['T'] / SG_DWELLINGS_MEAN:.2e}"),
    per_dwelling_rate_ENG=float(f"{eng['dwellings'] / ENG_T / ENG_DWELLINGS_MEAN:.2e}"),
    per_dwelling_rate_ratio_ENG_over_SG=round((eng["dwellings"] / ENG_T / ENG_DWELLINGS_MEAN)
                                              / (sg_res["events"] / sg_res["T"] / SG_DWELLINGS_MEAN), 1),
    note="The 2012-2018 and without-largest-event counts are read from the register "
         "(residential premises; the event removed is the four-death lodging-house fire, one of two four-death fires).")

# ---- Singapore split-sample and trend checks (paper Section 4.14) ----
from scipy.stats import binomtest, chi2, norm
from scipy.optimize import minimize
from fs_model import nll, zsum

reg = [r for r in csv.DictReader(io.open(os.path.join(REPO, "data", "fatal_incident_register.csv"), encoding="utf-8"))
       if r["building_scope"] == "building" and 2012 <= int(r["year"]) <= 2025]


def bins_of(rows):
    c = {}
    for r in rows:
        c[int(r["fatalities"])] = c.get(int(r["fatalities"]), 0) + 1
    return [(n, n, k) for n, k in sorted(c.items())]


h1 = [r for r in reg if int(r["year"]) <= 2020]; h2 = [r for r in reg if int(r["year"]) >= 2021]
b1, ll1 = fit_b(bins_of(h1), 4); b2, ll2 = fit_b(bins_of(h2), 4); b12, ll12 = fit_b(bins_of(reg), 4)
lrt = 2 * (ll1 + ll2 - ll12)
pred_mean = zsum(b1, 4, 1) / zsum(b1, 4, 0)
real_mean = sum(int(r["fatalities"]) for r in h2) / len(h2)
years = np.arange(2012, 2026); counts = np.array([sum(1 for r in reg if int(r["year"]) == y) for y in years])
x = years - years.mean()


def pois_nll(th):
    mu = np.exp(th[0] + th[1] * x)
    return float(np.sum(mu - counts * np.log(mu)))


fit1 = minimize(pois_nll, [np.log(counts.mean()), 0.0], method="BFGS")
fit0 = float(np.sum(counts.mean() - counts * np.log(counts.mean())))
trend_p = float(chi2.sf(2 * (fit0 - fit1.fun), 1))
bu = out["Singapore_buildings"]["b_untruncated"]; ci_b = out["Singapore_buildings"]["b_ci95"]
lam_sg = len(reg) / 14
out["singapore_checks"] = dict(
    residential_events=SG_N, residential_multi=SG_K, residential_events_2012_2018=SG18_N, residential_multi_2012_2018=SG18_K,
    events_2012_2020=len(h1), events_2021_2025=len(h2),
    b_2012_2020_obs_support=round(b1, 2), b_2021_2025_obs_support=round(b2, 2),
    equality_lrt_p=round(float(chi2.sf(lrt, 1)), 2),
    predicted_mean_deaths_per_event_2021_2025=round(pred_mean, 2), realised_mean=round(real_mean, 2),
    last5_of_14_years_events=int(counts[-5:].sum()),
    binomial_split_p_two_sided=round(float(binomtest(int(counts[-5:].sum()), len(reg), 5 / 14).pvalue), 3),
    loglinear_poisson_trend_p=round(trend_p, 2),
    RP10_untruncated_yr=round(1 / (lam_sg * survival(bu, 10 ** 4, 10))),
    RP10_at_b_interval_edges_yr=[round(1 / (lam_sg * survival(ci_b[0], 10 ** 4, 10))),
                                 round(1 / (lam_sg * survival(ci_b[1], 10 ** 4, 10)))],
    note="Split-sample fits at observed support (Nmax = 4). The five-year split was chosen after "
         "inspecting the data; the log-linear trend test is the unbiased check.")

io.open(os.path.join(RESULTS, "revision_bases.json"), "w", encoding="utf-8").write(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
