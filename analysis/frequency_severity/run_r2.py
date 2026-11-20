"""Second-round revision analysis (R2): comparison by number of deaths on matched exposure.

The first-round comparison used the share of fatal fires with two or more deaths. That share is
the ratio of two rates, so this script compares the rates themselves:

    lambda_1  rate of fatal fires with exactly one death
    lambda_2  rate of fatal fires with two or more deaths
    S_2(n)    P(N >= n | N >= 2), the shape of the curve above one death

Each is computed for Singapore residential fires and for English dwelling fires on stated,
matched windows and denominators, with England re-cut from record-level data to the scope of
each other published dataset. Nothing is pooled across datasets.

Reads   data/fatal_incident_register.csv            Singapore register
        data/home_office_fire_related_fatalities.csv  England, one row per fire-related death
        data/england_fatal_fires_official.csv         England, official counts of fatal fires and deaths
        data/exposure_comparison.csv, data/population_by_age.csv, data/revision_bases.csv
Writes  results/r2_comparison.json
Run:    python analysis/frequency_severity/run_r2.py
"""
import csv
import io
import json
import os
from collections import Counter, defaultdict

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import beta, chi2, fisher_exact
from scipy.stats.contingency import odds_ratio as exact_odds_ratio

from fs_model import REPO, RESULTS

DATA = os.path.join(REPO, "data")
MULTI = "Multiple fatality incident"
KEY = ["FINANCIAL_YEAR", "HOUR_BANDS", "INCIDENT_LOCATION_TYPE", "ACCIDENTAL_OR_DELIBERATE", "ALARM_SYSTEM",
       "CAUSE_OF_FIRE", "SOURCE_OF_IGNITION", "FIRE_START_LOCATION"]
MONTHS = {m: i for i, m in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split(), 1)}


def rd(name, enc="utf-8-sig"):
    return list(csv.DictReader(io.open(os.path.join(DATA, name), encoding=enc)))


# ----------------------------------------------------------------------------- statistics helpers
def cp(k, n, level=0.95):
    """exact (Clopper-Pearson) interval for a binomial proportion"""
    a = (1 - level) / 2
    lo = float(beta.ppf(a, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta.ppf(1 - a, k + 1, n - k)) if k < n else 1.0
    return lo, hi


def pois_ci(k, level=0.95):
    """exact (Garwood) interval for a Poisson count"""
    a = (1 - level) / 2
    lo = float(chi2.ppf(a, 2 * k) / 2) if k > 0 else 0.0
    hi = float(chi2.ppf(1 - a, 2 * (k + 1)) / 2)
    return lo, hi


def rate(k, exposure_millions):
    lo, hi = pois_ci(k)
    return dict(events=int(k), exposure_millions=round(exposure_millions, 3), rate=round(k / exposure_millions, 4),
                ci95=[round(lo / exposure_millions, 4), round(hi / exposure_millions, 4)])


def rate_ratio(k1, e1, k2, e2):
    """(k1/e1) / (k2/e2) with the exact conditional interval (binomial on k1 given k1 + k2)"""
    lo, hi = cp(k1, k1 + k2)
    f = lambda p: (p / (1 - p)) * e2 / e1 if p < 1 else float("inf")
    return dict(ratio=round((k1 / e1) / (k2 / e2), 6), ci95=[round(f(lo), 6), round(f(hi), 6)])


def tail_bins(dist, nmin=2, nmax=None):
    return [(n, c) for n, c in sorted(dist.items()) if n >= nmin and (nmax is None or n <= nmax)]


def fit_tail(bins, n0=2, cap=10 ** 4):
    """discrete power law on n0..cap for exact or interval bins [(lo, hi, count)]; returns b, loglik, profile CI"""
    i = np.arange(n0, cap + 1, dtype=float)

    def nll(b):
        w = i ** (-b); z = w.sum(); ll = 0.0
        for lo, hi, c in bins:
            p = w[int(lo) - n0:int(min(hi, cap)) - n0 + 1].sum() / z
            ll += c * np.log(p)
        return -ll

    r = minimize_scalar(nll, bounds=(1.05, 9.0), method="bounded")
    grid = np.linspace(1.05, 9.0, 1600)
    lls = np.array([-nll(b) for b in grid])
    keep = grid[lls >= -r.fun - chi2.ppf(0.95, 1) / 2]
    return float(r.x), float(-r.fun), [float(keep.min()), float(keep.max())]


def tail_survival(b, n, n0=2, cap=10 ** 4):
    if n > cap:
        return 0.0
    i = np.arange(n0, cap + 1, dtype=float); w = i ** (-b)
    return float(w[int(n) - n0:].sum() / w.sum())


def exact(dist):
    return [(n, n, c) for n, c in sorted(dist.items()) if n >= 2]


# ----------------------------------------------------------------------------- exposure
EXPO = defaultdict(dict)
for r in rd("exposure_comparison.csv"):
    EXPO[(r["jurisdiction"], r["measure"])][int(r["year"])] = float(r["value"])


def mean_of(jur, measure, years):
    s = EXPO[(jur, measure)]
    return float(np.mean([s[y] for y in years]))


def eng_fy_dwellings(fy):
    """dwelling stock at mid financial year: mean of the 31 March counts at its start and end"""
    y = int(fy[:4]); s = EXPO[("England", "dwellings")]
    return (s[y] + s[y + 1]) / 2 if (y + 1) in s else s[y]


def eng_fy_population(fy):
    return EXPO[("England", "population_total")][int(fy[:4])]      # mid-year estimate inside the financial year


def eng_households(year):
    """households: linear between the 2011 and 2021 censuses, extended at the same slope"""
    s = EXPO[("England", "households")]
    return s[2011] + (s[2021] - s[2011]) * (year - 2011) / 10.0


def eng_one_person(year):
    s = EXPO[("England", "households_one_person")]
    return s[2011] + (s[2021] - s[2011]) * (year - 2011) / 10.0


# ----------------------------------------------------------------------------- Singapore register
REG = [r for r in rd("fatal_incident_register.csv") if r["building_scope"] == "building"]
for r in REG:
    r["y"], r["m"], r["n"] = int(r["year"]), MONTHS.get(r["month"][:3], 6), int(r["fatalities"])
NOT_IN_STOCK = {"SGF-2014-04", "SGF-2015-01"}       # commercial shophouse lodging; converted workers' accommodation
COMMON_AREA_WORKERS = {"SGF-2014-01"}               # two security officers, lift lobby of a condominium


def is_res(r):
    p = r["premises_type"].lower()
    return ("hdb" in p or "residen" in p or "condo" in p) or r["incident_id"] == "SGF-2014-04"


def sg_sel(rows, basis):
    if basis == "buildings":
        return rows
    out = [r for r in rows if is_res(r)]
    if basis in ("dwelling_stock", "household_units"):
        out = [r for r in out if r["incident_id"] not in NOT_IN_STOCK]
    if basis == "household_units":
        out = [r for r in out if r["incident_id"] not in COMMON_AREA_WORKERS]
    return out


def dist_of(rows):
    return dict(sorted(Counter(r["n"] for r in rows).items()))


def counts(dist):
    n = sum(dist.values()); one = dist.get(1, 0)
    return dict(fatal_fires=n, single=one, multi=n - one, ge3=sum(c for k, c in dist.items() if k >= 3),
                ge4=sum(c for k, c in dist.items() if k >= 4), ge5=sum(c for k, c in dist.items() if k >= 5),
                deaths=sum(k * c for k, c in dist.items()), deaths_in_multi=sum(k * c for k, c in dist.items() if k >= 2),
                distribution={str(k): v for k, v in dist.items()})


SG_W = dict(
    primary=dict(label="2012-2025", rows=[r for r in REG if 2012 <= r["y"] <= 2025], years=list(range(2012, 2026))),
    overlap=dict(label="April 2012 to March 2024", years=list(range(2012, 2024)),
                 rows=[r for r in REG if (2012, 4) <= (r["y"], r["m"]) <= (2024, 3)]),
    review=dict(label="2012-2018", rows=[r for r in REG if 2012 <= r["y"] <= 2018], years=list(range(2012, 2019))),
)

# ----------------------------------------------------------------------------- England
HO = [{k: (v or "").strip() for k, v in r.items()} for r in rd("home_office_fire_related_fatalities.csv", "utf-8")]
OFF = rd("england_fatal_fires_official.csv")
FYS = sorted({r["FINANCIAL_YEAR"] for r in HO})
ENG_W = dict(primary=FYS[2:16], overlap=FYS[2:14], review=FYS[0:9], all=FYS)
assert ENG_W["primary"][0] == "2012/13" and ENG_W["primary"][-1] == "2025/26" and ENG_W["review"][-1] == "2018/19"


def is_grenfell(r):
    return (r["FINANCIAL_YEAR"] == "2017/18" and r["HOUR_BANDS"] == "00-01" and r["FATALITY_TOTAL"] == MULTI
            and r["INCIDENT_LOCATION_TYPE"] == "Dwellings" and r["FIRE_START_LOCATION"] == "Kitchen"
            and r["ALARM_SYSTEM"] == "Alarm Present but did not raise alarm")


assert sum(1 for r in HO if is_grenfell(r)) == 71


def eng_dist(sel):
    """deaths-per-fire distribution of a selection of fatality rows.
    Single-death fires are flagged in the source. Multiple-death fires are rebuilt by grouping the
    flagged rows on the fields that describe the incident."""
    d = Counter()
    d[1] = sum(1 for r in sel if r["FATALITY_TOTAL"] != MULTI)
    for size in Counter(tuple(r[k] for k in KEY) for r in sel if r["FATALITY_TOTAL"] == MULTI).values():
        d[size] += 1
    return dict(sorted(d.items()))


CUT = {
    "all fires": lambda r: True,
    "buildings": lambda r: r["INCIDENT_LOCATION_TYPE"] in ("Dwellings", "Other Buildings"),
    "dwellings": lambda r: r["INCIDENT_LOCATION_TYPE"] == "Dwellings",
    "dwellings, accidental": lambda r: r["INCIDENT_LOCATION_TYPE"] == "Dwellings" and r["ACCIDENTAL_OR_DELIBERATE"] == "Recorded as Accidental",
    "dwellings, deliberate": lambda r: r["INCIDENT_LOCATION_TYPE"] == "Dwellings" and r["ACCIDENTAL_OR_DELIBERATE"] == "Recorded as Deliberate",
    "other buildings": lambda r: r["INCIDENT_LOCATION_TYPE"] == "Other Buildings",
    "road vehicles": lambda r: r["INCIDENT_LOCATION_TYPE"] == "Road Vehicles",
    "outdoors": lambda r: r["INCIDENT_LOCATION_TYPE"] == "Other Outdoors",
}


def eng(cut, fys):
    return eng_dist([r for r in HO if r["FINANCIAL_YEAR"] in fys and CUT[cut](r)])


out = dict(method="rates by number of deaths on matched exposure; nothing pooled across datasets")

# ---- 1. validation of the reconstruction against the official counts of fatal fires -------------
val = dict(cells=0, cells_agree=0, by_location={})
for loc in ("Dwellings", "Other Buildings", "Road Vehicles", "Other Outdoors"):
    off_m = rec_m = dev = 0
    for a in OFF:
        if a["location"] != loc:
            continue
        sel = [r for r in HO if r["FINANCIAL_YEAR"] == a["financial_year"] and r["INCIDENT_LOCATION_TYPE"] == loc]
        d = eng_dist(sel)
        m = sum(c for k, c in d.items() if k >= 2)
        ok = (m == int(a["events_multiple"]) and d.get(1, 0) == int(a["events_single"]) and len(sel) == int(a["deaths"]))
        val["cells"] += 1; val["cells_agree"] += int(ok)
        off_m += int(a["events_multiple"]); rec_m += m; dev += abs(m - int(a["events_multiple"]))
    val["by_location"][loc] = dict(official_multiple_fatality_fires=off_m, reconstructed=rec_m, sum_abs_difference=dev)
out["england_reconstruction_check"] = val

# ---- 2. Singapore and England by number of deaths, three windows --------------------------------
def sg_block(wkey, basis):
    w = SG_W[wkey]; c = counts(dist_of(sg_sel(w["rows"], basis))); T = len(w["years"])
    dw = mean_of("Singapore", "dwellings", w["years"]) * T / 1e6
    po = mean_of("Singapore", "population_total", w["years"]) * T / 1e6
    hh = mean_of("Singapore", "households", w["years"]) * T / 1e6
    hh2 = (mean_of("Singapore", "households", w["years"]) - mean_of("Singapore", "households_one_person", w["years"])) * T / 1e6
    c.update(window=w["label"], years=T, exposure=dict(dwelling_years_M=round(dw, 3), person_years_M=round(po, 2),
             household_years_M=round(hh, 3), multi_person_household_years_M=round(hh2, 3)))
    return c


def eng_block(wkey, cut="dwellings"):
    fys = ENG_W[wkey]; c = counts(eng(cut, fys)); T = len(fys)
    dw = sum(eng_fy_dwellings(f) for f in fys) / 1e6
    po = sum(eng_fy_population(f) for f in fys) / 1e6
    hh = sum(eng_households(int(f[:4]) + 0.5) for f in fys) / 1e6
    hh2 = sum(eng_households(int(f[:4]) + 0.5) - eng_one_person(int(f[:4]) + 0.5) for f in fys) / 1e6
    c.update(window="%s to %s" % (fys[0], fys[-1]), years=T, exposure=dict(dwelling_years_M=round(dw, 3), person_years_M=round(po, 2),
             household_years_M=round(hh, 3), multi_person_household_years_M=round(hh2, 3)))
    return c


def compare(s, e, classes=("single", "multi", "ge3", "ge4", "fatal_fires"), denoms=("dwelling_years_M", "person_years_M")):
    res = {}
    for dn in denoms:
        res[dn] = {}
        for cl in classes:
            res[dn][cl] = dict(singapore=rate(s[cl], s["exposure"][dn]), england=rate(e[cl], e["exposure"][dn]),
                               ratio=rate_ratio(s[cl], s["exposure"][dn], e[cl], e["exposure"][dn]))
    # identity: odds ratio of the multi-fatality share = (rate ratio at two or more) / (rate ratio at one)
    orr = (s["multi"] / s["single"]) / (e["multi"] / e["single"])
    res["identity"] = dict(share_singapore=round(s["multi"] / s["fatal_fires"], 4), share_england=round(e["multi"] / e["fatal_fires"], 4),
                           share_ci_singapore=[round(x, 3) for x in cp(s["multi"], s["fatal_fires"])],
                           share_ci_england=[round(x, 3) for x in cp(e["multi"], e["fatal_fires"])],
                           odds_ratio=round(orr, 2),
                           odds_ratio_ci95=[round(float(x), 2) for x in exact_odds_ratio(np.array([[s["multi"], s["single"]], [e["multi"], e["single"]]])).confidence_interval(0.95)],
                           ratio_multi_over_ratio_single_dwellings=round(res["dwelling_years_M"]["multi"]["ratio"]["ratio"] / res["dwelling_years_M"]["single"]["ratio"]["ratio"], 2),
                           fisher_p=float("%.1e" % fisher_exact([[s["multi"], s["single"]], [e["multi"], e["single"]]])[1]))
    return res


cmpx = {}
for wkey in ("primary", "overlap", "review"):
    e = eng_block(wkey)
    cmpx[wkey] = dict(england_dwellings=e)
    for basis in ("residential", "dwelling_stock", "household_units"):
        s = sg_block(wkey, basis)
        cmpx[wkey]["singapore_" + basis] = s
        cmpx[wkey]["compare_" + basis] = compare(s, e)
out["singapore_vs_england"] = cmpx

# further denominators (primary window). The primary Singapore basis is fires in premises counted in the
# dwelling stock, so that numerator and denominator match; the wider residential basis is kept beside it.
s, e = cmpx["primary"]["singapore_dwelling_stock"], cmpx["primary"]["england_dwellings"]
s_res = cmpx["primary"]["singapore_residential"]
out["primary_basis"] = "Singapore: fatal fires in premises counted in the dwelling stock, 2012-2025; England: dwelling fires, 2012/13 to 2025/26"
out["denominator_sensitivity"] = {
    dn: {cl: rate_ratio(s[cl], s["exposure"][dn], e[cl], e["exposure"][dn]) for cl in ("single", "multi")}
    for dn in ("dwelling_years_M", "person_years_M", "household_years_M")}
out["denominator_sensitivity_residential_basis"] = {
    dn: {cl: rate_ratio(s_res[cl], s_res["exposure"][dn], e[cl], e["exposure"][dn]) for cl in ("single", "multi")}
    for dn in ("dwelling_years_M", "person_years_M", "household_years_M")}
out["denominator_sensitivity"]["multi_person_household_years_M"] = dict(
    multi=rate_ratio(s["multi"], s["exposure"]["multi_person_household_years_M"], e["multi"], e["exposure"]["multi_person_household_years_M"]),
    one_person_share_of_households=dict(
        singapore=round(mean_of("Singapore", "households_one_person", SG_W["primary"]["years"]) / mean_of("Singapore", "households", SG_W["primary"]["years"]), 3),
        england=round(np.mean([eng_one_person(y + 0.5) / eng_households(y + 0.5) for y in range(2012, 2026)]), 3)))
res_pop = mean_of("Singapore", "population_resident", SG_W["primary"]["years"]) * 14 / 1e6
out["denominator_sensitivity"]["resident_person_years_M"] = {cl: rate_ratio(s[cl], res_pop, e[cl], e["exposure"]["person_years_M"]) for cl in ("single", "multi")}

# sensitivities of the numerators (primary window, per dwelling)
sens = {}
dwS, dwE = s["exposure"]["dwelling_years_M"], e["exposure"]["dwelling_years_M"]
rows_p = sg_sel(SG_W["primary"]["rows"], "dwelling_stock")
no_hom = counts(dist_of([r for r in rows_p if r["incident_id"] != "SGF-2012-02"]))
sens["without_homicide_by_fire_SGF_2012_02"] = dict(single=rate_ratio(no_hom["single"], dwS, e["single"], dwE))
eacc = eng_block("primary", "dwellings, accidental")
sens["england_accidental_only"] = dict(single=rate_ratio(s["single"], dwS, eacc["single"], dwE), multi=rate_ratio(s["multi"], dwS, eacc["multi"], dwE),
                                       england_single=eacc["single"], england_multi=eacc["multi"])
OFFICIAL_DEATHS = dict(portal_series=57, with_reported_2024_revision=59)
reg_all = [r for r in rd("fatal_incident_register.csv") if 2012 <= int(r["year"]) <= 2025]
all_recorded_deaths = sum(int(r["fatalities"]) for r in reg_all)
sens["capture_bounds"] = dict(
    note="upper bounds on the Singapore fatal dwelling-fire rate: every death in the official series (any type of fire), "
         "or every death in the register including boundary events, taken as a separate single-death dwelling fire",
    official_deaths_2012_2025=OFFICIAL_DEATHS, register_deaths_including_boundary_events=all_recorded_deaths,
    bound_official=rate_ratio(OFFICIAL_DEATHS["with_reported_2024_revision"], dwS, e["fatal_fires"], dwE),
    bound_all_recorded=rate_ratio(all_recorded_deaths, dwS, e["fatal_fires"], dwE))
sens["missed_single_death_fires_needed_for_equal_single_rate"] = int(round(e["single"] / dwE * dwS - s["single"]))
out["numerator_sensitivity"] = sens

# ---- 3. how far the counts move under other definitions of a fire death (England, dwellings) -----
def definitional(fys):
    dw = [r for r in HO if r["FINANCIAL_YEAR"] in fys and r["INCIDENT_LOCATION_TYPE"] == "Dwellings"]

    def after(keep):
        single = sum(1 for r in dw if r["FATALITY_TOTAL"] != MULTI and keep(r))
        g = Counter(tuple(r[k] for k in KEY) for r in dw if r["FATALITY_TOTAL"] == MULTI and keep(r))
        return single, sum(1 for v in g.values() if v >= 2), sum(1 for v in g.values() if v == 1)

    b1, b2, _ = after(lambda r: True)
    F = {
        "deaths after leaving the scene alive excluded": lambda r: r["FATALITY_CIRCUMSTANCES_DESCRIPTION"] != "Alive on leaving scene, but died later",
        "deaths with cause recorded as other or not known excluded": lambda r: r["FATALITY_CAUSE"] not in ("Other", "Not known"),
        "deliberate fires excluded": lambda r: r["ACCIDENTAL_OR_DELIBERATE"] != "Recorded as Deliberate",
        "both death exclusions": lambda r: r["FATALITY_CIRCUMSTANCES_DESCRIPTION"] != "Alive on leaving scene, but died later" and r["FATALITY_CAUSE"] not in ("Other", "Not known"),
        "both death exclusions and deliberate fires excluded": lambda r: (r["FATALITY_CIRCUMSTANCES_DESCRIPTION"] != "Alive on leaving scene, but died later"
                                                                          and r["FATALITY_CAUSE"] not in ("Other", "Not known") and r["ACCIDENTAL_OR_DELIBERATE"] != "Recorded as Deliberate"),
    }
    res = dict(base=dict(single=b1, multi=b2))
    for name, keep in F.items():
        s1, m, dem = after(keep)
        # a multiple-death fire left with one qualifying death is counted as a single-death fire
        res[name] = dict(single=s1 + dem, multi=m, single_change_pct=round(100 * ((s1 + dem) / b1 - 1), 1), multi_change_pct=round(100 * (m / b2 - 1), 1))
    return res


out["definition_sensitivity_england_dwellings"] = dict(primary=definitional(ENG_W["primary"]), review=definitional(ENG_W["review"]))
strict = out["definition_sensitivity_england_dwellings"]["primary"]["both death exclusions and deliberate fires excluded"]
out["numerator_sensitivity"]["england_strictest_definition"] = dict(
    note="England restricted to accidental dwelling fires and to deaths at the scene with a recorded fire cause; Singapore unrestricted",
    england_single=strict["single"], england_multi=strict["multi"],
    single=rate_ratio(s["single"], dwS, strict["single"], dwE), multi=rate_ratio(s["multi"], dwS, strict["multi"], dwE))

# ---- 4. two kinds of fatal dwelling fire: victims of single- and multiple-death fires (England) ---
dw_all = [r for r in HO if r["INCIDENT_LOCATION_TYPE"] == "Dwellings"]
SV = [r for r in dw_all if r["FATALITY_TOTAL"] != MULTI]
MV = [r for r in dw_all if r["FATALITY_TOTAL"] == MULTI and not is_grenfell(r)]
NIGHT = {"00-01", "01-02", "02-03", "03-04", "04-05", "05-06"}


def age_lo(a):
    return None if a == "Age Not Known" else (90 if a == "Aged 90 or over" else int(a.split()[1]))


def share(sel, test):
    k = sum(1 for r in sel if test(r))
    return dict(k=k, n=len(sel), share=round(k / len(sel), 3), ci95=[round(x, 3) for x in cp(k, len(sel))])


TESTS = {
    "aged under 15": lambda r: age_lo(r["VICTIM_AGE"]) is not None and age_lo(r["VICTIM_AGE"]) < 15,
    "aged 15 to 39": lambda r: age_lo(r["VICTIM_AGE"]) is not None and 15 <= age_lo(r["VICTIM_AGE"]) < 40,
    "aged 40 to 64": lambda r: age_lo(r["VICTIM_AGE"]) is not None and 40 <= age_lo(r["VICTIM_AGE"]) < 65,
    "aged 65 or over": lambda r: age_lo(r["VICTIM_AGE"]) is not None and age_lo(r["VICTIM_AGE"]) >= 65,
    "in the room of origin when the fire started": lambda r: r["VICTIM_LOCATION_START"] == "Room, cabin or compartment of fire origin",
    "elsewhere in the building when the fire started": lambda r: r["VICTIM_LOCATION_START"].startswith("Different floor"),
    "found in the room of origin": lambda r: r["VICTIM_LOCATION_FOUND"] == "Room, cabin or compartment of fire origin",
    "found elsewhere in the building": lambda r: r["VICTIM_LOCATION_FOUND"].startswith("Different floor"),
    "fire reported between 00:00 and 06:00": lambda r: r["HOUR_BANDS"] in NIGHT,
    "ignition by smokers' materials": lambda r: r["SOURCE_OF_IGNITION"] == "Smokers' materials",
    "fire recorded as deliberate": lambda r: r["ACCIDENTAL_OR_DELIBERATE"] == "Recorded as Deliberate",
    "death by burns alone": lambda r: r["FATALITY_CAUSE"] == "Burns",
    "overcome by gas or smoke": lambda r: r["FATALITY_CAUSE"] == "Overcome by gas or smoke",
    "left the scene alive and died later": lambda r: r["FATALITY_CIRCUMSTANCES_DESCRIPTION"] == "Alive on leaving scene, but died later",
    "no alarm present": lambda r: r["ALARM_SYSTEM"] == "Alarm Absent",
}
out["two_processes_england_dwellings"] = dict(
    window="2010/11 to 2025/26", note="one row per death; the 71 Grenfell Tower records are set aside",
    victims_single=len(SV), victims_multiple=len(MV),
    shares={k: dict(single=share(SV, t), multiple=share(MV, t)) for k, t in TESTS.items()})
hours = sorted({r["HOUR_BANDS"] for r in HO})
out["two_processes_england_dwellings"]["by_hour"] = dict(
    hours=hours, single=[round(sum(1 for r in SV if r["HOUR_BANDS"] == h) / len(SV), 4) for h in hours],
    multiple=[round(sum(1 for r in MV if r["HOUR_BANDS"] == h) / len(MV), 4) for h in hours])
BANDS5 = ["0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
          "65-69", "70-74", "75-79", "80-84", "85+"]


def band5(a):
    lo = age_lo(a)
    return None if lo is None else ("85+" if lo >= 85 else "%d-%d" % (lo, lo + 4))


out["two_processes_england_dwellings"]["by_age"] = dict(
    bands=BANDS5, single=[sum(1 for r in SV if band5(r["VICTIM_AGE"]) == b) for b in BANDS5],
    multiple=[sum(1 for r in MV if band5(r["VICTIM_AGE"]) == b) for b in BANDS5])

# ---- 5. age standardisation of the single-death rate (indirect, England as the standard) ---------
AGE = defaultdict(dict)
for r in rd("population_by_age.csv"):
    AGE[(r["jurisdiction"], r["basis"])][(int(r["year"]), r["age_group"])] = float(r["value"])


def standardise(multi_flag, eng_fys, sg_years, sg_basis):
    sel = [r for r in HO if r["INCIDENT_LOCATION_TYPE"] == "Dwellings" and r["FINANCIAL_YEAR"] in eng_fys
           and (r["FATALITY_TOTAL"] == MULTI) == multi_flag and not is_grenfell(r)]
    known = [r for r in sel if band5(r["VICTIM_AGE"])]
    scale = len(sel) / len(known)                              # deaths of unknown age spread in proportion
    exp_total = 0.0; table = []
    for b in BANDS5:
        d = sum(1 for r in known if band5(r["VICTIM_AGE"]) == b) * scale
        py_e = sum(AGE[("England", "mid-year estimates")][(int(f[:4]), b)] for f in eng_fys)
        if sg_basis == "total":
            py_s = 0.0
            for y in sg_years:                                 # UN age shares applied to the official total
                tot = sum(AGE[("Singapore", "total population (UN WPP 2024)")][(y, bb)] for bb in BANDS5)
                py_s += AGE[("Singapore", "total population (UN WPP 2024)")][(y, b)] / tot * EXPO[("Singapore", "population_total")][y]
        else:
            py_s = sum(AGE[("Singapore", "resident population")][(y, b)] for y in sg_years)
        exp_total += d / py_e * py_s
        table.append(dict(age=b, england_deaths=round(d, 1), england_rate_per_M=round(d / py_e * 1e6, 3), singapore_person_years=round(py_s)))
    return exp_total, table, len(sel)


std = {}
obs_single = s["single"]
obs_multi_deaths = s["deaths_in_multi"]
for basis in ("total", "resident"):
    e1, tab1, n1 = standardise(False, ENG_W["primary"], SG_W["primary"]["years"], basis)
    e2, tab2, n2 = standardise(True, ENG_W["primary"], SG_W["primary"]["years"], basis)
    lo1, hi1 = pois_ci(obs_single); lo2, hi2 = pois_ci(obs_multi_deaths)
    std[basis] = dict(
        single_death_fires=dict(observed=obs_single, expected_at_english_age_specific_rates=round(e1, 1),
                                ratio=round(obs_single / e1, 4), ci95=[round(lo1 / e1, 4), round(hi1 / e1, 4)], english_deaths_used=n1),
        deaths_in_multiple_death_fires=dict(observed=obs_multi_deaths, expected_at_english_age_specific_rates=round(e2, 1),
                                            ratio=round(obs_multi_deaths / e2, 3),
                                            note="deaths in one fire are not independent, so no interval is given", english_deaths_used=n2))
    if basis == "total":
        std["age_table_single"] = tab1
py65 = lambda jur, basis, yrs: sum(AGE[(jur, basis)][(y, b)] for y in yrs for b in BANDS5[13:]) / sum(AGE[(jur, basis)][(y, b)] for y in yrs for b in BANDS5)
std["share_aged_65_or_over"] = dict(
    england=round(py65("England", "mid-year estimates", [int(f[:4]) for f in ENG_W["primary"]]), 3),
    singapore_total=round(py65("Singapore", "total population (UN WPP 2024)", SG_W["primary"]["years"]), 3),
    singapore_resident=round(py65("Singapore", "resident population", SG_W["primary"]["years"]), 3))
std["crude_ratio_per_person_single"] = out["singapore_vs_england"]["primary"]["compare_residential"]["person_years_M"]["single"]["ratio"]
out["age_standardisation"] = std

# ---- 5b. standardisation for household composition (indirect, England as the standard) -----------
# The Home Office review tabulates fatal dwelling fires by type of household for 2010/11 to 2018/19.
# A household of one person can have only a single-death fire (a visitor apart), so the table splits
# the English single-death fires between one-person households and households of two or more.
HH = rd("england_fatal_dwelling_fires_by_household.csv", "utf-8")
hh_total = sum(int(r["fatal_dwelling_fires"]) for r in HH)
hh_lone = sum(int(r["fatal_dwelling_fires"]) for r in HH if r["persons_in_household"] == "one")
hh_unstated = sum(int(r["fatal_dwelling_fires"]) for r in HH if r["persons_in_household"] == "not stated")
rev = cmpx["review"]["england_dwellings"]                       # record-level counts for the same nine years
fy9 = ENG_W["review"]
hy_lone_e = sum(eng_one_person(int(f[:4]) + 0.5) for f in fy9) / 1e6
hy_multi_e = sum(eng_households(int(f[:4]) + 0.5) - eng_one_person(int(f[:4]) + 0.5) for f in fy9) / 1e6
yrs_sg = SG_W["primary"]["years"]
hy_lone_s = sum(EXPO[("Singapore", "households_one_person")][y] for y in yrs_sg) / 1e6
hy_multi_s = sum(EXPO[("Singapore", "households")][y] - EXPO[("Singapore", "households_one_person")][y] for y in yrs_sg) / 1e6


def hh_standard(lone_fires, s=s):
    single_in_multi = hh_total - lone_fires - rev["multi"]      # single-death fires in households of two or more
    r_lone, r_m1, r_m2 = lone_fires / hy_lone_e, single_in_multi / hy_multi_e, rev["multi"] / hy_multi_e
    exp_single = r_lone * hy_lone_s + r_m1 * hy_multi_s
    exp_multi = r_m2 * hy_multi_s
    lo1, hi1 = pois_ci(s["single"]); lo2, hi2 = pois_ci(s["multi"])
    return dict(english_rates_per_M_household_years=dict(one_person=round(r_lone, 2), two_or_more_single_death=round(r_m1, 2), two_or_more_multiple_death=round(r_m2, 3)),
                rate_ratio_one_person_to_two_or_more=round(r_lone / (r_m1 + r_m2), 2),
                expected_single=round(exp_single, 1), observed_single=s["single"], ratio_single=round(s["single"] / exp_single, 4),
                ratio_single_ci95=[round(lo1 / exp_single, 4), round(hi1 / exp_single, 4)],
                expected_multi=round(exp_multi, 1), observed_multi=s["multi"], ratio_multi=round(s["multi"] / exp_multi, 4),
                ratio_multi_ci95=[round(lo2 / exp_multi, 4), round(hi2 / exp_multi, 4)],
                expected_share_of_fires_with_two_or_more=round(exp_multi / (exp_multi + exp_single), 3),
                english_share_given_household_of_two_or_more=round(rev["multi"] / (single_in_multi + rev["multi"]), 3))


known = hh_total - hh_unstated
out["household_standardisation"] = dict(
    source="Home Office (2022) in-depth review, Table 38; fatal dwelling fires 2010/11 to 2018/19",
    english_fatal_dwelling_fires=hh_total, in_one_person_households=hh_lone, share_one_person=round(hh_lone / hh_total, 3),
    type_other_or_not_known=hh_unstated, english_one_person_share_of_households=round(hy_lone_e / (hy_lone_e + hy_multi_e), 3),
    singapore_one_person_share_of_resident_households=round(hy_lone_s / (hy_lone_s + hy_multi_s), 3),
    unstated_counted_as_two_or_more=hh_standard(hh_lone),
    unstated_spread_in_proportion=hh_standard(hh_lone + hh_unstated * hh_lone / known),
    residential_basis_unstated_counted_as_two_or_more=hh_standard(hh_lone, s_res),
    crude_ratio_per_household=out["denominator_sensitivity"]["household_years_M"],
    note="English rates are for 2010/11 to 2018/19 and are applied to Singapore resident households of 2012 to 2025")

# ---- 6. the shape of the curve above one death -------------------------------------------------
def shape(dist, label, within=None):
    t = {k: v for k, v in dist.items() if k >= 2 and (within is None or k <= within)}
    m = sum(t.values()); two = t.get(2, 0)
    b, ll, ci = fit_tail(exact(t))
    return dict(label=label, multi=m, two=two, ge3=m - two, ge4=sum(v for k, v in t.items() if k >= 4), ge5=sum(v for k, v in t.items() if k >= 5),
                share_exactly_two=round(two / m, 3), ci95=[round(x, 3) for x in cp(two, m)], b2=round(b, 2), b2_ci95=[round(x, 2) for x in ci],
                loglik=round(ll, 3), distribution={str(k): v for k, v in sorted(t.items())})


def lrt_common(dists):
    free = sum(fit_tail(exact(d))[1] for d in dists)
    pooled = Counter()
    for d in dists:
        for k, v in d.items():
            if k >= 2:
                pooled[k] += v
    bc, llc, ci = fit_tail(exact(pooled))
    stat = 2 * (free - llc)
    return dict(b2_common=round(bc, 2), ci95=[round(x, 2) for x in ci], lrt=round(stat, 2), df=len(dists) - 1, p=round(float(chi2.sf(stat, len(dists) - 1)), 3))


ordinary = lambda d: {k: v for k, v in d.items() if k <= 10}        # sets aside Grenfell Tower (71) and the Shoreham crash (11)
sh = {}
for cut in CUT:
    sh["England, " + cut] = shape(ordinary(eng(cut, ENG_W["all"])), "England 2010/11 to 2025/26, " + cut)
sg_res_d = dist_of(sg_sel(SG_W["primary"]["rows"], "residential")); sg_bld_d = dist_of(SG_W["primary"]["rows"])
sh["Singapore, residential"] = shape(sg_res_d, "Singapore 2012-2025, residential")
sh["Singapore, buildings"] = shape(sg_bld_d, "Singapore 2012-2025, building fires")
sh["Singapore, dwelling stock"] = shape(dist_of(sg_sel(SG_W["primary"]["rows"], "dwelling_stock")), "Singapore 2012-2025, premises in the dwelling stock")
sh["United States, residential buildings"] = dict(
    label="United States 2021-2023, residential buildings, shares among reported fires", share_exactly_two=0.768,
    note="published percentage (76.8% two deaths, 14.3% three, 4.7% four, 4.2% five or more); counts not published")
BASES = defaultdict(list)
for r in rd("revision_bases.csv", "utf-8"):
    BASES[r["dataset"]].append((int(r["bin_lo"]), int(r["bin_hi"]), int(r["events"])))
for key, lab in (("Australia_preventable_residential", "Australia FY2003/04 to 2016/17, residential, intent excluded"),
                 ("NewZealand_unintentional_domestic", "New Zealand 2007-2014, domestic, intent excluded")):
    bins = [x for x in BASES[key] if x[0] >= 2]
    m = sum(c for _, _, c in bins); two = sum(c for lo, hi, c in bins if lo == 2 and hi == 2)
    b, ll, ci = fit_tail(bins)
    sh[key] = dict(label=lab, multi=m, two=two, ge3=m - two, share_exactly_two=round(two / m, 3), ci95=[round(x, 3) for x in cp(two, m)],
                   b2=round(b, 2), b2_ci95=[round(x, 2) for x in ci], loglik=round(ll, 3))
out["shape_above_one_death"] = sh
e_dw = ordinary(eng("dwellings", ENG_W["all"])); e_acc = ordinary(eng("dwellings, accidental", ENG_W["all"])); e_del = ordinary(eng("dwellings, deliberate", ENG_W["all"]))
tests = {}
tests["singapore_residential_vs_england_dwellings"] = dict(
    fisher_p_two_vs_three_or_more=round(float(fisher_exact([[sh["Singapore, residential"]["two"], sh["Singapore, residential"]["ge3"]],
                                                            [sh["England, dwellings"]["two"], sh["England, dwellings"]["ge3"]]])[1]), 3),
    common_exponent=lrt_common([sg_res_d, e_dw]))
tests["england_locations"] = lrt_common([ordinary(eng(c, ENG_W["all"])) for c in ("dwellings", "other buildings", "road vehicles", "outdoors")])
tests["england_dwellings_accidental_vs_deliberate"] = dict(
    fisher_p_two_vs_three_or_more=round(float(fisher_exact([[sh["England, dwellings, accidental"]["two"], sh["England, dwellings, accidental"]["ge3"]],
                                                            [sh["England, dwellings, deliberate"]["two"], sh["England, dwellings, deliberate"]["ge3"]]])[1]), 3),
    common_exponent=lrt_common([e_acc, e_del]))
for key in ("Australia_preventable_residential", "NewZealand_unintentional_domestic"):
    a = sh[key]
    tests[key + "_vs_england_dwellings_accidental"] = dict(
        fisher_p_two_vs_three_or_more=round(float(fisher_exact([[a["two"], a["ge3"]], [sh["England, dwellings, accidental"]["two"], sh["England, dwellings, accidental"]["ge3"]]])[1]), 3))
out["shape_tests"] = tests
# continuation probabilities h(n) = P(N >= n+1 | N >= n), English dwelling fires, all years, ordinary fires
cont = {}
for n in (1, 2, 3, 4, 5):
    full = ordinary(eng("dwellings", ENG_W["all"]))
    at = sum(v for k, v in full.items() if k >= n); nxt = sum(v for k, v in full.items() if k >= n + 1)
    cont[str(n)] = dict(at_least_n=at, at_least_n_plus_1=nxt, h=round(nxt / at, 3), ci95=[round(x, 3) for x in cp(nxt, at)])
out["continuation_probabilities_england_dwellings"] = cont

# ---- 7. one-parameter law against the two-part law ----------------------------------------------
def one_part_ll(dist, cap=10 ** 4):
    i = np.arange(1, cap + 1, dtype=float)

    def nll(b):
        w = i ** (-b); z = w.sum()
        return -sum(c * np.log(w[k - 1] / z) for k, c in dist.items())

    r = minimize_scalar(nll, bounds=(1.05, 9.0), method="bounded")
    return float(r.x), float(-r.fun)


def two_part(dist):
    n = sum(dist.values()); m = n - dist.get(1, 0)
    b1, ll1 = one_part_ll(dist)
    b2, ll_t, ci = fit_tail(exact(dist))
    ll2 = (n - m) * np.log((n - m) / n) + m * np.log(m / n) + ll_t
    stat = 2 * (ll2 - ll1)
    q_implied = 1 - 1 / float(np.sum(np.arange(1, 10 ** 4 + 1, dtype=float) ** (-b1)))
    return dict(fatal_fires=n, multi=m, b_one_parameter=round(b1, 2), share_implied_by_one_parameter=round(q_implied, 4), share_observed=round(m / n, 4),
                b2=round(b2, 2), b2_ci95=[round(x, 2) for x in ci], lrt=round(stat, 2),
                p=round(float(chi2.sf(stat, 1)), 5))       # five places: a value such as 0.17498 must not be stored as 0.175 and then read as 0.18


out["two_part_law"] = {
    "Singapore, building fires": two_part(sg_bld_d),
    "Singapore, residential": two_part(sg_res_d),
    "England, dwellings 2010/11 to 2025/26": two_part(ordinary(eng("dwellings", ENG_W["all"]))),
    "England, dwellings 2010/11 to 2018/19": two_part(ordinary(eng("dwellings", ENG_W["review"]))),
    "England, all fires 2010/11 to 2025/26": two_part(ordinary(eng("all fires", ENG_W["all"]))),
}

# ---- 8. what the English shape implies for Singapore, and for Grenfell Tower ---------------------
b2_e = sh["England, dwellings"]["b2"]; b2_ci = sh["England, dwellings"]["b2_ci95"]
lam2_sg_res = s["multi"] / 14.0
lam2_sg_bld = counts(sg_bld_d)["multi"] / 14.0
scen = {}
for cap, lab in ((8, "cap of 8"), (17, "cap of 17"), (10 ** 4, "no cap")):
    scen[lab] = {}
    for n0 in (5, 10, 30):
        sv = tail_survival(b2_e, n0, cap=cap)
        scen[lab]["return_period_%d_or_more_yr" % n0] = (round(1 / (lam2_sg_bld * sv), 1) if sv > 0 else None)
    scen[lab]["at_interval_edges_5_or_more_yr"] = [round(1 / (lam2_sg_bld * tail_survival(b, 5, cap=cap)), 1) for b in b2_ci]
out["singapore_tail_with_english_shape"] = dict(
    note="Singapore's own rate of multiple-death building fires with the shape of English dwelling fires above one death",
    multiple_death_fires_per_year_buildings=round(lam2_sg_bld, 3), multiple_death_fires_per_year_residential=round(lam2_sg_res, 3),
    english_b2=b2_e, english_b2_ci95=b2_ci, scenarios=scen,
    probability_no_fire_of_5_or_more_in_14_years=round(float(np.exp(-14 * lam2_sg_bld * tail_survival(b2_e, 5))), 2))
# monitoring: Bayes factor for a cap at the observed maximum (four deaths) against no cap, after y quiet years,
# under the three-part law with the English shape; the shape is integrated over a lognormal fitted to its interval
zq = np.linspace(-6, 6, 2001); wq = np.exp(-zq ** 2 / 2); wq /= wq.sum()
sd_b2 = (np.log(b2_ci[1]) - np.log(b2_ci[0])) / (2 * 1.96)
s5 = np.array([tail_survival(float(b), 5) for b in np.exp(np.log(b2_e) + sd_b2 * zq)])
out["singapore_tail_with_english_shape"]["monitoring_bayes_factor_cap4_vs_no_cap"] = {
    "after_%dyr" % y: round(1 / float(np.sum(wq * (1 - s5) ** (lam2_sg_bld * y))), 2) for y in (10, 20, 30, 50)}
own_b2, _, own_ci = fit_tail(exact(sg_bld_d))
out["singapore_tail_with_english_shape"]["singapore_own_shape"] = dict(
    b2=round(own_b2, 2), ci95=[round(x, 2) for x in own_ci],
    return_period_5_or_more_yr=round(1 / (lam2_sg_bld * tail_survival(own_b2, 5)), 1),
    return_period_10_or_more_yr=round(1 / (lam2_sg_bld * tail_survival(own_b2, 10)), 1))
gren = {}
for wkey in ("review", "all"):
    d = ordinary(eng("dwellings", ENG_W[wkey]))
    b2g, _, cig = fit_tail(exact(d))
    lam2 = sum(v for k, v in d.items() if k >= 2) / len(ENG_W[wkey])
    gren[wkey] = dict(window="%s to %s" % (ENG_W[wkey][0], ENG_W[wkey][-1]), ordinary_multiple_death_fires=sum(v for k, v in d.items() if k >= 2),
                      per_year=round(lam2, 2), b2=round(b2g, 2), b2_ci95=[round(x, 2) for x in cig],
                      P_70_or_more_given_two_or_more=float("%.2e" % tail_survival(b2g, 70)),
                      return_period_yr=int(round(1 / (lam2 * tail_survival(b2g, 70)), -2)),
                      return_period_yr_flat_edge=int(round(1 / (lam2 * tail_survival(cig[0], 70)), -2)),
                      P_at_least_one_in_window=float("%.2e" % (1 - np.exp(-len(ENG_W[wkey]) * lam2 * tail_survival(b2g, 70)))),
                      P_at_least_one_in_window_flat_edge=float("%.2e" % (1 - np.exp(-len(ENG_W[wkey]) * lam2 * tail_survival(cig[0], 70)))))
out["grenfell_under_shape_above_one_death"] = gren

# ---- 9. other published datasets, each against England re-cut to its scope ----------------------
def expo_mean(jur, measure, years):
    sset = EXPO.get((jur, measure))
    if not sset:
        return None
    ys = [y for y in years if y in sset]
    return float(np.mean([sset[y] for y in ys])) if ys else None


def eng_pair(cut, fys):
    c = counts(eng(cut, fys))
    c["exposure"] = dict(dwelling_years_M=sum(eng_fy_dwellings(f) for f in fys) / 1e6, person_years_M=sum(eng_fy_population(f) for f in fys) / 1e6,
                         household_years_M=sum(eng_households(int(f[:4]) + 0.5) for f in fys) / 1e6)
    c["window"] = "%s to %s" % (fys[0], fys[-1])
    return c


pairs = {}
PAIR_DEF = [
    ("Australia", "Australia_preventable_residential", 14, list(range(2004, 2018)), "dwellings, accidental", FYS[0:7]),
    ("New Zealand", "NewZealand_unintentional_domestic", 8, list(range(2007, 2015)), "dwellings, accidental", FYS[0:5]),
]
for jur, key, T, years, cut, fys in PAIR_DEF:
    bins = BASES[key]
    single = sum(c for lo, hi, c in bins if lo == 1); multi = sum(c for lo, hi, c in bins if lo >= 2)
    hh = expo_mean(jur, "households", years); po = expo_mean(jur, "population_total", years)
    ep = eng_pair(cut, fys)
    p = dict(scope_of_source=key, england_recut=cut, england_window=ep["window"], source_years=T, single=single, multi=multi,
             england_single=ep["single"], england_multi=ep["multi"])
    if hh:
        p["per_household"] = dict(single=rate_ratio(single, hh * T / 1e6, ep["single"], ep["exposure"]["household_years_M"]),
                                  multi=rate_ratio(multi, hh * T / 1e6, ep["multi"], ep["exposure"]["household_years_M"]),
                                  rate_single=rate(single, hh * T / 1e6), rate_multi=rate(multi, hh * T / 1e6),
                                  england_rate_single=rate(ep["single"], ep["exposure"]["household_years_M"]),
                                  england_rate_multi=rate(ep["multi"], ep["exposure"]["household_years_M"]))
    if po:
        p["per_person"] = dict(single=rate_ratio(single, po * T / 1e6, ep["single"], ep["exposure"]["person_years_M"]),
                               multi=rate_ratio(multi, po * T / 1e6, ep["multi"], ep["exposure"]["person_years_M"]))
    pairs[jur] = p
# United States: a scaled national estimate and rounded shares; no interval is attached to it
US = dict(fatal_fires_per_year=1900, share_multi=0.13, share_multi_rounding=[0.125, 0.135])
hh_us = expo_mean("United States", "households", [2018, 2019, 2020]); po_us = expo_mean("United States", "population_total", [2018, 2019, 2020])
ep = eng_pair("dwellings", FYS[8:11])
if hh_us:
    es = ep["single"] / ep["exposure"]["household_years_M"]; em = ep["multi"] / ep["exposure"]["household_years_M"]
    us_s = US["fatal_fires_per_year"] * (1 - US["share_multi"]) / (hh_us / 1e6); us_m = US["fatal_fires_per_year"] * US["share_multi"] / (hh_us / 1e6)
    pairs["United States"] = dict(
        scope_of_source="residential buildings, all causes, 2018-2020; scaled national estimate", england_recut="dwellings", england_window=ep["window"],
        england_single=ep["single"], england_multi=ep["multi"],
        per_household=dict(rate_single=round(us_s, 3), rate_multi=round(us_m, 3), england_rate_single=round(es, 3), england_rate_multi=round(em, 3),
                           ratio_single=round(us_s / es, 2), ratio_multi=round(us_m / em, 2),
                           ratio_multi_rounding=[round(US["fatal_fires_per_year"] * q / (hh_us / 1e6) / em, 2) for q in US["share_multi_rounding"]],
                           england_multi_ratio_ci_only=[round(x / (ep["multi"] / 1.0) if False else x, 2) for x in
                                                        [us_m / (pois_ci(ep["multi"])[1] / ep["exposure"]["household_years_M"]),
                                                         us_m / (pois_ci(ep["multi"])[0] / ep["exposure"]["household_years_M"])]]))
hh_us2 = expo_mean("United States", "households", [2021, 2022, 2023])
ep2 = eng_pair("dwellings", FYS[11:14])
if hh_us2:
    em2 = ep2["multi"] / ep2["exposure"]["household_years_M"]
    pairs["United States"]["second_period_2021_2023"] = dict(
        note="direct national estimate of 200 multiple-death fires a year, published rounded to the nearest 100, stated as 11% of fatal fires",
        england_window=ep2["window"], england_multi=ep2["multi"], england_rate_multi=round(em2, 3),
        rate_multi=round(200 / (hh_us2 / 1e6), 3), rate_multi_rounding=[round(150 / (hh_us2 / 1e6), 3), round(250 / (hh_us2 / 1e6), 3)],
        ratio_multi=round(200 / (hh_us2 / 1e6) / em2, 2), ratio_multi_rounding=[round(150 / (hh_us2 / 1e6) / em2, 2), round(250 / (hh_us2 / 1e6) / em2, 2)])
pairs["United States"]["per_household"]["ratio_multi_range_both_periods"] = [
    round(min(pairs["United States"]["per_household"]["ratio_multi_rounding"][0], pairs["United States"].get("second_period_2021_2023", {}).get("ratio_multi_rounding", [9])[0]), 2),
    round(max(pairs["United States"]["per_household"]["ratio_multi_rounding"][1], pairs["United States"].get("second_period_2021_2023", {}).get("ratio_multi_rounding", [0, 0])[1]), 2)]
out["scope_matched_pairs"] = pairs

# ---- 9b. checks added after independent review ---------------------------------------------------
rev_checks = {}
# (i) p-value on the same (central) basis as the exact interval of the ratio of ratios
for wkey in ("primary", "overlap", "review"):
    for basis in ("residential", "dwelling_stock", "household_units"):
        sb, eb = cmpx[wkey]["singapore_" + basis], cmpx[wkey]["england_dwellings"]
        tab = [[sb["multi"], sb["single"]], [eb["multi"], eb["single"]]]
        pc = min(1.0, 2 * min(fisher_exact(tab, alternative="greater")[1], fisher_exact(tab, alternative="less")[1]))
        cmpx[wkey]["compare_" + basis]["identity"]["p_central"] = float("%.2g" % pc)
# (ii) the screen for premises outside the dwelling stock applied to single-death fires as well
UNVERIFIED_SINGLE = {"SGF-2016-02", "SGF-2015-05"}       # premises recorded only as residential, in shophouse areas
sym = counts(dist_of([r for r in sg_sel(SG_W["primary"]["rows"], "dwelling_stock") if r["incident_id"] not in UNVERIFIED_SINGLE]))
rev_checks["symmetric_screen"] = dict(single=sym["single"], multi=sym["multi"], removed=sorted(UNVERIFIED_SINGLE),
                                      single_ratio=rate_ratio(sym["single"], dwS, e["single"], dwE), multi_ratio=rate_ratio(sym["multi"], dwS, e["multi"], dwE))
# (iii) household type and age together: one person over pensionable age, one person under it, two or more
lone_old = next(int(r["fatal_dwelling_fires"]) for r in HH if r["household_type"].startswith("Lone person over"))
lone_young = next(int(r["fatal_dwelling_fires"]) for r in HH if r["household_type"].startswith("Lone person under"))


def eng_lone_old(year):
    sset = EXPO[("England", "one_person_households_aged_65_plus")]
    return sset[2011] + (sset[2021] - sset[2011]) * (year - 2011) / 10.0


hy_old_e = sum(eng_lone_old(int(f[:4]) + 0.5) for f in fy9) / 1e6
hy_old_s = sum(EXPO[("Singapore", "one_person_households_aged_65_plus")][y] for y in yrs_sg) / 1e6
single_in_multi = hh_total - hh_lone - rev["multi"]
exp_joint = (lone_old / hy_old_e * hy_old_s + lone_young / (hy_lone_e - hy_old_e) * (hy_lone_s - hy_old_s) + single_in_multi / hy_multi_e * hy_multi_s)
lo1, hi1 = pois_ci(s["single"])
strict_p = out["definition_sensitivity_england_dwellings"]["primary"]["both death exclusions and deliberate fires excluded"]
base_p = out["definition_sensitivity_england_dwellings"]["primary"]["base"]
f1 = base_p["single"] / strict_p["single"]; f2 = base_p["multi"] / strict_p["multi"]
hs = out["household_standardisation"]
rev_checks["joint_household_and_age"] = dict(
    english_rates_per_M_household_years=dict(one_person_over_pension_age=round(lone_old / hy_old_e, 2), one_person_under=round(lone_young / (hy_lone_e - hy_old_e), 2),
                                             two_or_more_single_death=round(single_in_multi / hy_multi_e, 2)),
    singapore_one_person_aged_65_plus_share_of_households=round(hy_old_s / (hy_lone_s + hy_multi_s), 3),
    england_one_person_aged_65_plus_share_of_households=round(hy_old_e / (hy_lone_e + hy_multi_e), 3),
    expected_single=round(exp_joint, 1), ratio_single=round(s["single"] / exp_joint, 4), ratio_single_ci95=[round(lo1 / exp_joint, 4), round(hi1 / exp_joint, 4)])
rev_checks["all_allowances_compounded"] = dict(
    note="household and age standardisation multiplied by the effect of the strictest English definition (accidental fires, deaths at the scene, recorded cause) "
         "applied to England alone; an upper bound on the single-death ratio, since the definition is not applied to Singapore",
    definition_factor_single=round(f1, 3), definition_factor_multi=round(f2, 3),
    single_ratio=round(s["single"] / exp_joint * f1, 4), single_ratio_ci95=[round(lo1 / exp_joint * f1, 4), round(hi1 / exp_joint * f1, 4)],
    single_ratio_unstated_spread=round(hs["unstated_spread_in_proportion"]["ratio_single"] * f1, 4),
    single_ratio_unstated_spread_ci95=[round(x * f1, 4) for x in hs["unstated_spread_in_proportion"]["ratio_single_ci95"]],
    multi_ratio=round(hs["unstated_counted_as_two_or_more"]["ratio_multi"] * f2, 4),
    multi_ratio_ci95=[round(x * f2, 4) for x in hs["unstated_counted_as_two_or_more"]["ratio_multi_ci95"]])
# (iv) return period of a fire of five or more deaths with the uncertainty of the rate as well as of the shape
rng = np.random.default_rng(20261120)
k2 = counts(sg_bld_d)["multi"]
lam_draw = rng.gamma(k2 + 0.5, 1 / 14.0, 200000)
b_draw = np.exp(rng.normal(np.log(b2_e), sd_b2, 200000))
grid_b = np.linspace(2.0, 8.0, 601); grid_s5 = np.array([tail_survival(float(b), 5) for b in grid_b]); grid_s10 = np.array([tail_survival(float(b), 10) for b in grid_b])
rp5 = 1 / (lam_draw * np.interp(b_draw, grid_b, grid_s5)); rp10 = 1 / (lam_draw * np.interp(b_draw, grid_b, grid_s10))
rev_checks["return_period_with_rate_uncertainty"] = dict(
    five_or_more_yr=[round(float(np.percentile(rp5, q)), 0) for q in (2.5, 50, 97.5)],
    ten_or_more_yr=[round(float(np.percentile(rp10, q)), 0) for q in (2.5, 50, 97.5)],
    rate_interval_alone_five_or_more_yr=[round(1 / (pois_ci(k2)[1] / 14.0 * tail_survival(b2_e, 5)), 0), round(1 / (pois_ci(k2)[0] / 14.0 * tail_survival(b2_e, 5)), 0)])
# (v) sizes of multiple-death fires: reconstruction against the one published table (all fires, 2010/11 to 2018/19)
pub = {int(lo): c for lo, hi, c in BASES["England_all_fires_within_regime"]}
rec = ordinary(eng("all fires", ENG_W["review"]))
rev_checks["size_check_all_fires_review_window"] = dict(
    published={str(k): v for k, v in sorted(pub.items())}, reconstructed={str(k): v for k, v in sorted(rec.items())},
    b2_published=round(fit_tail(exact(pub))[0], 2), b2_published_ci95=[round(x, 2) for x in fit_tail(exact(pub))[2]],
    b2_reconstructed=round(fit_tail(exact(rec))[0], 2), b2_reconstructed_ci95=[round(x, 2) for x in fit_tail(exact(rec))[2]])
# (vi) Grenfell Tower left in the fit of English dwelling fires, 2010/11 to 2018/19
with_g = eng("dwellings", ENG_W["review"])
bg, _ = one_part_ll(with_g)
ig = np.arange(1, 10 ** 4 + 1, dtype=float); wg = ig ** (-bg); p70 = float(wg[69:].sum() / wg.sum())
rev_checks["grenfell_included_in_fit"] = dict(fatal_dwelling_fires=sum(with_g.values()), b_single_law=round(bg, 2), P_70_or_more_per_fire=float("%.2e" % p70),
                                              P_any_of_the_fires_70_or_more=float("%.2e" % (1 - (1 - p70) ** sum(with_g.values()))))
# (vii) fourteen years without a fire of five or more deaths, under the single law with no cap
from fs_model import fit_b as _fit_b, singapore_bins as _sgb, survival as _surv
_b = _fit_b(_sgb("buildings"), 10 ** 4)[0]
rev_checks["P_no_fire_of_5_or_more_in_14_years_single_law"] = round(float(np.exp(-counts(sg_bld_d)["fatal_fires"] * _surv(_b, 10 ** 4, 5))), 2)
# (viii) the separation of the two rates after the allowances: the multiple-death ratio over the single-death ratio
from scipy.stats import binom as _binom


def hh_expected(lone_fires):
    """expected Singapore fires (one death; two or more) at English rates by household type, unrounded"""
    sim = hh_total - lone_fires - rev["multi"]
    return lone_fires / hy_lone_e * hy_lone_s + sim / hy_multi_e * hy_multi_s, rev["multi"] / hy_multi_e * hy_multi_s


def ror_after(exp_single, exp_multi):
    """(k2 / exp_multi) / (k1 / exp_single); exact interval conditional on k1 + k2, English expectations taken as fixed"""
    k1, k2 = s["single"], s["multi"]
    lo, hi = cp(k2, k1 + k2)
    f = lambda p: (p / (1 - p)) * exp_single / exp_multi
    p0 = exp_multi / (exp_single + exp_multi)
    pc = min(1.0, 2 * min(float(_binom.sf(k2 - 1, k1 + k2, p0)), float(_binom.cdf(k2, k1 + k2, p0))))
    return dict(expected_single=round(exp_single, 1), expected_multi=round(exp_multi, 2), expected_share_two_or_more=round(p0, 3),
                ratio_of_ratios=round((k2 / exp_multi) / (k1 / exp_single), 2), ci95=[round(f(lo), 2), round(f(hi), 2)], p_central=float("%.2g" % pc))


e1a, e2a = hh_expected(hh_lone)
e1b, e2b = hh_expected(hh_lone + hh_unstated * hh_lone / known)
assert abs(s["single"] / (exp_joint / f1) - rev_checks["all_allowances_compounded"]["single_ratio"]) < 1e-3
assert abs(s["multi"] / (e2a / f2) - rev_checks["all_allowances_compounded"]["multi_ratio"]) < 1e-3
sep = dict(
    note="ratio of the multiple-death rate ratio to the single-death rate ratio (Equation 8) after each allowance; exact interval and central p "
         "conditional on Singapore's total, with the English expectations taken as fixed",
    household=ror_after(e1a, e2a),
    household_unstated_spread=ror_after(e1b, e2b),
    household_and_age=ror_after(exp_joint, e2a),
    compounded=ror_after(exp_joint / f1, e2a / f2),
    compounded_unstated_spread=ror_after(e1b / f1, e2b / f2))


def with_english_count_error(n=400000, seed=20261121):
    """the compounded ratios with the count errors of the English tables as well as Singapore's (Jeffreys gamma and beta draws)"""
    g = np.random.default_rng(seed)
    lo_ = g.gamma(lone_old + 0.5, 1, n); ly_ = g.gamma(lone_young + 0.5, 1, n)
    sm_ = g.gamma(single_in_multi + 0.5, 1, n); mm_ = g.gamma(rev["multi"] + 0.5, 1, n)
    x1 = lo_ / hy_old_e * hy_old_s + ly_ / (hy_lone_e - hy_old_e) * (hy_lone_s - hy_old_s) + sm_ / hy_multi_e * hy_multi_s
    x2 = mm_ / hy_multi_e * hy_multi_s
    q1 = g.beta(strict_p["single"] + 0.5, base_p["single"] - strict_p["single"] + 0.5, n)      # strict counts are a subset of the base counts
    q2 = g.beta(strict_p["multi"] + 0.5, base_p["multi"] - strict_p["multi"] + 0.5, n)
    k1_ = g.gamma(s["single"] + 0.5, 1, n); k2_ = g.gamma(s["multi"] + 0.5, 1, n)
    single = k1_ / (x1 * q1); multi = k2_ / (x2 * q2)
    pct = lambda x: [round(float(np.percentile(x, q)), 2) for q in (2.5, 50, 97.5)]
    return dict(note="2.5, 50 and 97.5 percentiles; all counts drawn, exposures fixed; an approximate check on the exact intervals, which carry Singapore's count error only",
                single_ratio=pct(single), multi_ratio=pct(multi), ratio_of_ratios=pct(multi / single),
                P_ratio_of_ratios_above_1=round(float(np.mean(multi / single > 1)), 3), P_single_ratio_below_1=round(float(np.mean(single < 1)), 4))


sep["compounded_with_english_count_error"] = with_english_count_error()
rev_checks["separation_after_allowances"] = sep
rev_checks["per_person_multi_ratio_range"] = sorted(round(cmpx[w]["compare_" + b]["person_years_M"]["multi"]["ratio"]["ratio"], 2)
                                                    for w in ("primary", "overlap", "review") for b in ("residential", "dwelling_stock", "household_units"))
rev_checks["per_dwelling_multi_ratio_range"] = sorted(round(cmpx[w]["compare_" + b]["dwelling_years_M"]["multi"]["ratio"]["ratio"], 2)
                                                      for w in ("primary", "overlap", "review") for b in ("residential", "dwelling_stock", "household_units"))
out["review_checks"] = rev_checks

# ---- 9c. what incomplete capture and an empty tail can and cannot show --------------------------
# (i) capture: if a single-death fire enters a record with probability c1 and a multiple-death fire
# with probability c2 whatever its size, the recorded odds of the share are the true odds times
# c2/c1; the rate of multiple-death fires depends on c2 alone and the shape above one death is unchanged.
_or = cmpx["primary"]["compare_dwelling_stock"]["identity"]["odds_ratio"]
_u = 0
while s["multi"] / (s["single"] + s["multi"] + _u) >= e["multi"] / (e["single"] + e["multi"]):
    _u += 1
rev_checks["capture_tipping"] = dict(
    note="capture of single-death fires relative to multiple-death fires, as a fraction of England's, that would alone produce the crude "
         "separation of the two rate ratios; and the unrecorded single-death fires that would bring Singapore's share down to England's",
    relative_capture_needed=round(1 / _or, 3), unrecorded_single_death_fires_needed=_u,
    single_ratio_if_one_fifth_of_single_death_fires_unrecorded=round(
        cmpx["primary"]["compare_dwelling_stock"]["dwelling_years_M"]["single"]["ratio"]["ratio"] / 0.8, 3),
    share_with_them=round(s["multi"] / (s["single"] + s["multi"] + _u), 4), england_share=round(e["multi"] / (e["single"] + e["multi"]), 4))
# (ii) the share as a monitoring signal: halving the single-death fires of the register, nothing else changed
_c = counts(sg_bld_d)
_deaths = sum(k * v for k, v in sg_bld_d.items())
_half = _c["single"] // 2
rev_checks["share_as_monitoring_signal"] = dict(
    single=_c["single"], multi=_c["multi"], deaths=_deaths, share=round(_c["multi"] / _c["fatal_fires"], 3),
    after_halving_single_death_fires=dict(single=_c["single"] - _half, deaths=_deaths - _half,
                                          share=round(_c["multi"] / (_c["fatal_fires"] - _half), 3)))
# (iii) statements about the tail that need no model
_T = 14.0
rev_checks["model_free_tail"] = dict(
    fires_of_5_or_more=sum(v for k, v in sg_bld_d.items() if k >= 5), years=_T,
    one_sided_95pct_upper_rate_per_yr=round(-np.log(0.05) / _T, 3), return_period_at_least_yr=round(_T / -np.log(0.05), 1),
    likelihood_factor_for_1pct_beyond_the_maximum=round(0.99 ** _c["fatal_fires"], 3), fatal_fires=_c["fatal_fires"],
    note="the likelihood of a record with no fire above four deaths is multiplied by (1 - eps)^K when probability eps is moved "
         "above four deaths, wherever it is placed")

# ---- 9d. sizes of English multiple-death fires bounded from the official counts alone -----------
# Table FIRE0502 gives fatal fires and deaths for every quarter, fire and rescue service and location
# type. With the number of multiple-death fires of each year and location (official fatal fires minus
# the flagged single-death fires), these counts bound the number of fires that reach n deaths whatever
# the allocation of deaths to fires. No grouping of records is used.
CELLS = rd("england_fire0502_cells.csv", "utf-8")
LOCNAME = {"dwellings": "Dwellings", "other buildings": "Other Buildings", "road vehicles": "Road Vehicles", "outdoors": "Other Outdoors"}


def cell_options(k, d, n):
    """a cell of k fatal fires and d deaths: {r: (fewest, most)} fires with n or more deaths when r fires have two or more"""
    ex = d - k
    if ex == 0:
        return {0: (0, 0)}
    opts = {}
    for r in range(1, min(k, ex) + 1):
        x = ex - r                      # deaths left after one for every fire and a second for each multiple-death fire
        opts[r] = (0 if x <= r * (n - 3) else 1, min(r, x // (n - 2)))
    return opts


def stratum_bounds(cells, m_total, n):
    """fewest and most fires with n or more deaths over every allocation with m_total multiple-death fires"""
    lo, hi = {0: 0}, {0: 0}
    for k, d in cells:
        nlo, nhi = {}, {}
        for r, (a, b) in cell_options(k, d, n).items():
            for t in lo:
                if t + r <= m_total:
                    nlo[t + r] = min(nlo.get(t + r, 10 ** 9), lo[t] + a)
                    nhi[t + r] = max(nhi.get(t + r, -1), hi[t] + b)
        lo, hi = nlo, nhi
    assert m_total in lo, "no allocation of deaths to fires matches the official totals"
    return lo[m_total], hi[m_total]


def margin_bounds(locs, fys):
    res = dict(multiple_death_fires=0, fatal_fires=0, deaths=0)
    tot = {n: [0, 0] for n in (3, 4, 5)}
    strata = 0; outside = []
    for loc in locs:
        for fy in fys:
            a = next(x for x in OFF if x["location"] == LOCNAME[loc] and x["financial_year"] == fy)
            cells = [(int(c_["fatal_fires"]), int(c_["deaths"])) for c_ in CELLS if c_["location"] == loc and c_["financial_year"] == fy]
            assert sum(k for k, d in cells) == int(a["events"]) and sum(d for k, d in cells) == int(a["deaths"]), (loc, fy)
            m_ = int(a["events_multiple"])
            res["multiple_death_fires"] += m_; res["fatal_fires"] += int(a["events"]); res["deaths"] += int(a["deaths"])
            # the reconstruction of this one year and location, to be compared with the bounds of the same stratum
            rec_s = eng_dist([r for r in HO if r["FINANCIAL_YEAR"] == fy and r["INCIDENT_LOCATION_TYPE"] == LOCNAME[loc]])
            strata += 1
            for n in (3, 4, 5):
                lo_, hi_ = stratum_bounds(cells, m_, n)
                tot[n][0] += lo_; tot[n][1] += hi_
                got = sum(v for k, v in rec_s.items() if k >= n)
                if not lo_ <= got <= hi_:
                    outside.append(dict(location=loc, financial_year=fy, threshold=n, bounds=[lo_, hi_], reconstruction=got))
    res["strata_checked"] = strata
    res["strata_with_reconstruction_outside_bounds"] = outside
    cutname = "dwellings" if locs == ["dwellings"] else ("buildings" if locs == ["dwellings", "other buildings"] else "all fires")
    rec = eng(cutname, fys)             # grouping reconstruction on the same years, Grenfell Tower included
    M_ = res["multiple_death_fires"]
    res["bounds"] = {"%d or more" % n: tot[n] for n in (3, 4, 5)}
    res["reconstruction"] = {"%d or more" % n: sum(v for k, v in rec.items() if k >= n) for n in (3, 4, 5)}
    res["reconstruction_multiple_death_fires"] = sum(v for k, v in rec.items() if k >= 2)
    res["reconstruction_inside_bounds"] = all(tot[n][0] <= res["reconstruction"]["%d or more" % n] <= tot[n][1] for n in (3, 4, 5))
    res["share_exactly_two_bounds"] = [round(1 - tot[3][1] / M_, 3), round(1 - tot[3][0] / M_, 3)]
    res["share_exactly_two_reconstruction"] = round(1 - res["reconstruction"]["3 or more"] / res["reconstruction_multiple_death_fires"], 3)
    return res


out["england_bounds_from_official_counts"] = {
    "note": "sharp bounds for each threshold separately; one allocation need not attain them at every threshold together. "
            "Grenfell Tower is included in the counts of 2017/18.",
    "dwellings, 2012/13 to 2023/24": margin_bounds(["dwellings"], ENG_W["overlap"]),
    "buildings, 2012/13 to 2023/24": margin_bounds(["dwellings", "other buildings"], ENG_W["overlap"]),
    "dwellings, 2010/11 to 2018/19": margin_bounds(["dwellings"], ENG_W["review"]),
    "dwellings, 2010/11 to 2025/26": margin_bounds(["dwellings"], ENG_W["all"]),
    "buildings, 2010/11 to 2025/26": margin_bounds(["dwellings", "other buildings"], ENG_W["all"]),
    "all fires, 2010/11 to 2018/19": margin_bounds(["dwellings", "other buildings", "road vehicles", "outdoors"], ENG_W["review"]),
}

# ---- 9e. position of victims: checks on the comparison of single- and multiple-death fires -------
_origin = lambda r: r["VICTIM_LOCATION_START"] == "Room, cabin or compartment of fire origin"
_away = lambda r: r["VICTIM_LOCATION_START"].startswith("Different floor")
_pairs = [g for g in Counter(tuple(r[k] for k in KEY) for r in MV).items() if g[1] == 2]
_pair_keys = {g[0] for g in _pairs}
_two = [r for r in MV if tuple(r[k] for k in KEY) in _pair_keys]
_no1718 = lambda sel: [r for r in sel if r["FINANCIAL_YEAR"] != "2017/18"]
out["two_processes_england_dwellings"]["position_checks"] = dict(
    note="position when the fire started; 'elsewhere' is a different floor, room, cabin or compartment from the fire origin; "
         "every victim stays in the denominator",
    single=dict(n=len(SV), origin=sum(map(_origin, SV)), elsewhere=sum(map(_away, SV)),
                other_or_not_known=len(SV) - sum(map(_origin, SV)) - sum(map(_away, SV)),
                share_elsewhere=round(sum(map(_away, SV)) / len(SV), 3),
                share_elsewhere_if_every_other_or_unknown_position_is_counted=round(1 - sum(map(_origin, SV)) / len(SV), 3)),
    multiple=dict(n=len(MV), origin=sum(map(_origin, MV)), elsewhere=sum(map(_away, MV)), share_elsewhere=round(sum(map(_away, MV)) / len(MV), 3)),
    fires_with_exactly_two_deaths=dict(fires=len(_pairs), victims=len(_two), elsewhere=sum(map(_away, _two)),
                                       share_elsewhere=round(sum(map(_away, _two)) / len(_two), 3),
                                       note="rebuilt by grouping, so conditional on the grouping; every fire has the same weight"),
    without_2017_18=dict(single=round(sum(map(_away, _no1718(SV))) / len(_no1718(SV)), 3),
                         multiple=round(sum(map(_away, _no1718(MV))) / len(_no1718(MV)), 3)))

# ---- 9f. a test of the single law where the answer is known --------------------------------------
# English dwelling fires, Grenfell Tower set aside. The number of fires with 3, 4 and 5 or more deaths is
# known from the reconstruction and, without any grouping, inside bounds from the official counts.
# (a) each form fitted and checked on five windows; (b) fitted on one half of the record and checked on the
# other; (c) the same multiple-death fires with and without the single-death fires of one-person households;
# (d) records of Singapore's size drawn from the English record.
from scipy.special import zeta as _zeta
from scipy.stats import poisson as _poisson

DWO = [r for r in HO if r["INCIDENT_LOCATION_TYPE"] == "Dwellings" and not is_grenfell(r)]


def _ge(d, n):
    return sum(v for k, v in d.items() if k >= n)


def fit_single_law(d):
    """one power law on 1, 2, 3, ...: P(N = n) proportional to n^-b (no cap)"""
    ns = np.array(sorted(d), dtype=float); cs = np.array([d[int(n)] for n in ns], dtype=float)
    return float(minimize_scalar(lambda b: b * float((cs * np.log(ns)).sum()) + cs.sum() * np.log(_zeta(b, 1)),
                                 bounds=(1.05, 9.0), method="bounded").x)


def fit_shape_above_one(d, upper=12.0):
    """power law on 2, 3, 4, ... fitted to the multiple-death fires alone (no cap)"""
    ns = np.array(sorted(k for k in d if k >= 2), dtype=float); cs = np.array([d[int(n)] for n in ns], dtype=float)
    return float(minimize_scalar(lambda b: b * float((cs * np.log(ns)).sum()) + cs.sum() * np.log(_zeta(b, 2)),
                                 bounds=(1.05, upper), method="bounded").x)


def predict_tail(d, thresholds=(3, 4, 5)):
    """fires expected at each threshold under each form fitted to d"""
    K_ = sum(d.values()); M_ = _ge(d, 2)
    b_ = fit_single_law(d); b2_ = fit_shape_above_one(d)
    return dict(fatal_fires=round(K_, 1), single_death_fires=round(d[1], 1), multiple_death_fires=M_, share=round(M_ / K_, 4),
                b_single=round(b_, 2), b2=round(b2_, 2),
                single_law={str(n): round(K_ * float(_zeta(b_, n) / _zeta(b_, 1)), 1) for n in thresholds},
                three_part={str(n): round(M_ * float(_zeta(b2_, n) / _zeta(b2_, 2)), 1) for n in thresholds})


def ordinary_bounds(fys, thresholds=(3, 4, 5)):
    """dwelling fires with n or more deaths, from official counts alone; the one fire of 71 deaths taken out of its cell"""
    tot = {n: [0, 0] for n in thresholds}
    for fy in fys:
        a = next(x for x in OFF if x["location"] == "Dwellings" and x["financial_year"] == fy)
        cells = [(int(c_["fatal_fires"]), int(c_["deaths"])) for c_ in CELLS if c_["location"] == "dwellings" and c_["financial_year"] == fy]
        m_ = int(a["events_multiple"])
        if fy == "2017/18":
            big = [i for i, (k, d) in enumerate(cells) if d - k >= 70]
            assert len(big) == 1, big
            cells[big[0]] = (cells[big[0]][0] - 1, cells[big[0]][1] - 71); m_ -= 1
        cells = [c_ for c_ in cells if c_[0] > 0]
        for n in thresholds:
            lo_, hi_ = stratum_bounds(cells, m_, n)
            tot[n][0] += lo_; tot[n][1] += hi_
    return {str(n): tot[n] for n in thresholds}


def against_bounds(pred, bounds):
    """inside, below or above the bounds; and how often chance alone puts a count that far from the expectation"""
    res = {}
    for n, (lo_, hi_) in bounds.items():
        x = pred[n]
        if lo_ <= x <= hi_:
            res[n] = dict(position="inside")
        elif x < lo_:
            res[n] = dict(position="below", chance=round(float(_poisson.sf(lo_ - 1, x)), 3))
        else:
            res[n] = dict(position="above", chance=round(float(_poisson.cdf(hi_, x)), 3))
    return res


def lr_single_within_three_part(d):
    """likelihood-ratio test of the single law inside the three-part form (one degree of freedom)"""
    ns = np.array(sorted(d), dtype=float); cs = np.array([d[int(n)] for n in ns], dtype=float)
    b_ = fit_single_law(d)
    ll1 = float(-(b_ * (cs * np.log(ns)).sum() + cs.sum() * np.log(_zeta(b_, 1))))
    K_ = cs.sum(); M_ = float(sum(v for k, v in d.items() if k >= 2)); q_ = M_ / K_
    b2_ = fit_shape_above_one(d)
    m_ns = ns[ns >= 2]; m_cs = cs[ns >= 2]
    ll3 = float((K_ - M_) * np.log(1 - q_) + M_ * np.log(q_) - (b2_ * (m_cs * np.log(m_ns)).sum() + M_ * np.log(_zeta(b2_, 2))))
    return float("%.2g" % chi2.sf(max(0.0, 2 * (ll3 - ll1)), 1))


KT = dict(note="English dwelling fires, Grenfell Tower set aside; fires expected at 3, 4 and 5 or more deaths under a single power law "
               "fitted to all fatal fires and under the three-part form; 'chance' is the Poisson probability of a count at least as far "
               "from the expectation as the nearer official bound")
kt_windows = {"2010/11 to 2025/26": FYS, "2010/11 to 2018/19": FYS[:9], "2012/13 to 2023/24": FYS[2:14],
              "2010/11 to 2017/18": FYS[:8], "2018/19 to 2025/26": FYS[8:]}
KT["fitted_and_checked_on_the_same_window"] = {}
for name_, fys_ in kt_windows.items():
    d_ = eng_dist([r for r in DWO if r["FINANCIAL_YEAR"] in fys_])
    p_ = predict_tail(d_); bw_ = ordinary_bounds(fys_)
    p_.update(reconstruction={str(n): _ge(d_, n) for n in (3, 4, 5)}, official_bounds=bw_,
              single_law_against_bounds=against_bounds(p_["single_law"], bw_), three_part_against_bounds=against_bounds(p_["three_part"], bw_))
    KT["fitted_and_checked_on_the_same_window"][name_] = p_
assert KT["fitted_and_checked_on_the_same_window"]["2010/11 to 2018/19"]["b_single"] == 4.41
assert KT["fitted_and_checked_on_the_same_window"]["2010/11 to 2018/19"]["b2"] == 3.87
kw_ = KT["fitted_and_checked_on_the_same_window"]
KT["windows_single_law_below_bounds"] = sum(1 for v in kw_.values() if any(x["position"] == "below" for x in v["single_law_against_bounds"].values()))
KT["windows_three_part_outside_bounds"] = sum(1 for v in kw_.values() if any(x["position"] != "inside" for x in v["three_part_against_bounds"].values()))


def split_test(fit_fys, check_fys):
    df_ = eng_dist([r for r in DWO if r["FINANCIAL_YEAR"] in fit_fys]); dc_ = eng_dist([r for r in DWO if r["FINANCIAL_YEAR"] in check_fys])
    b_ = fit_single_law(df_); b2_ = fit_shape_above_one(df_)
    Kc, Mc = sum(dc_.values()), _ge(dc_, 2)
    s1 = {str(n): round(Kc * float(_zeta(b_, n) / _zeta(b_, 1)), 1) for n in (3, 4, 5)}
    s3 = {str(n): round(Mc * float(_zeta(b2_, n) / _zeta(b2_, 2)), 1) for n in (3, 4, 5)}
    bw_ = ordinary_bounds(check_fys)
    return dict(fitted_on=fit_fys[0] + " to " + fit_fys[-1], checked_on=check_fys[0] + " to " + check_fys[-1], b_single=round(b_, 2), b2=round(b2_, 2),
                single_law=s1, three_part_given_the_later_count_of_multiple_death_fires=s3,
                reconstruction={str(n): _ge(dc_, n) for n in (3, 4, 5)}, official_bounds=bw_,
                single_law_against_bounds=against_bounds(s1, bw_), three_part_against_bounds=against_bounds(s3, bw_))


KT["fitted_on_first_half_checked_on_second"] = split_test(FYS[:8], FYS[8:])
KT["fitted_on_second_half_checked_on_first"] = split_test(FYS[8:], FYS[:8])

# (c) the same multiple-death fires with and without the single-death fires of one-person households.
# A one-person household can have only a single-death fire, so leaving those households out removes
# single-death fires and no multiple-death fire. Home Office review, Table 38, 2010/11 to 2018/19.
rev_d = eng_dist([r for r in DWO if r["FINANCIAL_YEAR"] in FYS[:9]])
m_review = _ge(rev_d, 2) + 1                           # the review's fatal fires include Grenfell Tower
hh_cases = {"all households": rev_d[1],
            "households of two or more, unstated type counted with them": hh_total - hh_lone - m_review,
            "households of two or more, unstated type spread in proportion": hh_total - (hh_lone + hh_unstated * hh_lone / known) - m_review}
KT["households"] = dict(window="2010/11 to 2018/19", reconstruction={str(n): _ge(rev_d, n) for n in (3, 4, 5)},
                        official_bounds=ordinary_bounds(FYS[:9]), cases={})
for name_, singles_ in hh_cases.items():
    d_ = dict(rev_d); d_[1] = singles_
    p_ = predict_tail(d_)
    p_["single_law_against_bounds"] = against_bounds(p_["single_law"], KT["households"]["official_bounds"])
    p_["p_single_law_within_three_part"] = lr_single_within_three_part(d_)
    KT["households"]["cases"][name_] = p_
# the premise that every multiple-death fire is in a household of two or more: what if one in ten is not
# (a visitor in a one-person household), in proportion at every size
sens = {}
for name_, singles_ in list(hh_cases.items())[1:]:
    moved = 0.1 * _ge(rev_d, 2)
    d_ = {k: (v * 0.9 if k >= 2 else singles_ + moved) for k, v in rev_d.items()}
    b_ = fit_single_law(d_); b2_ = fit_shape_above_one(d_)
    sens[name_] = dict(single_law_5_or_more=round(sum(d_.values()) * float(_zeta(b_, 5) / _zeta(b_, 1)), 1),
                       three_part_5_or_more=round(_ge(d_, 2) * float(_zeta(b2_, 5) / _zeta(b2_, 2)), 1), b_single=round(b_, 2))
KT["households"]["if_one_in_ten_multiple_death_fires_is_in_a_one_person_household"] = sens
KT["households"]["edition_note"] = dict(review_table_fatal_dwelling_fires=hh_total, record_level_fatal_dwelling_fires=int(sum(rev_d.values())) + 1,
                                        note="the household rows take single-death fires from the review table; both totals include Grenfell Tower")
c_ = KT["households"]["cases"]
KT["households"]["single_law_5_or_more_ratio_to_all_households"] = [
    round(c_[k]["single_law"]["5"] / c_["all households"]["single_law"]["5"], 1) for k in list(c_)[1:]]
# the curve between those points: the single law as single-death fires are removed and nothing else changes
keep_grid = [round(x, 3) for x in np.geomspace(1.0, 0.1, 31)]
curve = []
for f_ in keep_grid:
    d_ = dict(rev_d); d_[1] = rev_d[1] * f_
    b_ = fit_single_law(d_)
    curve.append(dict(kept=f_, b_single=round(b_, 3), single_law_5_or_more=round(sum(d_.values()) * float(_zeta(b_, 5) / _zeta(b_, 1)), 2)))
KT["households"]["curve_as_single_death_fires_are_removed"] = curve
d_ = dict(rev_d); d_[1] = rev_d[1] * cmpx["primary"]["compare_dwelling_stock"]["dwelling_years_M"]["single"]["ratio"]["ratio"]
KT["households"]["at_singapores_single_death_rate"] = dict(
    kept=round(cmpx["primary"]["compare_dwelling_stock"]["dwelling_years_M"]["single"]["ratio"]["ratio"], 3), b_single=round(fit_single_law(d_), 2),
    singapore_b_single_dwelling_stock=round(fit_single_law({int(k): v for k, v in cmpx["primary"]["singapore_dwelling_stock"]["distribution"].items()}), 2))

# (d) records of Singapore's size, simulated two ways. In both a record is 42 fatal fires in fourteen years.
#   English shape: each fire has one death with probability 34/42; otherwise its size is drawn from the English
#                  multiple-death dwelling fires. The single law is then not the true form.
#   Single law:    every size is drawn from one power law, with the exponent fitted to Singapore's own record
#                  on the matched basis. The single law is then the true form.
# The return period of a fire of five or more deaths is known in each case, and both forms are fitted to every record.
full_d = eng_dist(DWO)
sizes_ = np.array([k for k, v in full_d.items() if k >= 2 for _ in range(v)])
n1_, n2_, yrs_ = s["single"], s["multi"], 14.0
n_all = n1_ + n2_
rng_kt = np.random.default_rng(20261122)
pc = lambda x, ps: [round(float(np.percentile(x, p)), 2) for p in ps]


def fit_both(draw):
    d_ = Counter(int(x) for x in draw); d_.setdefault(1, 0)
    b_ = fit_single_law(d_)
    r1 = yrs_ / (len(draw) * float(_zeta(b_, 5) / _zeta(b_, 1)))
    multi_ = draw[draw >= 2]
    if len(multi_) == 0 or (multi_ == 2).all():
        return b_, r1, float("inf")                        # no multiple-death fire above two deaths: no finite estimate
    b2_ = fit_shape_above_one(d_, upper=40.0)
    return b_, r1, yrs_ / (len(multi_) * float(_zeta(b2_, 5) / _zeta(b2_, 2)))


def simulate(draw_record, truth, draws=4000):
    res = [fit_both(draw_record()) for _ in range(draws)]
    b1 = np.array([r[0] for r in res]); r1 = np.array([r[1] for r in res]); r3 = np.array([r[2] for r in res])
    within = lambda x: round(float(np.mean((x > truth / 2) & (x < truth * 2))), 2)
    return dict(fatal_fires=n_all, years=yrs_, draws=draws, true_return_period_5_or_more_yr=round(truth, 1),
                single_law_exponent_pct_2_5_50_97_5=[round(float(np.percentile(b1, p)), 2) for p in (2.5, 50, 97.5)],
                single_law_return_period_yr_pct_2_5_25_50_75_97_5=pc(r1, (2.5, 25, 50, 75, 97.5)),
                three_part_return_period_yr_pct_2_5_25_50_75=pc(r3, (2.5, 25, 50, 75)),
                three_part_share_of_records_with_no_finite_estimate=round(float(np.mean(np.isinf(r3))), 3),
                share_within_a_factor_of_two_of_the_truth=dict(single_law=within(r1), three_part=within(r3)))


def record_english_shape():
    m_ = rng_kt.binomial(n_all, n2_ / n_all)
    return np.concatenate([np.ones(n_all - m_, dtype=int), rng_kt.choice(sizes_, m_, replace=True)]).astype(int)


true_rp = yrs_ / (n2_ * _ge(full_d, 5) / _ge(full_d, 2))
KT["singapore_sized_english_records"] = simulate(record_english_shape, true_rp)
KT["singapore_sized_english_records"].update(expected_single_death_fires=n1_, expected_multiple_death_fires=n2_)

b_true = fit_single_law({int(k): v for k, v in cmpx["primary"]["singapore_dwelling_stock"]["distribution"].items()})
support = np.arange(1, 10 ** 4 + 1); w_ = support ** (-b_true); w_ = w_ / w_.sum()
true_rp_b = yrs_ / (n_all * float(_zeta(b_true, 5) / _zeta(b_true, 1)))
KT["records_drawn_from_a_single_power_law"] = simulate(lambda: rng_kt.choice(support, n_all, replace=True, p=w_), true_rp_b)
KT["records_drawn_from_a_single_power_law"]["exponent"] = round(b_true, 2)
out["known_answer_test"] = KT

# ---- 10. values for the table of fitted laws ----------------------------------------------------
from fs_model import fit_b, gof_pboot, profile_ci


def single_law_row(dist):
    bins = [(k, k, v) for k, v in sorted(dist.items())]
    b, _ = fit_b(bins, 10 ** 4)
    ci = profile_ci(bins, 10 ** 4, grid=np.linspace(0.5, 8, 3000))
    nmax = max(dist)
    _, p = gof_pboot(bins, nmax, fit_b(bins, nmax)[0])
    return dict(fatal_fires=sum(dist.values()), b=round(b, 2), ci95=[round(x, 2) for x in ci], gof_p=round(p, 2))


t3 = {"England, dwellings 2010/11 to 2025/26": single_law_row(ordinary(eng("dwellings", ENG_W["all"]))),
      "England, dwellings 2010/11 to 2018/19 (shape)": shape(ordinary(eng("dwellings", ENG_W["review"])), "England 2010/11 to 2018/19, dwellings")}
for key in ("England_all_fires_within_regime", "China_all_fires_within_regime"):
    bins = [x for x in BASES[key] if x[0] >= 2]
    m = sum(c for _, _, c in bins); two = sum(c for lo, hi, c in bins if lo == 2 and hi == 2)
    b, ll, ci = fit_tail(bins)
    t3[key] = dict(multi=m, two=two, share_exactly_two=round(two / m, 3), b2=round(b, 2), b2_ci95=[round(x, 2) for x in ci])
i2 = np.arange(2, 10 ** 4 + 1, dtype=float)
from scipy.optimize import brentq
t3["United States, implied b2 from published share"] = round(float(brentq(lambda b: 2.0 ** (-b) / np.sum(i2 ** (-b)) - 0.768, 1.01, 14)), 2)
t3["Singapore own shape, return period 30 or more yr"] = round(1 / (lam2_sg_bld * tail_survival(own_b2, 30)), 0)
yc = [sum(1 for r in SG_W["primary"]["rows"] if r["y"] == y) for y in range(2012, 2026)]
t3["Singapore annual counts"] = yc
t3["Singapore variance to mean of annual counts"] = round(float(np.var(yc, ddof=1) / np.mean(yc)), 2)
out["table_extras"] = t3

io.open(os.path.join(RESULTS, "r2_comparison.json"), "w", encoding="utf-8").write(json.dumps(out, indent=1))
print("wrote results/r2_comparison.json")
