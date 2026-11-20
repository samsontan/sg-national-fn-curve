"""England, observed deaths per fire by location, from Home Office record-level data.

Source: Home Office, Fire statistics incident level datasets, "Fire-related fatalities dataset"
(England, year ending March 2026 edition; one row per fire-related fatality; Open Government
Licence v3.0). Saved here as data/home_office_fire_related_fatalities.csv.

The dataset has no incident identifier. It flags each fatality as belonging to a single- or
multiple-fatality incident. Single-fatality fires are therefore counted directly. Multiple-
fatality fires are reconstructed by grouping the flagged rows on the incident-level fields
(year, hour band, location type, intent, alarm, cause, ignition source, start location).
The reconstruction is validated against the Home Office in-depth review's published all-fires
distribution for 2010/11-2018/19.

Writes results/england_incident_level.json.
Run:  python analysis/frequency_severity/england_incident_level.py
"""
import csv
import io
import json
import os
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import beta, fisher_exact

from fs_model import REPO, RESULTS, fit_b, gof_pboot, profile_ci, survival

SRC = os.path.join(REPO, "data", "home_office_fire_related_fatalities.csv")
KEY = ["FINANCIAL_YEAR", "HOUR_BANDS", "INCIDENT_LOCATION_TYPE", "ACCIDENTAL_OR_DELIBERATE", "ALARM_SYSTEM",
       "CAUSE_OF_FIRE", "SOURCE_OF_IGNITION", "FIRE_START_LOCATION"]
WINDOW = ["2010/11", "2011/12", "2012/13", "2013/14", "2014/15", "2015/16", "2016/17", "2017/18", "2018/19"]
T, POP, DWELLINGS = 9.0, 54.35, 23_635_441
PUBLISHED_ALL_FIRES = {1: 2213, 2: 108, 3: 17, 4: 8, 5: 5, 6: 3, 11: 1, 71: 1}   # in-depth review, Table 2

rows = [{k: (v or "").strip() for k, v in r.items()} for r in csv.DictReader(io.open(SRC, encoding="utf-8"))]
MULTI = "Multiple fatality incident"


def cp(k, n):
    return [round(float(beta.ppf(.025, k, n - k + 1)) if k else 0.0, 3), round(float(beta.ppf(.975, k + 1, n - k)), 3)]


def incidents(sel):
    """size distribution {deaths per fire: number of fires} for a selection of fatality rows"""
    dist = Counter()
    dist[1] = sum(1 for r in sel if r["FATALITY_TOTAL"] != MULTI)
    groups = Counter(tuple(r[k] for k in KEY) for r in sel if r["FATALITY_TOTAL"] == MULTI)
    for size in groups.values():
        dist[size] += 1
    return dict(sorted(dist.items()))


def summary(sel):
    d = incidents(sel)
    n = sum(d.values()); multi = n - d.get(1, 0)
    return dict(deaths=len(sel), fatal_fires=n, multi_fatality_fires=multi, distribution=d,
                mfi=round(multi / n, 4), mfi_ci95=cp(multi, n))


win = [r for r in rows if r["FINANCIAL_YEAR"] in WINDOW]
out = {"source": "Home Office fire-related fatalities dataset, England, year ending March 2026 edition",
       "window": "2010/11-2018/19", "rows_in_window": len(win)}
out["all_fires_reconstructed"] = incidents(win)
out["all_fires_published"] = PUBLISHED_ALL_FIRES
by_loc = {loc: summary([r for r in win if r["INCIDENT_LOCATION_TYPE"] == loc])
          for loc in ["Dwellings", "Other Buildings", "Road Vehicles", "Other Outdoors"]}
out["by_location"] = by_loc
dw = [r for r in win if r["INCIDENT_LOCATION_TYPE"] == "Dwellings"]
bld = [r for r in win if r["INCIDENT_LOCATION_TYPE"] in ("Dwellings", "Other Buildings")]
out["buildings"] = summary(bld)
out["dwellings_by_intent"] = {a: summary([r for r in dw if r["ACCIDENTAL_OR_DELIBERATE"] == a])
                              for a in ("Recorded as Accidental", "Recorded as Deliberate")}
d = by_loc["Dwellings"]
out["dwelling_rates"] = dict(per_million_population_yr=round(d["fatal_fires"] / T / POP, 2),
                             per_dwelling_yr=float(f"{d['fatal_fires'] / T / DWELLINGS:.2e}"))
# sensitivity: the reconstruction has two more four-death and two fewer two-death fires than the published
# all-fires table, i.e. two pairs of two-death fires share every grouping field. Splitting them:
out["dwellings_mfi_if_two_fours_are_split"] = round((d["multi_fatality_fires"] + 2) / (d["fatal_fires"] + 2), 4)

# severity law of ordinary dwelling fires (Grenfell Tower, N = 71, set aside) and the Grenfell test
ordinary = [(n, n, c) for n, c in d["distribution"].items() if n <= 6]
b, _ = fit_b(ordinary, 10 ** 4)
ci = profile_ci(ordinary, 10 ** 4, grid=np.linspace(0.5, 8, 3000))
_, gof = gof_pboot(ordinary, 6, fit_b(ordinary, 6)[0])
lam = d["fatal_fires"] / T
out["dwelling_severity_law"] = dict(
    ordinary_fires=sum(c for _, _, c in ordinary), b=round(b, 2), ci95=[round(v, 2) for v in ci], gof_p=round(gof, 2),
    fatal_dwelling_fires_per_yr=round(lam, 1),
    grenfell_return_period_yr=int(round(1 / (lam * survival(b, 10 ** 4, 70)), -2)),
    grenfell_return_period_yr_flat_edge=int(round(1 / (lam * survival(ci[0], 10 ** 4, 70)), -2)))

# Singapore contrast against the observed English dwelling share
from fs_model import singapore_events
_sg = singapore_events("residential")
SG = dict(k=sum(1 for r in _sg if int(r["fatalities"]) >= 2), n=len(_sg))
sg_ci = cp(SG["k"], SG["n"])
out["singapore_contrast"] = dict(
    sg_mfi=round(SG["k"] / SG["n"], 3), sg_ci95=sg_ci, england_mfi=d["mfi"], england_ci95=d["mfi_ci95"],
    ratio_point=round(SG["k"] / SG["n"] / d["mfi"], 1), ratio_at_sg_lower_limit=round(sg_ci[0] / d["mfi"], 1),
    ratio_sg_lower_over_england_upper=round(sg_ci[0] / d["mfi_ci95"][1], 1),
    fisher_p=float(f"{fisher_exact([[SG['k'], SG['n'] - SG['k']], [d['multi_fatality_fires'], d['distribution'][1]]])[1]:.1e}"))

# where victims were found, and night-time share, by single / multiple fatality (dwellings)
def share(sel, field, value):
    return round(sum(1 for r in sel if r[field] == value) / len(sel), 3)


night = {"00-01", "01-02", "02-03", "03-04", "04-05", "05-06"}
sv = [r for r in dw if r["FATALITY_TOTAL"] != MULTI]
mv = [r for r in dw if r["FATALITY_TOTAL"] == MULTI and not (r["FINANCIAL_YEAR"] == "2017/18" and r["HOUR_BANDS"] == "00-01"
      and r["FIRE_START_LOCATION"] == "Kitchen" and r["ALARM_SYSTEM"] == "Alarm Present but did not raise alarm")]
ORIGIN = "Room, cabin or compartment of fire origin"
ELSE = "Different floor, room, cabin or compartment than of fire origin"
out["dwelling_victims"] = dict(
    note="Multiple-fatality rows exclude the 71 Grenfell Tower records.",
    single=dict(n=len(sv), found_in_room_of_origin=share(sv, "VICTIM_LOCATION_FOUND", ORIGIN),
                found_elsewhere_in_building=share(sv, "VICTIM_LOCATION_FOUND", ELSE),
                night_00_06=round(sum(1 for r in sv if r["HOUR_BANDS"] in night) / len(sv), 3)),
    multiple=dict(n=len(mv), found_in_room_of_origin=share(mv, "VICTIM_LOCATION_FOUND", ORIGIN),
                  found_elsewhere_in_building=share(mv, "VICTIM_LOCATION_FOUND", ELSE),
                  night_00_06=round(sum(1 for r in mv if r["HOUR_BANDS"] in night) / len(mv), 3)))

# annual series, dwellings, 2010/11 to latest year
series = {}
for y in sorted({r["FINANCIAL_YEAR"] for r in rows}):
    s = summary([r for r in rows if r["FINANCIAL_YEAR"] == y and r["INCIDENT_LOCATION_TYPE"] == "Dwellings"])
    series[y] = dict(fatal_fires=s["fatal_fires"], multi=s["multi_fatality_fires"], mfi=s["mfi"], ci95=s["mfi_ci95"])
out["dwellings_annual"] = series
allyrs = summary([r for r in rows if r["INCIDENT_LOCATION_TYPE"] == "Dwellings"])
out["dwellings_2010_11_to_2025_26"] = dict(fatal_fires=allyrs["fatal_fires"], multi=allyrs["multi_fatality_fires"],
                                           mfi=allyrs["mfi"], ci95=allyrs["mfi_ci95"])

io.open(os.path.join(RESULTS, "england_incident_level.json"), "w", encoding="utf-8").write(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "dwellings_annual"}, indent=1))
