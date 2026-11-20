"""Rewrite the Singapore rows of the two bins files from the register.

The register (data/fatal_incident_register.csv) is the single source of the Singapore counts.
data/international_severity_bins.csv and data/revision_bases.csv carry a copy of those counts
beside the other datasets; this script refreshes the copy after any change to the register.

Run:  python analysis/frequency_severity/sync_singapore_bins.py
"""
import csv
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
REG = os.path.join(REPO, "data", "fatal_incident_register.csv")


def bins(basis):
    c = {}
    for r in csv.DictReader(io.open(REG, encoding="utf-8-sig")):
        if r["building_scope"].strip().lower() != "building" or not (2012 <= int(r["year"]) <= 2025):
            continue
        p = r["premises_type"].lower()
        res = ("hdb" in p or "residen" in p or "condo" in p) or r["incident_id"] == "SGF-2014-04"
        if basis == "residential" and not res:
            continue
        n = int(r["fatalities"])
        c[n] = c.get(n, 0) + 1
    return sorted(c.items())


def rewrite(path, key_col, mapping):
    rows = list(csv.DictReader(io.open(path, encoding="utf-8")))
    cols = list(rows[0].keys())
    out, done = [], set()
    for r in rows:
        k = r[key_col]
        if k in mapping:
            if k not in done:
                for n, cnt in bins(mapping[k]):
                    out.append(dict(r, bin_lo=str(n), bin_hi=str(n), events=str(cnt)))
                done.add(k)
            continue
        out.append(r)
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(os.path.basename(path), {k: bins(v) for k, v in mapping.items()})


rewrite(os.path.join(REPO, "data", "international_severity_bins.csv"), "jurisdiction", {"Singapore": "buildings"})
rewrite(os.path.join(REPO, "data", "revision_bases.csv"), "dataset",
        {"Singapore_residential": "residential", "Singapore_buildings": "buildings"})
