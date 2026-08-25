"""
FN-curve analysis for BUILDING fires in Singapore (national scope) - v4.

Round-3 revision of fn_sg_national_v3.py, addressing the adversarial attack file
("B. Attacks on the topic 1" 2026-07-14) and the international-comparison expansion.
Changes over v3:
  - TAIL STABILITY block (attack point: "12 years may hold too few multi-fatality events
    to estimate the upper tail. A single event could decide whether the curve crosses a
    criterion line."). Four quantifications:
      (a) leave-one-out jackknife on the 47 in-window events: exact F(>=N) range per N;
      (b) one-more-event perturbations (one added N=4 event; one added N=5 event);
      (c) rolling 10-year windows over the full 2005-2025 register: F(>=1..3) ranges;
      (d) criterion-crossing decisiveness at per-dwelling scale: the number of events of
          size >= N the window would need to hold for the per-dwelling curve to reach
          each criterion line - the direct answer to "a single event could decide".
  - Event-strip figure (FN_sg_eventstrip_v4): every fatal building-envelope event
    2005-2025 as a timeline dot (y = N, colour = occupancy), primary window shaded.
    Shows the sparsity the tail-stability analysis quantifies.
  - Tail-stability figure (FN_sg_tailstability_v4): national curve + exact jackknife
    envelope + hypothetical one-more-event markers.
  - International comparison figure (FN_sg_international_v4): drawn ONLY if
    data/international_fn_curves.json exists (curves verified from primary/secondary
    literature reads; per-million-population normalisation). Absent file -> skipped.
  - v4 outputs; v1/v2/v3 outputs left intact.

All v3 core logic (register load, stepped FN, seeded bootstrap, criterion overlay,
window sensitivity, classification sensitivities) is retained unchanged so that every
v3-reported number reproduces identically from this script.
"""

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)

SEEDS = {
    "fatalities_national": 20260721,
    "fatalities_residential": 20260722,
    "fatalities_non_residential": 20260723,
    "fat_plus_inj_national": 20260724,
}
N_BOOT = 10_000

REGISTER_FILE = "fatal_incident_register.csv"
RECON_FILE = "scdf_fatality_reconciliation.csv"
EXPO_FILE = "exposure_national.csv"
INTL_FILE = "international_fn_curves.json"

CRITERIA = {
    "HSE R2P2 tolerability (1e-2/N)": {"intercept": 1e-2, "slope": -1.0},
    "HSE broadly acceptable (1e-4/N)": {"intercept": 1e-4, "slope": -1.0},
    "HK PHI unacceptable (1e-3/N)": {"intercept": 1e-3, "slope": -1.0},
    "HK PHI acceptable (1e-5/N)": {"intercept": 1e-5, "slope": -1.0},
}

PRIMARY_START_YEAR = 2012
PRIMARY_END_YEAR = 2025

RESIDENTIAL_TYPES = {"HDB", "HDB_RESIDENTIAL", "PRIVATE_RESIDENTIAL", "RESIDENTIAL"}


def read_csv(path):
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def to_int(x):
    try:
        return int(float(str(x).strip()))
    except (ValueError, TypeError):
        return None


def to_float(x):
    try:
        return float(str(x).replace(",", "").strip())
    except (ValueError, TypeError):
        return None


def get_col(row, *names):
    lowered = {k.strip().lower(): v for k, v in row.items() if k}
    for n in names:
        if n in lowered and str(lowered[n]).strip() != "":
            return lowered[n]
    return None


def norm_premises(raw):
    s = (raw or "").strip().upper().replace(" ", "_").replace("-", "_")
    if not s:
        return "UNKNOWN"
    if "HDB" in s:
        return "HDB_RESIDENTIAL"
    if "RESIDEN" in s or "DWELLING" in s or "CONDO" in s or "LANDED" in s or "SHOPHOUSE_RES" in s:
        return "PRIVATE_RESIDENTIAL"
    if "COMMERC" in s or "SHOP" in s or "OFFICE" in s or "HOTEL" in s or "RETAIL" in s:
        return "COMMERCIAL"
    if "INDUST" in s or "WAREHOUSE" in s or "FACTORY" in s or "WORKSHOP" in s:
        return "INDUSTRIAL"
    return s


def load_register():
    rows = read_csv(os.path.join(DATA, REGISTER_FILE))
    events = []
    n_dropped = 0
    for r in rows:
        scope = (get_col(r, "building_scope") or "building").strip().lower()
        if scope != "building":
            n_dropped += 1
            continue
        year = to_int(get_col(r, "year"))
        fat = to_int(get_col(r, "fatalities", "fatality_count", "deaths")) or 0
        inj_raw = get_col(r, "injuries", "injury_count", "injured")
        inj = to_int(inj_raw) or 0
        ev = {
            "id": (get_col(r, "incident_id") or "").strip(),
            "year": year, "fatalities": fat, "injuries": inj,
            "injuries_blank": inj_raw is None,
            "premises": norm_premises(get_col(r, "premises_type", "building_type", "premises")),
        }
        # [E1] Occupancy-basis primary classification (see v3 header).
        if ev["id"] == "SGF-2014-04":
            ev["premises"] = "PRIVATE_RESIDENTIAL"
        events.append(ev)
    if n_dropped:
        print(f"Dropped {n_dropped} non-building-scope register rows")
    return events


def is_residential(ev):
    return ev["premises"] in RESIDENTIAL_TYPES


def fn_curve(event_sizes, T):
    sizes = np.asarray([s for s in event_sizes if s >= 1])
    if sizes.size == 0:
        return np.array([]), np.array([])
    ns = np.arange(1, sizes.max() + 1)
    freq = np.array([(sizes >= n).sum() / T for n in ns])
    return ns, freq


def bootstrap_band(event_sizes, T, ns, seed):
    rng = np.random.default_rng(seed)
    sizes = np.asarray([s for s in event_sizes if s >= 1])
    if sizes.size == 0 or ns.size == 0:
        return None, None
    k = sizes.size
    boots = np.empty((N_BOOT, ns.size))
    for b in range(N_BOOT):
        k_b = rng.poisson(k)
        sample = rng.choice(sizes, size=k_b, replace=True) if k_b > 0 else np.array([])
        boots[b] = [(sample >= n).sum() / T if sample.size else 0.0 for n in ns]
    return np.percentile(boots, 2.5, axis=0), np.percentile(boots, 97.5, axis=0)


def plot_fn(series, title, outstem, xlabel, ylabel, criteria=CRITERIA, legend_loc="below",
            band_floor=None, extra_artists=None, abcb_envelope=False, n_line_max=100):
    """series = list of dicts: {ns, freq, lo, hi, label, color, band}.
    extra_artists: optional callback(ax) to draw additional elements before legend.
    legend_loc="below" places the legend outside the axes (round-3 M2/M4/B3: an
    in-axes legend obscured the observed curve's tail in the primary figure).
    Each observed curve's terminal value is extended by a short horizontal segment
    so the last step has visible width (round-3 A-M4); captions state the convention."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 5.0), dpi=300)
    n_max = 8
    for s in series:
        ns, freq = s["ns"], s["freq"]
        if not len(ns):
            continue
        n_max = max(n_max, max(ns))
        ax.step(ns, freq, where="post", color=s.get("color", "black"),
                lw=1.6, ls=s.get("ls", "-"), label=s["label"])
        # terminal-step extension: give the last value visible width
        ax.plot([ns[-1], ns[-1] * 1.25], [freq[-1], freq[-1]],
                color=s.get("color", "black"), lw=1.6, ls=s.get("ls", "-"))
        ax.plot(ns[-1] * 1.25, freq[-1], marker="|", ms=7, mew=1.6,
                color=s.get("color", "black"), ls="none")
        if s.get("band") and s.get("lo") is not None:
            floor = band_floor if band_floor is not None else 1e-9
            ax.fill_between(ns, np.clip(s["lo"], floor, None), np.clip(s["hi"], floor, None),
                            step="post", color=s.get("color", "0.5"), alpha=0.35,
                            label="95% bootstrap band (clipped at one observable event)"
                                  if band_floor is not None else "95% bootstrap band")
    if criteria:
        n_line = np.array([1, n_line_max])
        grey_styles = [("0.25", (0, (5, 2))), ("0.45", (0, (1, 1.5))),
                       ("0.35", (0, (5, 2, 1, 2))), ("0.6", (0, (3, 1, 1, 1, 1, 1)))]
        for (name, c), (col, dash) in zip(criteria.items(), grey_styles):
            ax.plot(n_line, c["intercept"] * n_line ** c["slope"], ls=dash, color=col,
                    lw=1.0, label=name)
        ax.text(0.02, 0.03,
                "Criterion lines are single-installation criteria,\n"
                "shown for scale reference only (Section 5.2)",
                transform=ax.transAxes, fontsize=6.5, color="0.25", va="bottom",
                zorder=20, clip_on=False,
                bbox=dict(facecolor="white", edgecolor="0.7", lw=0.5, alpha=0.85))
    if abcb_envelope:
        # ABCB draft NCC 2022 societal-risk envelope (England 2020 for ABCB, Tables
        # 2/A8.2b; tabulated N=5..1000, consistent with 1e-3/N^1.5 and 1e-5/N^1.5).
        # Per-BUILDING design envelope, never adopted; drawn for scale reference only.
        n_env = np.array([5, 1000])
        n_max = max(n_max, 1000)
        for a in (1e-3, 1e-5):
            ax.plot(n_env, a * n_env ** -1.5, ls=(0, (7, 3)), color="tab:green", lw=1.1)
        ax.plot([], [], ls=(0, (7, 3)), color="tab:green", lw=1.1,
                label="ABCB draft NCC 2022 envelope, per building\n"
                      "(1e-3/N$^{1.5}$ and 1e-5/N$^{1.5}$; not adopted)")
    if extra_artists:
        extra_artists(ax)
    ax.set_xscale("log")
    ax.set_xlim(0.9, n_max * 1.35)
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.8")
    if legend_loc == "below":
        ax.legend(fontsize=6.5, frameon=True, loc="upper center",
                  bbox_to_anchor=(0.5, -0.13), ncol=2)
    else:
        ax.legend(fontsize=6.5, frameon=True, loc=legend_loc)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGS, f"{outstem}.{ext}"), bbox_inches="tight")
    plt.close(fig)


def mean_metric(expo, needle, y_lo, y_hi):
    vals = [to_float(get_col(r, "value")) for r in expo
            if needle in (get_col(r, "metric") or "").lower()
            and to_int(get_col(r, "year")) is not None
            and y_lo <= to_int(get_col(r, "year")) <= y_hi
            and to_float(get_col(r, "value")) is not None]
    return float(np.mean(vals)) if vals else None


def undercount_factor(recon_rows, y_lo, y_hi):
    scdf = reg = 0
    years_missing_scdf = []
    for r in recon_rows:
        y = to_int(get_col(r, "year"))
        if y is None or not (y_lo <= y <= y_hi):
            continue
        s = to_int(get_col(r, "scdf_total_deaths", "scdf_deaths"))
        g = to_int(get_col(r, "register_deaths"))
        if s is None:
            years_missing_scdf.append(y)
            continue
        scdf += s
        reg += g or 0
    factor = (scdf / reg) if reg else None
    return factor, {"scdf_total_deaths": scdf, "register_deaths": reg,
                    "years_missing_scdf_total": years_missing_scdf}


# ---------------------------------------------------------------------------
# v4 additions
# ---------------------------------------------------------------------------

def tail_stability(events, events_all, T, dwell, criteria=CRITERIA):
    """Quantify how fragile each F(>=N) is to single events - the direct response to the
    adversarial point that "a single event could decide whether the curve crosses a
    criterion line"."""
    sizes = np.asarray([e["fatalities"] for e in events if e["fatalities"] >= 1])
    k = sizes.size
    out = {}

    # (a) Leave-one-out jackknife: removing one event of size >= N lowers the count at N
    # by exactly one; removing any other event leaves it unchanged. The LOO envelope is
    # therefore exact and needs no simulation. Report it per N with the relative swing.
    ns = np.arange(1, sizes.max() + 1)
    loo = []
    for n in ns:
        k_n = int((sizes >= n).sum())
        loo.append({
            "N": int(n), "count": k_n, "F": k_n / T,
            "F_loo_min": (k_n - 1) / T if k_n else 0.0,
            "F_loo_max": k_n / T,
            "relative_swing_pct": round(100.0 / k_n, 1) if k_n else None,
            "events_supporting": k_n,
        })
    out["leave_one_out"] = loo

    # (b) One-more-event perturbations: what a single additional in-window event does.
    add4 = np.append(sizes, 4)
    add5 = np.append(sizes, 5)
    out["one_more_event"] = {
        "add_N4": {"F_ge4": float((add4 >= 4).sum() / T),
                   "note": "one additional four-death event doubles F(>=4)"},
        "add_N5": {"F_ge5": float((add5 >= 5).sum() / T),
                   "note": "one five-death event would extend the curve to N=5 at 1/T"},
        "one_over_T": 1.0 / T,
    }

    # (c) Rolling 10-year windows over the full building-envelope register 2005-2025:
    # the range of F(>=1..3) across every start year. Uses events_all (building scope).
    years = sorted({e["year"] for e in events_all if e["year"]})
    y_lo, y_hi = min(years), max(years)
    win = 10
    rolls = []
    for start in range(y_lo, y_hi - win + 2):
        end = start + win - 1
        s = np.asarray([e["fatalities"] for e in events_all
                        if e["year"] and start <= e["year"] <= end and e["fatalities"] >= 1])
        rolls.append({"window": [start, end],
                      "F_ge1": float((s >= 1).sum() / win),
                      "F_ge2": float((s >= 2).sum() / win),
                      "F_ge3": float((s >= 3).sum() / win)})
    out["rolling_10yr_windows"] = rolls
    for key in ("F_ge1", "F_ge2", "F_ge3"):
        vals = [r[key] for r in rolls]
        out[f"rolling_{key}_range"] = [min(vals), max(vals)]

    # (d) Criterion-crossing decisiveness at per-dwelling scale: how many events of size
    # >= N would the window need to hold for the per-dwelling curve to REACH each line?
    # needed = line(N) * T * mean dwellings. If needed - observed >> 1, no single event
    # (or small number of events) can decide a crossing.
    dec = {}
    if dwell:
        for name, c in criteria.items():
            rows = []
            for n in ns:
                line_val = c["intercept"] * n ** c["slope"]
                needed = line_val * T * dwell
                observed = int((sizes >= n).sum())
                rows.append({"N": int(n), "line_per_dwelling_yr": line_val,
                             "events_needed": needed, "events_observed": observed,
                             "shortfall_events": needed - observed})
            dec[name] = rows
    out["per_dwelling_crossing_decisiveness"] = dec
    # Single-event crossing threshold (round-3 A-M3): one event contributes 1/(T*D) per
    # dwelling-year; a slope -1 line F = a/N falls below that at N > a*T*D, so a single
    # event LARGER than that N would by itself lift the per-dwelling point above the line.
    if dwell:
        out["single_event_crossing_threshold_N"] = {
            name: float(c["intercept"] * T * dwell) for name, c in criteria.items()
            if c["slope"] == -1.0}
    out["note"] = ("LOO envelope is exact (analytic). Per-dwelling decisiveness: number of "
                   "size->=N events the 14-yr window would need for the per-dwelling curve "
                   "to reach each line, vs the number observed. "
                   "single_event_crossing_threshold_N: event size above which ONE event "
                   "crosses the line at its own N (bounds the no-single-event claim).")
    return out


def plot_tail_stability(curves_national, tail, T, outstem="FN_sg_tailstability_v4"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ns, freq, blo, bhi = curves_national
    fig, ax = plt.subplots(figsize=(6.5, 5.0), dpi=300)
    # bootstrap band, floored
    floor = 1.0 / T
    ax.fill_between(ns, np.clip(blo, floor, None), np.clip(bhi, floor, None), step="post",
                    color="0.5", alpha=0.25, label="95% bootstrap band (clipped at 1/T)")
    ax.step(ns, freq, where="post", color="black", lw=1.8,
            label="Observed FN, building fires 2012-2025")
    # exact leave-one-out envelope as whiskers at each N
    loo = tail["leave_one_out"]
    for row in loo:
        n, lo, hi = row["N"], row["F_loo_min"], row["F_loo_max"]
        lo_plot = max(lo, floor * 0.5) if lo > 0 else floor * 0.5
        ax.plot([n, n], [lo_plot, hi],
                color="tab:orange", lw=2.2, solid_capstyle="butt", zorder=5)
        if lo == 0:
            # true lower end is zero, unplottable on a log axis: arrowhead marks the
            # off-scale continuation (round-3 A-m1 / C-m10); caption states the convention
            ax.plot(n, lo_plot, marker="v", ms=5, color="tab:orange", zorder=5, ls="none")
    ax.plot([], [], color="tab:orange", lw=2.2,
            label="Exact leave-one-event-out envelope")
    # terminal-step extension for the observed curve (visible width for the last value);
    # kept short (x1.12) so it stays clear of the hypothetical N=5 marker at (5, 1/T)
    ax.plot([ns[-1], ns[-1] * 1.12], [freq[-1], freq[-1]], color="black", lw=1.8)
    ax.plot(ns[-1] * 1.12, freq[-1], marker="|", ms=7, mew=1.8, color="black", ls="none")
    # one-more-event markers
    ome = tail["one_more_event"]
    ax.plot(4, ome["add_N4"]["F_ge4"], marker="s", mfc="none", mec="tab:red", ms=7, ls="none",
            label="One added 4-death event (F(>=4) doubles)")
    ax.plot(5, ome["add_N5"]["F_ge5"], marker="o", mfc="none", mec="tab:red", ms=7, ls="none",
            label="One hypothetical 5-death event (curve extends)")
    ax.annotate("entire N=4 point rests\non a single event", xy=(4, loo[3]["F"]),
                xytext=(4.6, loo[3]["F"] * 3.2), fontsize=6.5, color="0.2",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.4"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.9, 9)
    ax.set_xlabel("N, fatalities per event")
    ax.set_ylabel("F, annual frequency of events with >= N fatalities (1/yr)")
    ax.set_title("Tail stability of the national FN curve under single-event perturbations",
                 fontsize=10)
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.8")
    ax.legend(fontsize=6.5, frameon=True, loc="lower left")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGS, f"{outstem}.{ext}"), bbox_inches="tight")
    plt.close(fig)


def plot_event_strip(events_all, outstem="FN_sg_eventstrip_v4"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 3.4), dpi=300)
    ax.axvspan(PRIMARY_START_YEAR - 0.5, PRIMARY_END_YEAR + 0.5, color="0.92", zorder=0,
               label="Primary window 2012-2025")
    # jitter events sharing (year, N) horizontally so all dots stay visible
    seen = {}
    for e in sorted(events_all, key=lambda x: (x["year"], -x["fatalities"])):
        if not e["year"] or e["fatalities"] < 1:
            continue
        key = (e["year"], e["fatalities"])
        j = seen.get(key, 0)
        seen[key] = j + 1
        x = e["year"] + (j - (0 if j == 0 else 0.5)) * 0.22 * (1 if j % 2 else -1)
        res = is_residential(e)
        ax.scatter(x, e["fatalities"], s=34 + 26 * (e["fatalities"] - 1),
                   color="tab:blue" if res else "tab:red", alpha=0.85, zorder=3,
                   edgecolors="white", linewidths=0.5)
    ax.scatter([], [], color="tab:blue", label="Residential (occupancy basis)")
    ax.scatter([], [], color="tab:red", label="Non-residential")
    for yr, n, txt, dx, dy in [(2014, 4, "Geylang shophouse-lodging (4)", -6.2, 0.25),
                               (2021, 3, "Tuas dust explosion (3)", -7.8, 0.45),
                               (2022, 3, "Bedok North (3)", -1.6, 0.62),
                               (2025, 3, "Hougang (3)", 0.4, 0.45)]:
        ax.annotate(txt, xy=(yr, n), xytext=(yr + dx, n + dy), fontsize=6.5, color="0.25",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5"))
    ax.set_xlim(2004, 2026.5)
    ax.set_ylim(0.4, 4.9)
    ax.set_yticks([1, 2, 3, 4])
    ax.set_xticks(range(2005, 2026, 5))
    ax.xaxis.set_minor_locator(__import__("matplotlib.ticker", fromlist=["MultipleLocator"]).MultipleLocator(1))
    ax.set_xlabel("Year")
    ax.set_ylabel("N, deaths in event")
    ax.set_title("Every fatal building-envelope fire event, Singapore 2005-2025", fontsize=10)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="0.85")
    ax.legend(fontsize=6.5, frameon=True, loc="upper left", ncol=3)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGS, f"{outstem}.{ext}"), bbox_inches="tight")
    plt.close(fig)


def plot_international(curves_national, T, pop_mean, outstem="FN_sg_international_v4"):
    """Per-million-population international comparison. Draws ONLY if the verified
    data file exists. Every foreign curve in the JSON must carry provenance fields;
    curves without them are refused (F159 discipline: no unverified series enters a
    figure)."""
    path = os.path.join(DATA, INTL_FILE)
    if not os.path.exists(path):
        print("No international_fn_curves.json - skipping international figure")
        return None
    with open(path, encoding="utf-8") as f:
        intl = json.load(f)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ns, freq, _, _ = curves_national
    fig, ax = plt.subplots(figsize=(6.8, 5.2), dpi=300)
    f_sg = freq / (pop_mean / 1e6)
    ax.step(ns, f_sg, where="post", color="black", lw=2.0, zorder=6,
            label=f"Singapore 2012-2025, building fires (this study)")
    ax.plot([ns[-1], ns[-1] * 1.25], [f_sg[-1], f_sg[-1]], color="black", lw=2.0, zorder=6)
    ax.plot(ns[-1] * 1.25, f_sg[-1], marker="|", ms=7, mew=2.0, color="black", ls="none",
            zorder=6)
    palette = ["tab:red", "tab:blue", "tab:green", "tab:purple", "tab:brown", "tab:olive"]
    markers = ["o", "s", "D", "^", "v", "P"]
    used = []
    for i, c in enumerate(intl.get("curves", [])):
        if not (c.get("source") and c.get("pop_millions") and c.get("N") and c.get("F")):
            print(f"REFUSED unverified/incomplete intl curve: {c.get('label')}")
            continue
        f_pm = np.asarray(c["F"], dtype=float) / c["pop_millions"]
        col = palette[i % len(palette)]
        if c.get("style") == "markers":
            # Merged-bin or point-resolved data: markers joined by a faint dotted line,
            # never a step (a step would assert F at unresolved intermediate N).
            ax.plot(c["N"], f_pm, ls=":", lw=0.8, color=col, zorder=3)
            ax.plot(c["N"], f_pm, ls="none", marker=markers[i % len(markers)], ms=4.5,
                    color=col, label=c["label"], zorder=4)
        else:
            ax.step(c["N"], f_pm, where="post", color=col, lw=1.4,
                    ls=c.get("ls", "--"), label=c["label"], zorder=4)
            ax.plot([c["N"][-1], c["N"][-1] * 1.25], [f_pm[-1], f_pm[-1]],
                    color=col, lw=1.4, ls=c.get("ls", "--"), zorder=4)
            ax.plot(c["N"][-1] * 1.25, f_pm[-1], marker="|", ms=6, mew=1.4, color=col,
                    ls="none", zorder=4)
        used.append(c["label"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("N, fatalities per event")
    ax.set_ylabel("F per million population, events with >= N fatalities (1/yr/million)")
    ax.set_title("National fire FN curves, normalised per million population\n"
                 "(scopes differ by source; see note)", fontsize=10)
    ax.text(0.02, 0.02, intl.get("figure_note", ""), transform=ax.transAxes, fontsize=6.0,
            color="0.25", va="bottom",
            bbox=dict(facecolor="white", edgecolor="0.7", lw=0.5, alpha=0.85))
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.8")
    ax.legend(fontsize=6.5, frameon=True, loc="upper right")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGS, f"{outstem}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    return used


def main():
    events_all_full = load_register()  # all building-envelope rows, 2005-2025
    events = [e for e in events_all_full
              if e["year"] and PRIMARY_START_YEAR <= e["year"] <= PRIMARY_END_YEAR]
    T = PRIMARY_END_YEAR - PRIMARY_START_YEAR + 1
    print(f"Window {PRIMARY_START_YEAR}-{PRIMARY_END_YEAR} (T={T}); "
          f"register events in window: {len(events)} of {len(events_all_full)}")

    results = {"window": [PRIMARY_START_YEAR, PRIMARY_END_YEAR], "T_years": T,
               "n_events": len(events), "bootstrap_draws": N_BOOT, "seeds": SEEDS,
               "register_source": REGISTER_FILE,
               "geylang_primary_classification": "residential (occupancy basis, E1)",
               "n_blank_injury_rows": sum(1 for e in events if e["injuries_blank"]),
               "premises_breakdown": {}}
    for e in events:
        results["premises_breakdown"][e["premises"]] = \
            results["premises_breakdown"].get(e["premises"], 0) + 1

    scopes = {
        "national": events,
        "residential": [e for e in events if is_residential(e)],
        "non_residential": [e for e in events if not is_residential(e)
                            and e["premises"] != "UNKNOWN"],
    }
    curves = {}
    for scope, sel in scopes.items():
        sizes = [e["fatalities"] for e in sel]
        ns, freq = fn_curve(sizes, T)
        blo, bhi = bootstrap_band(sizes, T, ns, SEEDS[f"fatalities_{scope}"])
        curves[scope] = (ns, freq, blo, bhi)
        results[f"FN_fatalities_{scope}"] = {
            "N": ns.tolist(), "F": freq.tolist(),
            "ci_lo": blo.tolist() if blo is not None else None,
            "ci_hi": bhi.tolist() if bhi is not None else None,
            "n_fatal_events": int(sum(1 for s in sizes if s >= 1)),
            "total_fatalities": int(sum(sizes)),
        }

    ns, freq, blo, bhi = curves["national"]
    plot_fn([{"ns": ns, "freq": freq, "lo": blo, "hi": bhi,
              "label": "Observed FN, building fires (this study)", "band": True}],
            f"Societal risk of building fires, Singapore "
            f"{PRIMARY_START_YEAR}-{PRIMARY_END_YEAR} (fatalities)",
            "FN_sg_fatalities_v4",
            xlabel="N, fatalities per event",
            ylabel="F, annual frequency of events with >= N fatalities (1/yr)",
            band_floor=1.0 / T)

    r_ns, r_f, r_lo, r_hi = curves["residential"]
    x_ns, x_f, x_lo, x_hi = curves["non_residential"]
    plot_fn([{"ns": r_ns, "freq": r_f, "lo": r_lo, "hi": r_hi,
              "label": "Residential", "color": "tab:blue", "band": False},
             {"ns": x_ns, "freq": x_f, "lo": x_lo, "hi": x_hi,
              "label": "Non-residential", "color": "tab:red", "band": False}],
            f"Residential vs non-residential building fires, "
            f"{PRIMARY_START_YEAR}-{PRIMARY_END_YEAR}",
            "FN_sg_scopes_v4",
            xlabel="N, fatalities per event",
            ylabel="F, annual frequency of events with >= N fatalities (1/yr)")

    expo = read_csv(os.path.join(DATA, EXPO_FILE))
    dwell = mean_metric(expo, "total national residential dwelling stock",
                        PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    pop_mean = mean_metric(expo, "total population",
                           PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    if dwell and r_ns.size:
        plot_fn([{"ns": r_ns, "freq": r_f / dwell,
                  "lo": r_lo / dwell, "hi": r_hi / dwell,
                  "label": "Residential FN per dwelling (this study)", "band": True}],
                f"Per-dwelling societal risk, residential building fires, "
                f"{PRIMARY_START_YEAR}-{PRIMARY_END_YEAR}",
                "FN_sg_perdwelling_v4",
                xlabel="N, fatalities per event",
                ylabel="F, frequency of events with >= N fatalities (per dwelling-year)",
                band_floor=1.0 / (T * dwell), abcb_envelope=True, n_line_max=1000)
        results["FN_perdwelling_residential"] = {
            "N": r_ns.tolist(),
            "F_per_dwelling_yr": [float(x) for x in (r_f / dwell)],
            "mean_dwellings": dwell,
        }

    recon = read_csv(os.path.join(DATA, RECON_FILE))
    factor, detail = undercount_factor(recon, PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    results["reconciliation"] = {"factor_scdf_over_register": factor, **detail}

    nfat_nat = results["FN_fatalities_national"]["n_fatal_events"]
    nfat_res = results["FN_fatalities_residential"]["n_fatal_events"]
    gfa = mean_metric(expo, "national total building gfa",
                      PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    results["normalisations"] = {
        "mean_national_dwellings": dwell,
        "mean_total_population": pop_mean,
        "fatal_event_freq_per_dwelling_yr": (nfat_nat / T / dwell) if dwell else None,
        "residential_fatal_event_freq_per_dwelling_yr": (nfat_res / T / dwell) if dwell else None,
        "n_residential_fatal_events": nfat_res,
        "mean_national_gfa": gfa,
        "fatal_event_freq_per_gfa_yr": (nfat_nat / T / gfa) if gfa else None,
    }

    windows = {"first_half": (PRIMARY_START_YEAR, PRIMARY_START_YEAR + T // 2 - 1),
               "second_half": (PRIMARY_START_YEAR + T // 2, PRIMARY_END_YEAR),
               "extended_register_2005": (2005, PRIMARY_END_YEAR)}
    ws = {}
    for name, (wlo, whi) in windows.items():
        tw = whi - wlo + 1
        pool = events if wlo >= PRIMARY_START_YEAR else [
            e for e in events_all_full if e["year"] and wlo <= e["year"] <= whi]
        sizes = [e["fatalities"] for e in pool if wlo <= e["year"] <= whi]
        ws[name] = {"window": [wlo, whi], "T_years": tw,
                    "n_fatal_events": int(sum(1 for s in sizes if s >= 1)),
                    "F_ge1_per_yr": round(sum(1 for s in sizes if s >= 1) / tw, 4),
                    "F_ge2_per_yr": round(sum(1 for s in sizes if s >= 2) / tw, 4)}
    results["window_sensitivity"] = ws

    results["expected_annual_fatalities"] = \
        results["FN_fatalities_national"]["total_fatalities"] / T

    def fn_points(sel, t=T):
        ns_, f_ = fn_curve([e["fatalities"] for e in sel], t)
        return {"N": ns_.tolist(), "F": [round(x, 4) for x in f_.tolist()],
                "n_fatal_events": len([e for e in sel if e["fatalities"] >= 1])}

    sens = {}
    reclass = [dict(e, premises="COMMERCIAL") if e["id"] == "SGF-2014-04" else e
               for e in events]
    sens["geylang_premises_class_commercial"] = {
        "residential": fn_points([e for e in reclass if is_residential(e)]),
        "non_residential": fn_points([e for e in reclass if not is_residential(e)
                                      and e["premises"] != "UNKNOWN"]),
    }
    excl = [e for e in events if e["id"] != "SGF-2021-05"]
    sens["tuas_process_excluded"] = {
        "national": fn_points(excl),
        "non_residential": fn_points([e for e in excl if not is_residential(e)
                                      and e["premises"] != "UNKNOWN"]),
    }
    delayed = [dict(e, fatalities=4) if e["id"] == "SGF-2022-02" else e for e in events]
    sens["bedok_delayed_death"] = {"national": fn_points(delayed)}
    # (d) Round-3 A-m3: the occupancy rule applied in the REVERSE direction. SGF-2014-01
    # (condominium common-area fire; both victims were on-duty security officers, a
    # working occupancy) moves to non-residential on the same occupancy basis that moves
    # Geylang to residential.
    reclass2 = [dict(e, premises="COMMERCIAL") if e["id"] == "SGF-2014-01" else e
                for e in events]
    sens["condo_commonarea_occupancy_nonres"] = {
        "residential": fn_points([e for e in reclass2 if is_residential(e)]),
        "non_residential": fn_points([e for e in reclass2 if not is_residential(e)
                                      and e["premises"] != "UNKNOWN"]),
    }
    results["classification_sensitivities"] = sens

    # Round-3 B-1: dwelling-stock-matched per-dwelling variant. The SingStat denominator
    # counts dwelling units; two residential-occupancy events occurred in converted
    # commercial lodging premises that are not verifiably inside that stock
    # (SGF-2014-04 Geylang shophouse-lodging, SGF-2015-01 converted workers accommodation).
    # This variant drops them from the numerator; the primary per-dwelling figure keeps
    # them and is relabelled occupancy-matched rather than scope-matched.
    if dwell:
        stock_matched = [e for e in scopes["residential"]
                         if e["id"] not in ("SGF-2014-04", "SGF-2015-01")]
        ns_sm, f_sm = fn_curve([e["fatalities"] for e in stock_matched], T)
        results["FN_perdwelling_dwellingstock_matched"] = {
            "N": ns_sm.tolist(),
            "F_per_dwelling_yr": [float(x) for x in (f_sm / dwell)],
            "n_events": len(stock_matched),
            "excluded_ids": ["SGF-2014-04", "SGF-2015-01"],
        }

    # Round-3 A-M2: resident-population variant of the per-million reading (the primary
    # uses total population = all persons present; SG's "resident" series is a
    # citizens-plus-PRs construct, narrower than the usual-residence concept behind the
    # ABS / Stats NZ ERP comparators).
    respop = mean_metric(expo, "resident population", PRIMARY_START_YEAR, PRIMARY_END_YEAR)
    if respop:
        results["per_million_sg_resident_variant"] = {
            "mean_resident_population": respop,
            "N": ns.tolist(),
            "F_per_million_residents": [float(x) for x in freq / (respop / 1e6)],
        }

    sizes_fi = [e["fatalities"] + e["injuries"] for e in events]
    ns2, f2 = fn_curve(sizes_fi, T)
    b2lo, b2hi = bootstrap_band(sizes_fi, T, ns2, SEEDS["fat_plus_inj_national"])
    results["FN_fat_plus_inj_national"] = {
        "N": ns2.tolist(), "F": f2.tolist(),
        "ci_lo": b2lo.tolist() if b2lo is not None else None,
        "ci_hi": b2hi.tolist() if b2hi is not None else None,
    }

    # ---- v4: tail stability + new figures ----
    tail = tail_stability(events, events_all_full, T, dwell)
    results["tail_stability"] = tail
    plot_tail_stability(curves["national"], tail, T)
    plot_event_strip(events_all_full)
    intl_used = plot_international(curves["national"], T, pop_mean)
    results["international_curves_plotted"] = intl_used

    # Per-million-population values for the manuscript text (same JSON the figure uses;
    # SG values from this run's national curve and the exposure file's population mean).
    intl_path = os.path.join(DATA, INTL_FILE)
    if os.path.exists(intl_path) and pop_mean:
        with open(intl_path, encoding="utf-8") as f:
            intl = json.load(f)
        per_million = {"Singapore (this study)": {
            "N": ns.tolist(),
            "F_per_million": [float(x) for x in freq / (pop_mean / 1e6)],
            "pop_millions": pop_mean / 1e6}}
        for c in intl.get("curves", []):
            if c.get("pop_millions") and c.get("N") and c.get("F"):
                per_million[c["label"]] = {
                    "N": c["N"],
                    "F_per_million": [f / c["pop_millions"] for f in c["F"]],
                    "pop_millions": c["pop_millions"]}
        results["international_per_million"] = per_million

    out = os.path.join(HERE, "fn_sg_results.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(json.dumps({k: v for k, v in results.items()
                      if not isinstance(v, dict) or len(str(v)) < 800}, indent=2))
    print(f"Wrote {out} + figures/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
