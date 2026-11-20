"""Figures of the second-round revision (R2). Writes PNGs to figures_revision/.

  framework.png          Fig. 1   the problem and the three-part construction
  worked_example.png     Fig. 2   estimation on the Singapore register
  matched_fn.png         Fig. 8   Singapore and England per dwelling, and the ratio by number of deaths
  pairs_forest.png       Fig. 9   each dataset against England re-cut to its scope
  two_kinds.png          Fig. 10  victims of single- and multiple-death dwelling fires, England
  grenfell.png           Fig. 11  English dwelling fires and Grenfell Tower
  sg_tail_scenarios.png  Fig. 12  Singapore tail scenarios
  monitoring.png         Fig. 13  monitoring chart
  design_chart.png       Fig. 14  placement chart and estimator
  known_answer.png       Fig. 15  the two forms tested where the answer is known

Run after run_all.py, run_revision.py, england_incident_level.py and run_r2.py:
    python analysis/frequency_severity/r2_figs.py
"""
import json
import os
from decimal import ROUND_HALF_UP, Decimal

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.optimize import brentq
from scipy.stats import chi2

from fs_model import REPO, RESULTS, fit_b, nll, profile_ci, singapore_bins, survival as S

OUT = os.path.join(REPO, "figures_revision")
os.makedirs(OUT, exist_ok=True)
J = lambda name: json.load(open(os.path.join(RESULTS, name), encoding="utf-8"))
R2, RB, EI, FITS, TAIL = J("r2_comparison.json"), J("revision_bases.json"), J("england_incident_level.json"), J("fits_per_jurisdiction.json"), J("singapore_tail.json")
plt.rcParams.update({"font.family": "serif", "font.size": 9.5, "axes.labelsize": 10, "legend.fontsize": 8,
                     "figure.dpi": 300, "axes.linewidth": 0.7})
C_SG, C_EN, C_AU, C_NZ, C_US = "#a63d24", "#2c5a6b", "#c98a3d", "#3f6b45", "#4b7d8f"


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


def pois(k):
    lo = chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
    return lo, chi2.ppf(0.975, 2 * (k + 1)) / 2


def tail_S(b, n, cap=10 ** 4):
    if n > cap:
        return 0.0
    i = np.arange(2, cap + 1, dtype=float); w = i ** (-b)
    return float(w[int(n) - 2:].sum() / w.sum())


# ---------------------------------------------------------------- Fig. 1 framework
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))
ax = axes[0]
Ngrid = np.arange(1, 40)
for b_, c, lab in [(2.3, C_SG, r"$b$ = 2.3 (flat law)"), (4.2, C_EN, r"$b$ = 4.2 (steep law)")]:
    ax.plot(Ngrid, 3.0 * np.array([S(b_, 10 ** 4, n) for n in Ngrid]), "-", color=c, lw=1.4, label=lab)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"Fatalities, $N$"); ax.set_ylabel(r"$F(\geq N)$ [yr$^{-1}$]")
ax.set_title("(a) Same headline rate, different tails", fontsize=9.5)
ax.annotate("identical\n" + r"$F(\geq 1)$", xy=(1, 3), xytext=(1.7, 5.2), fontsize=7.5, arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))
r1 = 3.0 * S(2.3, 10 ** 4, 10); r2 = 3.0 * S(4.2, 10 ** 4, 10)
ax.annotate(f"{r1 / r2:,.0f}" + r"$\times$ apart at $N$ = 10", xy=(10, np.sqrt(r1 * r2)), xytext=(1.15, 2.2e-5), fontsize=7.5,
            arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))
ax.legend(frameon=False, fontsize=7.2, loc="upper right"); ax.set_ylim(1e-6, 30); clean(ax)
ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)


def box(x, y, w, h, text, fc, fs=7.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.10", fc=fc, ec="0.35", lw=0.7))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)


box(0.3, 7.7, 9.4, 1.6, "National F-N curve\n" + r"$F(\geq 1)=\lambda_1+\lambda_2$,   $F(\geq N)=\lambda_2\,S_2(N)$ for $N\geq 2$", "#f2efea", 7.4)
box(0.3, 4.2, 2.85, 2.5, r"$\lambda_1$" + "\nfires with\none death\n\nrate per\nexposure", "#f0e0d9", 6.6)
box(3.55, 4.2, 2.85, 2.5, r"$\lambda_2$" + "\nfires with two\nor more deaths\n\nrate per\nexposure", "#efe6d8", 6.6)
box(6.8, 4.2, 2.85, 2.5, r"$S_2(N)$, $N_{max}$" + "\nshape above\none death,\nand its cap", "#dde8ec", 6.6)
box(0.3, 1.0, 9.4, 1.9, "Regime failure (dragon kings): cladding, egress defeat.\nA separate channel that the ordinary law does not describe.", "#e8e2ee", 6.4)
for xc in (1.72, 4.97, 8.22):
    ax.add_patch(FancyArrowPatch((xc, 7.65), (xc, 6.75), arrowstyle="-|>", mutation_scale=9, color="0.35", lw=0.8))
ax.add_patch(FancyArrowPatch((4.97, 4.1), (4.97, 3.05), arrowstyle="-|>", mutation_scale=9, color="0.35", lw=0.8))
ax.text(5.25, 3.45, "outside the ordinary law", fontsize=6.4, color="0.35")
ax.set_title("(b) The curve in three parts", fontsize=9.5)
fig.tight_layout(); fig.savefig(f"{OUT}/framework.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 2 worked example
sg = singapore_bins("buildings"); n_sg = sum(c for _, _, c in sg)
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
grid = np.linspace(1.0, 4.5, 400)
ll = np.array([-nll(b_, sg, 4) for b_ in grid])
bhat = FITS["Singapore"]["b_obs_support"]; ci = FITS["Singapore"]["ci"]
ax = axes[0]
ax.plot(grid, ll - ll.max(), "-", color=C_SG, lw=1.3)
ax.axhline(-1.921, color="0.5", lw=0.8, ls="--"); ax.axvline(bhat, color="0.35", lw=0.7)
for x in ci:
    ax.axvline(x, color="0.6", lw=0.6, ls=":")
ax.set_xlabel("Severity exponent $b$"); ax.set_ylabel("Profile log-likelihood")
ax.set_title("(a) Estimation: $\\hat{b}$ = %.2f [%.2f, %.2f]" % (bhat, ci[0], ci[1]), fontsize=9.5)
ax.text(3.45, -1.75, "95% cutoff", fontsize=7, color="0.4"); ax.set_ylim(-6, 0.4); clean(ax)
ax = axes[1]
Ns = [1, 2, 3, 4]
ax.plot(Ns, [sum(c for lo, hi, c in sg if hi >= n) / n_sg for n in Ns], "o", ms=6, color="k", mfc="w", mew=1.0, label="observed share")
ax.plot(Ns, [S(bhat, 4, n) for n in Ns], "s-", ms=4.5, color=C_SG, lw=1.1, label="fitted law")
ax.set_yscale("log"); ax.set_xticks(Ns)
ax.set_xlabel("Fatalities, $N$"); ax.set_ylabel("$P($severity$ \\geq N)$")
# the p value is a count out of 2,000 replicates and can end in 5 exactly (0.955); round it half up, as the tables do
ax.set_title("(b) Fit and data, Singapore (fit test p = %s)"
             % Decimal(repr(float(FITS["Singapore"]["gof_p"]))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP), fontsize=9.5)
ax.legend(frameon=False, fontsize=7.5); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/worked_example.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 8 matched comparison
P = R2["singapore_vs_england"]["primary"]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7), gridspec_kw={"width_ratios": [1.05, 1]})
ax = axes[0]
for key, col, mk, lab, off in (("england_dwellings", C_EN, "o", "England, dwelling fires", 1.0),
                               ("singapore_dwelling_stock", C_SG, "s", "Singapore, fires in the dwelling stock", 1.0)):
    d = {int(k): v for k, v in P[key]["distribution"].items() if int(k) <= 6}
    ex = P[key]["exposure"]["dwelling_years_M"]
    ns = sorted(d)
    ks = [sum(v for k, v in d.items() if k >= n) for n in ns]
    F = np.array(ks) / ex
    lo = np.array([pois(k)[0] for k in ks]) / ex; hi = np.array([pois(k)[1] for k in ks]) / ex
    ax.step(ns + [ns[-1] + 0.35], list(F) + [F[-1]], where="post", color=col, lw=1.2, alpha=0.75)
    ax.errorbar(np.array(ns) * (1.0 if key.startswith("eng") else 1.045), F, yerr=[F - lo, hi - F], fmt=mk, ms=5.5, color=col, mec="k", mew=0.4,
                elinewidth=0.8, capsize=2, label=lab)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks([1, 2, 3, 4, 5, 6]); ax.set_xticklabels(list("123456"))
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlim(0.9, 7.2); ax.set_ylim(2e-3, 90)
ax.set_xlabel(r"Fatalities per fire, $N$"); ax.set_ylabel(r"$F(\geq N)$ per million dwelling-years")
ax.set_title("(a) Curves per dwelling", fontsize=9); ax.legend(frameon=False, fontsize=6.8, loc="upper right"); clean(ax)
ax = axes[1]
rows = []
for cl, lab in (("single", "One death"), ("multi", "Two or more"), ("ge3", "Three or more")):
    for basis, blab, mk in (("residential", "residential premises", "s"), ("dwelling_stock", "premises in dwelling stock", "D"), ("household_units", "dwelling units only", "^")):
        if (cl == "single" and basis != "residential") or (cl == "ge3" and basis == "household_units"):
            continue
        rr = P["compare_" + basis]["dwelling_years_M"][cl]["ratio"]
        if cl == "single":
            # the count of single-death fires is the same on all three numerator bases
            assert all(P["compare_" + b_]["dwelling_years_M"]["single"]["ratio"]["ratio"] == rr["ratio"]
                       for b_ in ("residential", "dwelling_stock", "household_units"))
            blab = "all three numerators"
        rows.append((lab + ",\n" + blab, rr["ratio"], rr["ci95"], mk))
ypos = np.arange(len(rows))[::-1]
for (lab, r, ci, mk), y in zip(rows, ypos):
    ax.plot([max(ci[0], 0.09), min(ci[1], 7.5)], [y, y], "-", color=C_SG, lw=1.3)
    ax.plot(r, y, mk, color=C_SG, ms=5.5, mec="k", mew=0.4)
ax.axvline(1.0, color="0.4", lw=0.7, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.09, 7.5); ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_yticks(ypos); ax.set_yticklabels([r_[0] for r_ in rows], fontsize=6.4)
ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Singapore rate / England rate", fontsize=9)
ax.set_title("(b) Ratio by number of deaths", fontsize=9); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/matched_fn.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 9 scope-matched pairs
PR = R2["scope_matched_pairs"]
items = [("Singapore, dwelling stock, all causes\nagainst English dwelling fires, all causes", C_SG,
          P["compare_dwelling_stock"]["dwelling_years_M"]["single"]["ratio"], P["compare_dwelling_stock"]["dwelling_years_M"]["multi"]["ratio"]),
         ("Australia, residential, intent excluded\nagainst English accidental dwelling fires", C_AU, PR["Australia"]["per_household"]["single"], PR["Australia"]["per_household"]["multi"]),
         ("New Zealand, domestic, intent excluded\nagainst English accidental dwelling fires", C_NZ, PR["New Zealand"]["per_household"]["single"], PR["New Zealand"]["per_household"]["multi"])]
us = PR["United States"]["per_household"]
fig, ax = plt.subplots(figsize=(6.9, 3.5))
yy = np.arange(len(items) + 1)[::-1] * 1.0
for (lab, col, s1, s2), y in zip(items, yy):
    for rr, dy, mk, fc in ((s1, 0.16, "o", "w"), (s2, -0.16, "o", col)):
        ax.plot(rr["ci95"], [y + dy, y + dy], "-", color=col, lw=1.4)
        ax.plot(rr["ratio"], y + dy, mk, color=col, ms=6, mec="k" if fc != "w" else col, mew=0.6 if fc == "w" else 0.4, mfc=fc)
y = yy[-1]
ax.plot(us["ratio_single"], y + 0.16, "D", color=C_US, ms=5.5, mfc="w", mew=0.8)
ax.plot(us["england_multi_ratio_ci_only"], [y - 0.16, y - 0.16], ":", color=C_US, lw=1.2)
ax.plot(us["ratio_multi"], y - 0.16, "D", color=C_US, ms=5.5, mec="k", mew=0.4)
labs = [it[0] for it in items] + ["United States, residential buildings, all causes\nagainst English dwelling fires (scaled estimate, no interval)"]
ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=6.9)
ax.axvline(1.0, color="0.4", lw=0.7, ls="--"); ax.set_xscale("log"); ax.set_xlim(0.12, 9)
ax.set_xticks([0.15, 0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.15", "0.25", "0.5", "1", "2", "4", "8"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Rate relative to England on the same scope", fontsize=9)
ax.plot([], [], "o", color="0.3", mfc="w", ms=6, label="fires with one death")
ax.plot([], [], "o", color="0.3", ms=6, label="fires with two or more deaths")
ax.legend(frameon=False, fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 1.08)); ax.set_ylim(-0.6, len(items) + 0.6); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/pairs_forest.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 10 two kinds of fatal dwelling fire
T = R2["two_processes_england_dwellings"]
fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.25), gridspec_kw={"width_ratios": [1.15, 1.0, 1.15]})
ax = axes[0]
bands = T["by_age"]["bands"]; xs = np.arange(len(bands))
s_ = np.array(T["by_age"]["single"], float); m_ = np.array(T["by_age"]["multiple"], float)
ax.step(xs, s_ / s_.sum(), where="mid", color="#9fb7c2", lw=1.6, label="single-death fires")
ax.step(xs, m_ / m_.sum(), where="mid", color=C_EN, lw=1.6, label="multiple-death fires")
ax.set_xticks(xs[::4]); ax.set_xticklabels([b.replace("-", "–") for b in bands[::4]], fontsize=6.4)
ax.set_xlabel("Age of victim"); ax.set_ylabel("Share of victims"); ax.set_title("(a) Age", fontsize=9)
ax.legend(frameon=False, fontsize=6.6, loc="upper center"); ax.set_ylim(0, 0.2); clean(ax)
ax = axes[1]
h = np.arange(24)
ax.step(h, T["by_hour"]["single"], where="mid", color="#9fb7c2", lw=1.6)
ax.step(h, T["by_hour"]["multiple"], where="mid", color=C_EN, lw=1.6)
ax.set_xticks([0, 6, 12, 18, 23]); ax.set_xticklabels(["00", "06", "12", "18", "23"], fontsize=7)
ax.set_xlabel("Hour the fire was reported"); ax.set_title("(b) Time of day", fontsize=9); ax.set_ylim(0, 0.13); clean(ax)
ax = axes[2]
keys = [("in the room of origin when the fire started", "in room of origin\nat the start"), ("elsewhere in the building when the fire started", "elsewhere in the\nbuilding at the start"),
        ("overcome by gas or smoke", "overcome by\ngas or smoke"), ("death by burns alone", "burns alone"), ("ignition by smokers' materials", "smokers'\nmaterials")]
yb = np.arange(len(keys))[::-1]; wd = 0.36
ax.barh(yb + wd / 2, [T["shares"][k]["single"]["share"] for k, _ in keys], wd, color="#9fb7c2", ec="k", lw=0.4)
ax.barh(yb - wd / 2, [T["shares"][k]["multiple"]["share"] for k, _ in keys], wd, color=C_EN, ec="k", lw=0.4)
ax.set_yticks(yb); ax.set_yticklabels([l for _, l in keys], fontsize=6.4); ax.set_xlim(0, 0.8)
ax.set_xlabel("Share of victims"); ax.set_title("(c) Circumstances", fontsize=9); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/two_kinds.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 11 Grenfell
d9 = {int(k): v for k, v in EI["by_location"]["Dwellings"]["distribution"].items()}
T9 = 9.0
g = R2["grenfell_under_shape_above_one_death"]["review"]; law = EI["dwelling_severity_law"]
fig, ax = plt.subplots(figsize=(5.4, 4.0))
ns = [n for n in sorted(d9) if n <= 6]
ax.plot(ns, [sum(v for k, v in d9.items() if k >= n) / T9 for n in ns], "o", ms=5, color=C_EN, mec="k", mew=0.4, label="England, dwelling fires, 2010/11–2018/19")
ax.plot(71, 1 / T9, "*", ms=14, color=C_SG, mec="k", mew=0.5, label="Grenfell Tower ($N$ = 71)", zorder=6)
Ng = np.unique(np.round(np.logspace(0, np.log10(400), 80)).astype(int))
lam = law["fatal_dwelling_fires_per_yr"]
ax.plot(Ng, [lam * S(law["b"], 10 ** 4, n) for n in Ng], "-", color=C_EN, lw=1.1, label="single law from $N$ = 1 ($b$ = %.2f)" % law["b"])
Ng2 = Ng[Ng >= 2]
ax.plot(Ng2, [g["per_year"] * tail_S(g["b2"], n) for n in Ng2], "--", color="#7a4f9a", lw=1.1, label="law above one death ($b_2$ = %.2f)" % g["b2"])
ax.annotate("return period\nabout %s years (single law)\nabout %s years (law above one death)" % (format(int(round(law["grenfell_return_period_yr"], -3)), ","), format(g["return_period_yr"], ",")),
            xy=(71, 1 / T9), xytext=(60, 4.0), fontsize=7.3, ha="center", va="bottom", arrowprops=dict(arrowstyle="->", lw=0.7, color="0.3"))
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.9, 400); ax.set_ylim(1e-6, 500)
ax.set_xlabel("Number of fatalities, $N$"); ax.set_ylabel(r"$F(\geq N)$ [events $\cdot$ yr$^{-1}$]")
ax.legend(frameon=False, loc="lower left", fontsize=7.2); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/grenfell.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 12 Singapore tail scenarios
lam_sg = n_sg / 14.0
ST = R2["singapore_tail_with_english_shape"]
fig, ax = plt.subplots(figsize=(5.6, 4.1))
for nmax, col, ls, lab in [(8, C_NZ, "-", r"single law, $N_{max}$ = 8"), (17, C_SG, "-", r"single law, $N_{max}$ = 17"), (10 ** 4, C_EN, "-", "single law, no cap")]:
    b, _ = fit_b(sg, nmax)
    ng = np.arange(1, min(nmax, 80) + 1)
    ax.plot(ng, [lam_sg * S(b, nmax, n) for n in ng], ls, color=col, lw=1.2, label=lab + ", $b$ = %.2f" % (int(round(b, 3) * 100 + 0.5 + 1e-9) / 100))
ng = np.arange(2, 81); l2 = ST["multiple_death_fires_per_year_buildings"]
ax.fill_between(ng, [l2 * tail_S(ST["english_b2_ci95"][1], n) for n in ng], [l2 * tail_S(ST["english_b2_ci95"][0], n) for n in ng], color="#7a4f9a", alpha=0.13, lw=0)
ax.plot(ng, [l2 * tail_S(ST["english_b2"], n) for n in ng], "--", color="#7a4f9a", lw=1.4,
        label=r"three-part law, English shape above one death, $b_2$ = %.2f" % ST["english_b2"])
obs = [(n, sum(c for lo, hi, c in sg if hi >= n) / 14.0) for n in (1, 2, 3, 4)]
ax.plot([o[0] for o in obs], [o[1] for o in obs], "o", ms=6, color="k", mfc="w", mew=1.1, label="observed 2012–2025", zorder=6)
ax.axhline(0.02, color="0.7", lw=0.6, ls=":"); ax.text(1.05, 0.0225, "50-yr return period", fontsize=7, color="0.5")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.9, 80); ax.set_ylim(1e-4, 8)
ax.set_xlabel(r"Number of fatalities, $N$"); ax.set_ylabel(r"$F(\geq N)$ [events $\cdot$ yr$^{-1}$]")
ax.legend(frameon=False, fontsize=6.8, loc="lower left"); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/sg_tail_scenarios.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 13 monitoring
zq = np.linspace(-6, 6, 2001); wq = np.exp(-zq ** 2 / 2); wq /= wq.sum()       # the grid used in run_all.py


def post(bins_, nmax):
    b, _ = fit_b(bins_, nmax)
    ci_ = profile_ci(bins_, nmax, grid=np.linspace(1.05, 6.5, 1500))
    sd = (np.log(ci_[1]) - np.log(ci_[0])) / (2 * 1.96)
    return np.exp(np.log(b) + sd * zq)


bu = post(sg, 10 ** 4); b8 = post(sg, 8)
s5u = np.array([S(b, 10 ** 4, 5) for b in bu]); s58 = np.array([S(b, 8, 5) for b in b8])
yrs = np.arange(0, 41)
bf_sharp = np.array([1 / np.sum(wq * (1 - s5u) ** (lam_sg * y)) for y in yrs])
bf_8 = np.array([np.sum(wq * (1 - s58) ** (lam_sg * y)) / np.sum(wq * (1 - s5u) ** (lam_sg * y)) for y in yrs])
json.dump(dict(cap8_against_no_cap_per_decade_single_law=round(float(bf_8[10]), 2),
               note="the factors for a cap at the observed maximum are in singapore_tail.json and r2_comparison.json"),
          open(os.path.join(RESULTS, "monitoring_curve.json"), "w"), indent=1)
sd2 = (np.log(ST["english_b2_ci95"][1]) - np.log(ST["english_b2_ci95"][0])) / (2 * 1.96)
s5t = np.array([tail_S(float(b), 5) for b in np.exp(np.log(ST["english_b2"]) + sd2 * zq)])
bf_three = np.array([1 / np.sum(wq * (1 - s5t) ** (ST["multiple_death_fires_per_year_buildings"] * y)) for y in yrs])
fig, ax = plt.subplots(figsize=(5.4, 3.7))
ax.plot(yrs, bf_sharp, "-", color=C_SG, lw=1.4, label="single law: cap at observed maximum against no cap")
ax.plot(yrs, bf_three, "--", color="#7a4f9a", lw=1.4, label="three-part law: cap at observed maximum against no cap")
ax.plot(yrs, bf_8, "-", color=C_NZ, lw=1.4, label="single law: $N_{max}$ = 8 against no cap")
ax.axhline(10, color="0.6", lw=0.7, ls=":"); ax.text(0.6, 11.5, "strong evidence (factor 10)", fontsize=7, color="0.4")
ax.axhline(1, color="0.8", lw=0.6); ax.set_yscale("log")
ax.set_xlabel("Years of continued observation without a fire of five or more deaths")
ax.set_ylabel("Bayes factor for the capped model")
ax.legend(frameon=False, fontsize=6.8, loc="upper left"); ax.set_ylim(0.8, 60)
ax.text(2.0, 0.845, "one fire of five or more deaths falsifies the sharp cap, and its size locates the cap", fontsize=6.6, color="0.35"); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/monitoring.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 14 placement chart and estimator
SH = R2["shape_above_one_death"]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axes[0]
lam_g = np.logspace(-1.0, 0.6, 200); b_g = np.linspace(2.2, 6.6, 220)
L, B = np.meshgrid(lam_g, b_g)
S5 = np.array([tail_S(b_, 5) for b_ in b_g])
cs = ax.contour(L, B, np.log10(L * S5[:, None]), levels=np.log10([1e-3, 3e-3, 1e-2, 3e-2, 1e-1]), colors="0.6", linewidths=0.7)
ax.clabel(cs, fmt={lv: ("%.3f" % 10 ** lv).rstrip("0") for lv in cs.levels}, fontsize=6.3)
# the two points nearly coincide, so England is a filled circle labelled above and Singapore an open square labelled below
for key, col, blk, lab, off, mk in (("England, dwellings", C_EN, P["england_dwellings"], "England,\ndwelling fires", (1.12, 0.10), dict(fmt="o", ms=5)),
                                    ("Singapore, dwelling stock", C_SG, P["singapore_dwelling_stock"], "Singapore,\ndwelling stock", (1.12, -0.62),
                                     dict(fmt="s", ms=10, mfc="none"))):
    ex = blk["exposure"]["dwelling_years_M"]; k2 = blk["multi"]; x = k2 / ex
    xlo, xhi = pois(k2)[0] / ex, pois(k2)[1] / ex
    b2, ci2 = SH[key]["b2"], SH[key]["b2_ci95"]
    ax.errorbar(x, b2, xerr=[[x - xlo], [xhi - x]], yerr=[[b2 - ci2[0]], [min(ci2[1], 6.55) - b2]], color=col, mec=col if "mfc" in mk else "k",
                mew=1.2 if "mfc" in mk else 0.4, elinewidth=0.8, capsize=2.5, **mk)
    ax.annotate(lab, (x * off[0], b2 + off[1]), fontsize=7.0, color=col)
ax.set_xscale("log"); ax.set_xlim(0.1, 4); ax.set_ylim(2.2, 6.6)
ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 4]); ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "4"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel(r"Fires with two or more deaths, $\lambda_2$ [per million dwelling-years]", fontsize=8.3); ax.set_ylabel(r"Exponent above one death, $b_2$")
ax.set_title(r"(a) Placement chart: $F(\geq 5)$ per million dwelling-years", fontsize=8.6); clean(ax)
ax = axes[1]
p_g = np.linspace(0.50, 0.92, 300)


def b2_of_share(p2, cap):
    i = np.arange(2, cap + 1, dtype=float)
    return brentq(lambda b: 2.0 ** (-b) / np.sum(i ** (-b)) - p2, 1.01, 14)


for Nm, c, lab in [(8, C_NZ, r"$N_{max}$ = 8"), (17, C_SG, r"$N_{max}$ = 17"), (10 ** 4, C_EN, "no cap")]:
    ax.plot(p_g, [b2_of_share(p_, Nm) for p_ in p_g], "-", color=c, lw=1.2, label=lab)
for key, col, mk in (("England, dwellings", C_EN, "o"), ("Singapore, dwelling stock", C_SG, "s"), ("Australia_preventable_residential", C_AU, "P"),
                     ("NewZealand_unintentional_domestic", C_NZ, "X")):
    ax.plot(SH[key]["share_exactly_two"], SH[key]["b2"], mk, ms=5.5, color=col, mec="k", mew=0.3)
us_p2 = SH["United States, residential buildings"]["share_exactly_two"]
ax.plot(us_p2, b2_of_share(us_p2, 10 ** 4), "D", ms=5, color=C_US, mfc="w", mew=0.9)
ax.set_xlabel("Share with exactly two deaths\namong fires with two or more", fontsize=7.8); ax.set_ylabel(r"Implied exponent $b_2$")
ax.set_ylim(2.0, 7.0); ax.set_title("(b) Estimator from published counts", fontsize=9); ax.legend(frameon=False, fontsize=7, loc="upper left"); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/design_chart.png"); plt.close(fig)

# ---------------------------------------------------------------- Fig. 15 a test where the answer is known
KT = R2["known_answer_test"]; HHT = KT["households"]; SZ = KT["singapore_sized_english_records"]
C_ONE, C_THREE = "#333333", C_NZ
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5), gridspec_kw={"width_ratios": [1.12, 1]})
ax = axes[0]
cv = HHT["curve_as_single_death_fires_are_removed"]
lo_b, hi_b = HHT["official_bounds"]["5"]
ax.axhspan(lo_b, hi_b, color="0.86", lw=0, label="official counts allow %d to %d" % (lo_b, hi_b))
ax.plot([c["kept"] for c in cv], [c["single_law_5_or_more"] for c in cv], "-", color=C_ONE, lw=1.4, label="single power law")
three = HHT["cases"]["all households"]["three_part"]["5"]
ax.axhline(three, color=C_THREE, lw=1.4, ls="--", label="three-part form")
all_single = HHT["cases"]["all households"]["single_death_fires"]
labels = {"all households": ("all households", (0.93, 1.25), "left"),
          "households of two or more, unstated type counted with them": ("households of two or more\n(unstated type counted\nwith them)", (0.8, 8.7), "left"),
          "households of two or more, unstated type spread in proportion": ("households of two\nor more (unstated\ntype spread)", (0.2, 8.4), "left")}
for name, case in HHT["cases"].items():
    x = case["single_death_fires"] / all_single; y = case["single_law"]["5"]
    ax.plot(x, y, "o" if name == "all households" else "s", ms=6, color=C_ONE, mec="k", mew=0.4, zorder=5)
    txt, (tx, ty), ha = labels[name]
    ax.annotate(txt, (x, y), xytext=(tx, ty), fontsize=6.6, ha=ha, va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
sgx = HHT["at_singapores_single_death_rate"]["kept"]
ax.axvline(sgx, color=C_SG, lw=0.8, ls=":")
ax.text(sgx * 0.97, 0.45, "Singapore's rate of\nsingle-death fires", fontsize=6.4, color=C_SG, ha="left", va="bottom")
ax.set_xscale("log"); ax.set_xlim(1.08, 0.095); ax.set_ylim(0, 13.5)
ax.set_xticks([1, 0.5, 0.25, 0.1]); ax.set_xticklabels(["1", "0.5", "0.25", "0.1"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Single-death fires counted, as a share of the full record", fontsize=8.6)
ax.set_ylabel("Expected fires of five or more deaths", fontsize=9)
ax.set_title("(a) Same multiple-death fires, 2010/11 to 2018/19", fontsize=9)
ax.legend(frameon=False, fontsize=6.8, loc="upper left"); clean(ax)

ax = axes[1]
XMAX = 1500.0
SB = KT["records_drawn_from_a_single_power_law"]
groups = [(SZ, 3.0, "sizes of multiple-death fires\nfrom the English record"), (SB, 1.0, "all sizes from a single\npower law (exponent %.2f)" % SB["exponent"])]
yt, yl = [], []
for blk, y0, title in groups:
    p1 = blk["single_law_return_period_yr_pct_2_5_25_50_75_97_5"]; p3 = blk["three_part_return_period_yr_pct_2_5_25_50_75"]
    truth = blk["true_return_period_5_or_more_yr"]
    ya, yb = y0 + 0.3, y0 - 0.3
    ax.plot([truth, truth], [yb - 0.28, ya + 0.28], "-", color="0.25", lw=1.3)
    ax.text(truth * 1.07, ya + 0.33, "true value %d" % round(truth), fontsize=6.4, color="0.25", ha="left", va="center")
    ax.text(XMAX * 0.95, ya + 0.72, title, fontsize=6.6, color="0.2", ha="right", va="center", style="italic")
    ax.plot([p1[0], min(p1[4], XMAX * 0.9)], [ya, ya], "-", color=C_ONE, lw=1.0); ax.plot([p1[1], p1[3]], [ya, ya], "-", color=C_ONE, lw=3.2, solid_capstyle="butt")
    ax.plot(p1[2], ya, "o", ms=5.5, color="w", mec=C_ONE, mew=1.2, zorder=5)
    ax.plot([p3[0], XMAX * 0.8], [yb, yb], "-", color=C_THREE, lw=1.0); ax.plot([p3[1], p3[3]], [yb, yb], "-", color=C_THREE, lw=3.2, solid_capstyle="butt")
    ax.plot(p3[2], yb, "o", ms=5.5, color="w", mec=C_THREE, mew=1.2, zorder=5)
    ax.annotate("", xy=(XMAX * 0.97, yb), xytext=(XMAX * 0.7, yb), arrowprops=dict(arrowstyle="->", lw=1.0, color=C_THREE))
    ax.text(XMAX * 0.95, yb - 0.27, "no finite estimate in %d%%" % round(100 * blk["three_part_share_of_records_with_no_finite_estimate"]),
            fontsize=6.0, color=C_THREE, ha="right", va="center")
    yt += [ya, yb]; yl += ["single law", "three-part"]
ax.set_xscale("log"); ax.set_xlim(4, XMAX); ax.set_ylim(0.05, 4.35)
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=7.6)
ax.set_xticks([10, 30, 100, 300, 1000]); ax.set_xticklabels(["10", "30", "100", "300", "1,000"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Estimated return period of a fire of\nfive or more deaths [years]", fontsize=8.6)
ax.set_title("(b) Simulated records of 42 fatal fires", fontsize=9); clean(ax)
fig.tight_layout(); fig.savefig(f"{OUT}/known_answer.png"); plt.close(fig)

print("figures written to", OUT)
