# sg-national-fn-curve

Fatal-fire incident register, reconciliation, exposure series and analysis code for:

> **Societal Fire Risk in Singapore, 2012–2025: A National F–N Curve for Building Fires**
> Samson Tan, Teoh Teik Toe, Paul Joseph, Khalid A. M. Moinuddin
> Submitted to *Fire Safety Journal*, August 2026. Citation to be updated on acceptance.

Singapore publishes fire fatalities as annual aggregates, so the distribution of deaths
across incidents is not observable from the official series alone. This repository holds a
per-incident register assembled from public sources, its year-by-year reconciliation against
Singapore Civil Defence Force (SCDF) published totals for 2005–2025, the exposure
denominators used for normalisation, and the script that reproduces every number and figure
in the paper.

## Contents

| Path | What |
|---|---|
| `data/fatal_incident_register.csv` | The register. 58 incidents, 2005–2025, one row per fatal fire, each with a public source URL and access date |
| `data/scdf_fatality_reconciliation.csv` | Year-by-year register-to-SCDF reconciliation, including the years where the register exceeds the official total |
| `data/exposure_national.csv` | Dwelling-stock and population denominators |
| `data/frequency_series.csv` | SCDF fire-occurrence series used for frequency context |
| `data/international_fn_curves.json` | Digitised comparison data: China (Wang, Lu and Li), Australia, New Zealand |
| `data/register_provenance.md` | How each incident was located and verified, including the two-pass search that recovered SCDF's own itemised releases |
| `data/frequency_provenance.md` | Provenance of the occurrence series, including the 2019 definition break |
| `data/criterion_lines.md` | Formulae and anchors for the HSE, Dutch and Hong Kong reference lines |
| `analysis/fn_sg_national.py` | The analysis. Builds the F–N curves, bootstrap bands, leave-one-out envelope, normalisations and international comparison |
| `analysis/fn_sg_results.json` | Committed output of that script, so the reported numbers can be checked without rerunning |

## Headline numbers

For 2012–2025, building-envelope scope: 47 fatal events, 64 deaths, 41 of them residential.

| N | F(N) yr⁻¹ | 95% bootstrap band |
|---|---|---|
| ≥ 1 | 3.4 | 2.4–4.4 |
| ≥ 2 | 0.9 | 0.4–1.4 |
| ≥ 3 | 0.3 | 0.07–0.6 |
| ≥ 4 | 0.07 | 0–0.2 |

The curve terminates at N = 4. Area under the curve is 4.6 expected deaths per year.
Removing one event changes F(N ≥ 1) by 2% and F(N ≥ 4) by 100%: the first two steps are
stable to any single event, the terminal step is wholly event-dependent.

## Sources and re-identification

Every row carries its source URL. Sources are SCDF annual statistics releases, SCDF and
SingStat open datasets, Ministry of Home Affairs records, and press reports. Nothing here
is confidential or was obtained under restriction.

The `description` field carries the short factual account published at source. Fatal fires
are reported events, and the manuscript's own Table 5 prints area names and cause
categories, so the register is released as assembled rather than redacted. Users should
note that following the source URLs can still identify individuals, and should not do so.

## Reproducing

```
python analysis/fn_sg_national.py
```

Writes results to `fn_sg_results.json`. Compare against the committed copy in `analysis/`.

## Licence

MIT, see `LICENSE`. If you use the register, please cite the paper.

## Reproduction status

Verified 25 August 2026. Running `analysis/fn_sg_national.py` against the data in this
repository reproduces:

- every value in `analysis/fn_sg_results.json` (520 fields compared, zero differences)
- all six figures **byte-for-byte identical** to the images embedded in the submitted
  manuscript (SHA-1 match on each of `figures/*.png`)

```
figures/FN_sg_fatalities_v4.png      3021400214...  = manuscript Figure 1
figures/FN_sg_scopes_v4.png          bb2d071307...  = manuscript Figure 2
figures/FN_sg_perdwelling_v4.png     23fb767651...  = manuscript Figure 3
figures/FN_sg_eventstrip_v4.png      e5bd8f204a...  = manuscript Figure 4
figures/FN_sg_tailstability_v4.png   e2c11efaea...  = manuscript Figure 5
figures/FN_sg_international_v4.png   96a33dbf5f...  = manuscript Figure 6
```
