# Frequency series provenance — national SCDF fire occurrences, 2012-2025

## Primary sources used

1. **"Fire Occurrences, Annual"** — data.gov.sg dataset `d_808473a208220960f07a0b064ef16bde`, source agency Singapore Civil Defence Force (SCDF), maintained/published by SingStat. Full 1968-2025 time series pulled via the data.gov.sg public download API (`https://api-open.data.gov.sg/v1/public/api/datasets/d_808473a208220960f07a0b064ef16bde/poll-download`) and downloaded as CSV directly (not scraped from a rendered page). Last updated by the agency 2 Mar 2026 per the dataset page. Accessed 2026-07-11.
2. **SCDF news release "Emergency Medical Services, Fire and Enforcement Statistics 2024"** — `https://www.scdf.gov.sg/docs/default-source/media-room-(publications)/annual-statistics/fire-ems-and-enforcement-statistics-2024.pdf`. Accessed 2026-07-11 (used for narrative cross-checks and to source residential-cause counts).
3. **SCDF news release "Emergency Medical Services, Fire and Enforcement Statistics 2025"** — `https://www.scdf.gov.sg/docs/default-source/media-room-(publications)/annual-statistics/scdf-annual-statistics-on-ems-fire-enforcement-2025.pdf`. Accessed 2026-07-11 (extends coverage to 2025 and gives AMD/PMD sub-breakdown).
4. **"Fire Injuries And Fire Fatalities, Annual"** — data.gov.sg dataset `d_2c81b575edc555f6c8f0cb7e09c8df02`, same agency/publisher pattern as (1), same download method. Accessed 2026-07-11.

## Category structure and the confirmed 2019 break

SCDF/SingStat publish fire occurrences in four categories in dataset (1):
- `Total Number Of Fire Occurrences`
- `Fire Occurrences In Buildings` (subtotal of the next two)
  - `Fire Occurrences In Residential Buildings`
  - `Fire Occurrences In Non-Residential Buildings`
- `Fire Occurrences In Non-Building Structures`

**The dataset's own methodology note states, verbatim (as extracted from the dataset page):** "From 2019, data exclude fire incidents involving minor fires such as fires involving rubbish bins and rubbish chutes which are of a very low risk relative to other fires... data on fire occurrences before and from 2019 are not directly comparable."

This is a genuine, confirmed, agency-documented category/definition break, not an artifact of this compilation. It is visible directly in the numbers: Residential fires fall from 2,411 (2018) to 1,168 (2019), a 51.5% one-year drop that is definitional, not behavioural — rubbish-chute and rubbish-bin fires (previously the single largest residential sub-cause, ~1,150-1,550 incidents/year through 2012-2018 per the "Rubbish Chute/Bin" row in `Claw_Singapore_Fire_Occurrences_2012-2024.xlsx`) simply stop being counted in the official total from 2019 onward.

**Implication for the FN-curve: the 2012-2018 segment and the 2019-2025 segment of the "Total" and "Residential" series are two different measurement definitions and must not be pooled or fitted as a single continuous frequency trend.** Any Barrois-style or trend-fit treatment across the full 2012-2025 window needs either (a) two separate segments with the break called out, or (b) a reconstructed pre-2019 series with rubbish-bin/chute fires backed out (possible using the "Rubbish Chute/Bin" row values for 2012-2018, themselves independently confirmed as plausible in magnitude though not machine-source-verified this session — see `claw_xlsx_verification.md`), or (c) restricting the manuscript's headline frequency analysis to the 2019-2025 window only, which is internally consistent.

Non-Residential and Non-Building Structures categories are **not** affected by this specific break (rubbish bin/chute fires are a residential-premises phenomenon); their series can be read as continuous 2012-2025, subject to no other documented redefinition being found this session (none was; only the residential-fires rubbish-bin/chute exclusion is documented in the source methodology note).

No other category rename or reclassification was identified in the 2012-2025 window from the source's own metadata. The only other definitional wrinkle found is a labeling one, not a break in the headline series: SCDF's supplementary "Personal Mobility Device (PMD)" fire count was broadened/renamed to "Active Mobility Device (AMD)" (PABs + PMDs + PMAs combined) at some point between 2022 and 2024 in SCDF's own news releases; this affects only the PMD/AMD sub-series (not reproduced in `frequency_series_draft.csv`, which sticks to the four headline categories that are cleanly sourced), not the headline residential/non-residential/non-building/total counts.

## Year coverage achieved

`frequency_series_draft.csv` covers **2012-2025 inclusive (14 years)** for all four categories plus the buildings subtotal, i.e. full coverage of the requested 2012-2025 window including the 2025 data point (published by SCDF Feb 2026, captured here). This is a complete pull of the primary dataset's own published range for these years — not a subset or spot check.

## Coverage limits / things not resolved this session

- Cause-level breakdown of residential fires (cooking / electrical / discarded items / rubbish chute-bin) is only machine-source-confirmed for 2023-2024 (via the two SCDF news-release PDFs). A full 2012-2025 cause-level series was not reconstructed from a single authoritative machine-readable dataset — data.gov.sg does not appear to publish one; SCDF's cause breakdowns before ~2023 exist mainly as chart percentages in annual reports, not tabulated counts, and would need PDF-by-PDF extraction across ~13 annual statistics releases to rebuild properly. Flagged as a follow-up task, not attempted here given time budget.
- Non-residential fire cause breakdown (Commercial/Industrial/Social-communal split) was not traced to a primary source this session.
- PMD/AMD fires were not rebuilt as a clean series here (see `claw_xlsx_verification.md` for why the existing compilation's version is partly fabricated); a rebuild would need each annual SCDF statistics PDF from ~2019/2020 onward, since AMD/PMD reporting only started being a stable annual line item in SCDF's news releases from around then.
