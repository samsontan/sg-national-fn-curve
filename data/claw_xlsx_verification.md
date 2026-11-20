# Verification of three unverified Claw compilations

Source workbooks used during data assembly:
- `Claw_Singapore_Fire_Occurrences_2012-2024.xlsx` (Sheet1, A1:O24)
- `Claw_Singapore_GFA_2012-2023.xlsx` (Singapore GFA, A1:N32)
- `Claw_Singapore_Ignition_Frequency_2012-2023.xlsx` (Sheet2 A1:O61, Sheet3 A1:E34)

None of the three files carry embedded source URLs, access dates, or methodology notes — every cell had to be re-traced externally. Primary sources used for re-tracing:
- SCDF/SingStat dataset **"Fire Occurrences, Annual"**, data.gov.sg `d_808473a208220960f07a0b064ef16bde` — full 1968-2025 series pulled via the data.gov.sg download API, accessed 2026-07-11.
- SCDF/SingStat dataset **"Fire Injuries And Fire Fatalities, Annual"**, data.gov.sg `d_2c81b575edc555f6c8f0cb7e09c8df02` — full series pulled the same way, accessed 2026-07-11.
- SCDF **"Emergency Medical Services, Fire and Enforcement Statistics 2024"** news release PDF (scdf.gov.sg), accessed 2026-07-11.
- SCDF **"Emergency Medical Services, Fire and Enforcement Statistics 2025"** news release PDF (scdf.gov.sg), accessed 2026-07-11.
- SingStat/data.gov.sg dataset **"Residential Dwellings, Annual"**, `d_6c58a523a05b55836f383fa2a68d332d` — full series, accessed 2026-07-11.

---

## File 1: Claw_Singapore_Fire_Occurrences_2012-2024.xlsx

### Block A — "Fire occurrences in Singapore, by property type" (rows 1-5): Residential / Non-Residential / Non-building structures / Total, 2012-2024

**Verdict: VERIFIED — full match, all 13 years, all 4 rows.**

Cross-checked cell-by-cell against `d_808473a208220960f07a0b064ef16bde` (raw CSV pulled directly, not just spot-checked). Every one of the 52 values (13 years x 4 rows) matches the primary dataset exactly, e.g.:
- 2012: Residential 3184, Non-Res 582, Non-Building 719, Total 4485 — exact match.
- 2018: Residential 2411, Non-Res 515, Non-Building 959, Total 3885 — exact match.
- 2024: Residential 968, Non-Res 415, Non-Building 607, Total 1990 — exact match (also independently confirmed against the SCDF 2024 news-release PDF, para 2 and para 10: "residential fires... decreased by 0.2%... to 968 calls").

This block also matched the SCDF 2025 news release for the (out-of-range but useful) 2025 figures: 1051/471/528/2050.

**Critical finding — undisclosed category break in 2019, not flagged anywhere in the file:** the source dataset's own methodology note states "From 2019, data exclude fire incidents involving minor fires such as fires involving rubbish bins and rubbish chutes... data on fire occurrences before and from 2019 are not directly comparable." The steep drop visible in the file (Residential 2411 in 2018 -> 1168 in 2019) is this break, not a real 49% one-year fall in fires. The xlsx presents both segments as one unbroken series with no annotation. See `frequency_provenance.md` for detail.

### Block B — "Residential fires" cause breakdown (rows 8-11): Rubbish Chute/Bin, Unattended cooking fires, Discarded items fires, Electrical, 2012-2024

**Verdict: MIXED — Unattended cooking / Electrical VERIFIED for years with a primary match; Discarded items fires DISCREPANT (fabricated); Rubbish Chute/Bin UNTRACEABLE for 2019-2021.**

- **Unattended cooking fires**: 2023 = 384, 2024 = 335 in the file. SCDF 2024 news release para 11: "fires due to unattended cooking decreased by 49 calls (-12.8%) to 335 calls" in 2024, implying 2023 = 384. Exact match both years. VERIFIED for the years checked.
- **Electrical**: 2023 = 276, 2024 = 299 in the file. SCDF 2024 news release para 12: "electrical origin increased by 23 calls (8.3%) to 299 calls" in 2024, implying 2023 = 276. Exact match both years. VERIFIED for the years checked.
- **Discarded items fires**: file values are 707, 595, 535, 545, 550, 551, 552, 553, 554, 555, 99, 91, (blank) for 2012-2024. The run 545-546-...-555 (2015-2021) increases by exactly 1 every year for seven straight years — not a pattern any real fire-incident count produces. This is an interpolation/fabrication artifact, not sourced data. **DISCREPANT.** Separately, SCDF's own reporting does not publish a "discarded items" cause specifically for *residential* fires in recent years (its residential Chart 2/3 only breaks out cooking and electrical as the "top two causes"; "discarded items" appears only as a *non-building-premises* category — 97 calls in 2024, 107 in 2025 per the SCDF news releases — a structurally different bucket than "residential discarded items"). The file's 2024 cell for this row is blank, and the closest number for a same-named 2024 category (97, non-building) does not match any nearby file value, confirming the residential-specific series has no traceable primary source.
- **Rubbish Chute/Bin**: values run 1355, 1289, 1379, 1546, 1444, 1240, 1161, 1168, 1175, 1182, then blank for 2022-2024. Given the confirmed 2019 methodology break (rubbish bin/chute fires excluded from official totals from 2019), a continuing named "Rubbish Chute/Bin" count for 2019-2021 (1168, 1175, 1182 — itself a suspiciously smooth +7/+7 progression) could not be traced to any published SCDF table in the time available. **UNTRACEABLE for 2019-2021**; 2012-2018 plausible in magnitude but not independently confirmed against a machine-readable primary source this session.

### Block C — "Non-residential fires" and "Commercial / Industrial / Social-communal" breakdown (rows 14-17)

**Verdict: Non-residential Total row VERIFIED (matches Block A); Commercial/Industrial/Social-communal sub-rows UNTRACEABLE** — no data.gov.sg dataset or SCDF news-release table with this exact 3-way non-residential cause split was located this session. Flagged for follow-up, not confirmed fabricated (no internal red flags like Block B's arithmetic-sequence artifact were found here).

### Block D — "Personal Mobility Devices (PMD) fires" (row 20)

**Verdict: DISCREPANT for 2012-2019 (fabricated); plausible/partially VERIFIED for 2020-2024.**

File values: 31, 49, 74, 115, 31, 49, 74, 115, 68, 55, 42, 55, 67 for 2012-2024. The first four values (2012-2015: 31, 49, 74, 115) are an **exact repeat** of the next four (2016-2019: 31, 49, 74, 115) — a duplicated block, not independent yearly data. This is also substantively implausible: registrable PMDs (e-scooters) did not exist as a distinct Singapore regulatory/statistical category before LTA's Active Mobility Act registration scheme rolled out through 2019; SCDF's AMD/PMD fire-count reporting only becomes a named annual statistic from around 2020-2021 in its news releases. A "PMD fires" count of 31 for 2012 has no plausible source.
- 2023 = 55, 2024 = 67 in the file: matches the SCDF 2024 news release exactly (para 21: "Fires involving AMDs increased by 12 calls (21.8%) to 67 calls in 2024... compared to 55 calls" in 2023). Also note the file's category label "PMD" is stale — SCDF itself renamed/broadened this to "AMD" (Active Mobility Devices = PABs + PMDs + PMAs) from around 2023; the file's row is really tracking the AMD total, not PMD-only (PMD-only 2024 = 25, PMD-only 2023 = 18 per the same PDF). This is a category-label mismatch, not just a numeric one.
- 2020-2022 (68, 55, 42) could not be independently confirmed this session (no primary source located with that exact 3-year run); treat as unverified pending a check against SCDF's 2021/2022/2023 annual statistics releases.

### Block E — "STATISTICS FROM SCDF": Fire Injuries, Fire Fatalities, Fatality Rate (rows 22-24), 2012-2024

**Verdict: VERIFIED for 12 of 13 years; DISCREPANT for 2019 specifically, all three rows.**

Cross-checked cell-by-cell against `d_2c81b575edc555f6c8f0cb7e09c8df02` (raw CSV, full match not spot check):
- Fire Injuries: file vs primary — 2012:199/199, 2013:162/162, 2014:111/111, 2015:111/111, 2016:62/62, 2017:60/60, 2018:90/90, **2019: 141/142 (MISMATCH)**, 2020:184/184, 2021:193/193, 2022:171/171, 2023:81/81, 2024:80/80.
- Fire Fatalities: file vs primary — 2012:1/1, 2013:4/4, 2014:8/8, 2015:7/7, 2016:2/2, 2017:3/3, 2018:4/4, **2019: 2/1 (MISMATCH)**, 2020:1/1, 2021:4/4, 2022:8/8, 2023:3/3, 2024:5/5.
- Fatality Rate (per 100,000 population): file vs primary — matches every year except **2019: 0.04/0.02 (MISMATCH, consistent with the fatalities miscount above)**.

All other 12 years across all three rows are exact matches, including a two-decimal-place match on the rate — this row is genuinely well-sourced apart from the single 2019 transcription error.

---

## File 2: Claw_Singapore_GFA_2012-2023.xlsx

### "NUMBER OF DWELLING UNITS" block (rows 2-10)

**Verdict: VERIFIED — full match, 2012-2022, all rows.** Cross-checked cell-by-cell against `d_6c58a523a05b55836f383fa2a68d332d` (raw CSV): Total Residential Dwelling Units, Total HDB Flats, HDB 1-2/3/4/5-room, Condominiums and Other Apartments, Landed Properties, Others — every cell 2012-2022 matches exactly (e.g. 2012 Total = 1,216,217; HDB 1-2 Room = 61,234; Condo = 204,340 — all exact). The file's 2023 column (Total = 1,545,800) does **not** match the primary dataset's 2023 figure of 1,563,641 — see below, this is because 2023 in the file is not primary data at all but a growth-rate extrapolation (see Sheet3 of File 3).

### "RESIDENTIAL BUILDING FLOOR AREAS" block (rows 13-21)

**Verdict: DISCREPANT — not primary-sourced data; internally-generated proxy.**

Every 2012-2022 floor-area cell in this block is exactly reproduced by multiplying the corresponding dwelling-unit count (from the block above) by a fixed constant sitting in column B of the same row: HDB 1-2 Room x 45, HDB 3-Room x 70, HDB 4-Room x 90, HDB 5-Room x 110, Condominiums x 93.5, Landed x 250, Others x 300 (all sqm/unit). Verified programmatically: `unit_count x constant = floor_area_cell` to the exact integer for every one of the 77 cells checked (7 rows x 11 years, 2012-2022). No citation for these seven per-unit-type average-size assumptions exists anywhere in the file, and they were not traced to any BCA/URA/SingStat published average-unit-size table this session. **This block is not measured GFA — it is a derived proxy computed from the (verified) unit-count table using unsourced constants.** The 2023 column breaks the formula entirely (e.g. HDB 1-2 Room: 103,600 units x 45 = 4,662,000, but the file shows 4,859,004) — 2023 instead matches a straight-line "growth rate" extrapolation applied to the 2022 value (see File 3, Sheet3, which contains this exact working: "HDB 1-2... Growth Rate (%) 3.51... 2023 (Estimated) 4859004").

### "NON-RESIDENTIAL BUILDING FLOOR AREAS" block (rows 24-32)

**Verdict: UNTRACEABLE.** Hotel / Offices / Retail / Business Park / Multi-user factory / Single-user factory / Warehouse floor areas, 2012-2023. No internal-formula artifact was found (unlike the residential block), so these may be genuinely sourced, but no matching published series (SingStat Table Builder stock-of-commercial-and-industrial-space tables, BCA Annual Report, or URA REALIS) was located and confirmed within this session's search effort. Flagged for a dedicated follow-up pass against SingStat Table Builder (tables under the "Building and Construction" / "Real Estate" themes) or BCA's annual "Building" stock statistics before this block is used in the manuscript.

---

## File 3: Claw_Singapore_Ignition_Frequency_2012-2023.xlsx

### Sheet2: mirrors Files 1 and 2's tables plus a computed "IGNITION FREQUENCY [fire/m2 per year]" row for Residential and Non-Residential, 2012-2021 only (2022-2023 columns blank despite the sheet title claiming 2012-2021 and the header row listing through 2023)

**Verdict: Arithmetic internally consistent; DISCREPANT by inheritance (denominator problem).**

Spot-checked the ignition-frequency arithmetic itself: 2012 Residential = 3184 fires / 122,360,610 m² = 2.6021e-05, matching the cell's 2.60214459538899e-05 to full float precision. The division is done correctly. But because the floor-area denominator (File 2's residential GFA block) is a derived proxy with unsourced per-unit-size constants rather than measured floor area, and the numerator inherits File 1's undisclosed 2019 methodology break plus the Block B fabrication issues noted above, **the ignition-frequency series itself cannot be called verified** — it is correct arithmetic performed on an unverified/partly-fabricated input pair. It should not be used as-is; it needs to be rebuilt once (a) primary non-residential/residential GFA is sourced independently and (b) the fire-count numerator is rebuilt with the 2019 break explicitly segmented.

### Sheet3: "Category / 2021 / 2022 / Growth Rate (%) / 2023 (Estimated)" scratch table

**Verdict: This is the smoking gun for the whole GFA/dwelling-unit 2023 column** — it is explicitly labeled "Estimated" and shows the exact growth-rate formula (e.g. HDB 1-2 Room: 2021=4,536,540, 2022=4,695,750, growth=3.51%, 2023 (Est.)=4,859,004) that reproduces File 2's un-flagged "2023" column verbatim. File 2 presents this 2023 column as if it were on the same footing as 2012-2022 primary data; it is not — it is a same-workbook internal estimate carried forward without a label.

---

## Summary table

| Series | Verdict |
|---|---|
| Fire occurrences by premises type (Res/Non-Res/Non-Building/Total), 2012-2024 | **VERIFIED** (full match, all years) — but 2019 category break not disclosed in file |
| Residential cause: Unattended cooking | VERIFIED (2023-2024 confirmed; earlier years unconfirmed this session) |
| Residential cause: Electrical | VERIFIED (2023-2024 confirmed; earlier years unconfirmed this session) |
| Residential cause: Discarded items | **DISCREPANT** — fabricated (arithmetic-sequence artifact 2015-2021) |
| Residential cause: Rubbish Chute/Bin | UNTRACEABLE for 2019-2021; unconfirmed 2012-2018 |
| Non-residential Commercial/Industrial/Social-communal split | UNTRACEABLE |
| PMD/AMD fires | **DISCREPANT** for 2012-2019 (duplicated block, implausible); VERIFIED for 2023-2024; unconfirmed 2020-2022; category label stale (PMD -> AMD) |
| Fire Injuries / Fatalities / Fatality Rate | VERIFIED 12/13 years; **DISCREPANT for 2019** (injuries 141 vs 142, fatalities 2 vs 1, rate 0.04 vs 0.02) |
| Dwelling units by type, 2012-2022 | **VERIFIED** (full match) |
| Dwelling units, 2023 column | DISCREPANT — is a growth-rate estimate, not primary 2023 data (real 2023 primary total = 1,563,641 vs file's 1,545,800) |
| Residential floor areas | **DISCREPANT** — derived proxy (unit count x unsourced constant), not measured GFA |
| Non-residential floor areas | UNTRACEABLE |
| Ignition frequency (fire/m²/yr) | Arithmetic correct; **DISCREPANT by inheritance** (built on unverified GFA denominator + un-segmented fire-count numerator) |

**Bottom line: none of the three files should be used as-is in the manuscript.** The premises-type fire counts and the dwelling-unit counts are solid (now independently re-verified against primary machine-readable data, full series 2012-2025 pulled fresh into `frequency_series_draft.csv` / `exposure_national_draft.csv`). Everything downstream of those two tables in the existing files — cause breakdowns, PMD/AMD row, GFA, and ignition frequency — needs to be rebuilt from primary sources rather than reused.
