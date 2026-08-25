# Register Provenance: Fatal Building Fires in Singapore, 2005-2025

**Compiled:** 2026-07-11. All source URLs accessed 2026-07-11.
**Companion files:** `fatal_incident_register_draft.csv` (58 incidents, 78 deaths), `scdf_fatality_reconciliation.csv` (2001-2025).
**Scope:** every fatal structural/building fire in Singapore findable in public record. Vehicle-only, vegetation, and outdoor-only device fires excluded. PMD/battery fires inside dwellings or building common areas included (cause recorded). Two borderline inclusions are flagged in the `building_scope` column (2007 Pulau Ayer Chawan process plant; 2013 Jalan Samulun shipyard/vessel fire): both are inside SCDF's own fire-fatality series but are not conventional building envelopes, so the FN analysis should test sensitivity to their exclusion.

## 1. Search strategy actually used

1. **Official totals first.** The reconciliation target series came from the data.gov.sg dataset "Fire Injuries and Fire Fatalities, Annual" (SingStat, original source SCDF, dataset id `d_2c81b575edc555f6c8f0cb7e09c8df02`), read via the dataset viewer and the raw datastore API. Cross-checked against: MHA written Parliamentary reply of 18 Feb 2022 (2017-2021 table), SCDF news release of 11 Feb 2026 (2025 statistics), and year-specific press coverage.
2. **Per-incident register, two passes.**
   - Pass A (web news archaeology): six parallel year-band searches (2012-14, 2015-17, 2018-20, 2021-23, 2024-25, plus the totals sweep) across Straits Times, CNA, TODAY, Mothership, Yahoo SG, The Online Citizen, The Star, AsiaOne, coroner-inquiry coverage, MOM/WSH reports, and SCDF/MHA newsroom pages.
   - Pass B (primary documents): direct download and local text extraction (pdfplumber) of SCDF's annual fire statistics PDFs for every year 2005-2025 from the SCDF Annual Statistics publications index. This pass proved decisive: for most years SCDF's own release itemizes every fatal incident by location and date, something web fetch tools missed because the PDFs defeat remote text extraction.
3. **Reconciliation.** Per-year register sums checked against the official series; every mismatch investigated rather than smoothed over (see section 3).

## 2. Sources consulted (primary, load-bearing)

- SCDF annual fire statistics PDFs, 2005-2025 (all downloaded and text-extracted locally). Index: https://www.scdf.gov.sg/home/about-scdf/media-room/publications/annual-statistics
- data.gov.sg "Fire Injuries and Fire Fatalities, Annual": https://data.gov.sg/datasets/d_2c81b575edc555f6c8f0cb7e09c8df02/view (covers 2001-2025; last updated 2 Mar 2026 at access)
- MHA written reply to PQ on fire fatalities (18 Feb 2022, 2017-2021 table): https://www.mha.gov.sg/media-room/newsroom/written-reply-to-pq-on-fire-fatalities-over-the-last-five-years-and-extension-of-the-home-fire-alarm-device-assistance-scheme/
- MHA response to Adjournment Motion on Fire Safety in Our Homes (2025 total of six, named incidents): https://www.mha.gov.sg/media-room/newsroom/response-to-adjournment-motion-on-fire-safety-in-our-homes/
- SCDF news release, EMS/Fire/Enforcement Statistics 2025 (11 Feb 2026): https://www.scdf.gov.sg/home/about-scdf/media-room/latest-happenings/newsarticledetail/emergency-medical-services--fire-and-enforcement-statistics-2025
- Press corroboration per incident: Mothership, The Online Citizen, Yahoo News SG, The Star, Malay Mail, theindependent.sg, STOMP, Must Share News, TWC2, Dust Safety Science, bse.com.sg (URLs per row in the register CSV).

## 3. Key findings about the official series (headline material)

1. **SCDF's definition (first stated explicitly in the 2025 release, footnote 9):** "SCDF classifies fire fatalities as deaths resultant from the direct effects of burns, smoke inhalation, or both. Any fatality pending investigations are excluded from the annual figures released." This single footnote explains most of the reconciliation behaviour below.
2. **Initial-release totals are systematically provisional.** Documented upward revisions once investigations/coroner inquiries concluded: 2015 (4 to 7), 2016 (1 to 2), 2021 (3 to 4, per 2022 release footnote 4), 2022 (6 to 8, visible only in the data series), 2024 (5 to 7, per the Feb 2026 release as carried in press coverage). Any analysis using year-of-publication figures undercounts the most recent 12-18 months. **data.gov.sg still showed the stale 5 for 2024 at access date; treat 7 as authoritative.**
3. **The register never falls short of the official totals in any year 2005-2025.** This is the reverse of the expected undercount: bottom-up public-record compilation reaches or exceeds SCDF's official count every year. The three years where the register EXCEEDS the official total are classification effects, not register errors:
   - **2019 (register 3 vs official 1):** SCDF counts only the Bukit Batok PMD fatality. The Jalan Buroh LPG-facility worker death (MOM workplace incident) and a Jurong West flat death (police unnatural-death investigation) are outside the official series.
   - **2021 (register 7 vs official 4):** the 3-death Tuas potato-starch dust-explosion factory fire is excluded from SCDF's series (workplace/WSH classification); SCDF's four are residential incidents.
   - **2025 (register 8 vs official 6):** the press-confirmed Toa Payoh pair (Jul 2025, 2 deaths at scene) is absent from the official tally named by the Minister of State (Hougang 3 + River Valley 1 + Bukit Merah 2 = 6). No public explanation; plausibly a coronial finding of non-fire cause of death, or exclusion under the pending-investigations rule. Flagged, not reconciled away.
4. **2023 anomaly:** the SCDF 2023 statistics release text contains no fire-fatality figure at all. The 2023 total of 3 rests solely on the data.gov.sg series; the register's three incidents sum to it exactly.
5. **Series breaks to respect (never chain across):** (a) from 2019 SCDF excludes rubbish/rubbish-chute fires from total fire counts; (b) from 2024 the fire-injury count includes only hospital admissions conveyed by SCDF, so the injury series is not comparable across 2023/2024; (c) minor injury-figure discrepancy for 2021: SCDF 2021 release says 194 injuries, data.gov.sg says 193 (reconciliation CSV uses the dataset value).

## 4. Coverage limitations (documented per the coverage-bias lesson)

- **Coverage bias is real but bounded here.** Prosecuted or dramatic fires (Bedok North 2022, Geylang 2014, River Valley 2025) have deep press trails; quiet single-fatality fires often exist ONLY in SCDF's annual release itemization (e.g. Tampines St 12, Punggol Central and Gangsa Road 2018; Yishun St 22 2021; the 2020 sole fatality). Roughly a quarter of register rows rest on the SCDF PDF as sole source, with no independently locatable news coverage. Had the SCDF PDFs not been text-extracted, pure news archaeology would have found roughly 4 of 2018's 4, 0 of 2020's 1, and 0 of 2021's first three: the quiet-fire blind spot the failures lesson predicts.
- **Web search recency bias.** General search engines heavily favour 2024-2026 fire stories; 2012-2020 single-fatality incidents are close to invisible without the primary PDFs. Straits Times and CNA archives are paywalled or unfetchable, so several corroborations rest on secondary outlets (Mothership etc.).
- **2013 fatality split is approximate:** SCDF states 4 deaths across 3 named incidents without a per-incident split; the register assigns 2 to the landed-house fire by elimination (Marine Crescent confirmed as 1 by news; shipyard confirmed as 1).
- **2022 composition uncertain:** official total 8 equals the register sum, but SCDF itemized only Bedok North (3); whether the Dec 2022 Tuas gas-cylinder workplace death sits inside the official 8 is unverified. If it does not, one officially counted incident is missing from the register and one register row is an over-inclusion; net total unchanged.
- **Deaths that occur later from fire injuries** create attribution ambiguity: a fourth Bedok North (2022) victim died of her injuries in Mar 2024; which year SCDF credits that death to is not public.
- **2001-2004:** official totals exist (11, 0, 1, 7 deaths) but no per-incident itemization was pursued; pre-2005 register rows were not compiled in this pass. The 2001 figure (11) is anomalously high and worth a dedicated archival pass (NLB NewspaperSG) if the paper extends the window.
- **PDFs stored locally** in the session scratchpad only; re-download from the URLs in the CSVs if needed (all live at access date).

## 5. SCDF total-series coverage statement

Official annual fatality totals were FOUND for every year 2001-2025; none interpolated, none missing. Injuries found for 2003-2025 (2001-2002 not pulled). For 2005-2025 the per-incident register fully accounts for or exceeds every official total. Recommended paper stance: use the official series as the floor, the register as the event-resolved consequence set, and present the 2019/2021/2025 divergences as evidence that a bottom-up register plus the official series together bound the true toll, with the residual uncertainty confined to classification (workplace vs fire; pending vs concluded), not to undiscovered incidents.
