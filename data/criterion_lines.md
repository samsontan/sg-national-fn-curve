# Societal-risk (FN-curve) criterion lines - anchor points, sources, repurposing caveats

Compiled 2026-07-11 for the PMD/e-bike HDB FN-curve paper (Pipeline Topic 6). Every number below
carries a public source and access date per project rule (F010/F124). Journal-article DOIs were
verified against the Crossref API (`https://api.crossref.org/works/<DOI>`) on 2026-07-11; the
returned title/author/journal/year/volume/page fields are quoted in each citation below.

**Scope note used throughout:** all four criterion lines in this file were developed for a single
hazardous *installation* (a chemical plant, a land-use-planning consultation zone, a single
building) or, at most, a national flood-defence *system* (the Dutch case). This paper instead
benchmarks a *building-stock / population-level* FN-curve (PMD/e-bike fires aggregated across all
Singapore HDB blocks). That is a further step of repurposing beyond what most of the source
literature itself performs, and is flagged again in each caveat below and in the paper's own
Discussion section - it is not resolved by any precedent found in this search.

---

## (a) UK HSE R2P2 (2001) - "intolerability" line

**Formula / slope:** A single anchor point, extended by a conventionally-assumed slope of −1 on
log-log FN axes: F·N = constant = 10⁻² per year. This gives F(N) = 10⁻²/N.

**Anchor points**

| N (fatalities) | F, tolerability/intolerable threshold (per year) | F, "broadly acceptable" (per year, two orders of magnitude lower) |
|---|---|---|
| 1 | 1×10⁻² | 1×10⁻⁴ |
| 5 | 2×10⁻³ | 2×10⁻⁵ |
| 50 (R2P2's own anchor point) | 2×10⁻⁴ | 2×10⁻⁶ |
| 500 | 2×10⁻⁵ | 2×10⁻⁷ |

**Jurisdiction/scope:** UK land-use planning and regulatory permissioning around major hazard
(chemical/process) installations, administered by the Health and Safety Executive (HSE), including
downstream tools such as PADHI (Planning Advice for Developments near Hazardous Installations).

**Primary citation:** Health and Safety Executive (2001). *Reducing Risks, Protecting People: HSE's
Decision-Making Process* ("R2P2"). HSE Books, London.
Accessed via secondary literature 2026-07-11 (HSE's own hosted copy of R2P2 returned HTTP 404 at
`hse.gov.uk/risk/theory/r2p2.htm` when checked today; the document's content is quoted verbatim,
with page-level detail, in the peer-reviewed secondary source below).

**Secondary source used for the exact anchor point and the −1 slope convention (Crossref/venue
verified independently, non-DOI conference paper - IChemE symposium series):**
Vince, I. (2011). "Societal Risk Criteria in Land Use Planning – The Scale of 'Scale Aversion'."
*Symposium Series No. 156, Hazards XXII*, IChemE, pp. 408–410. Accessed via WebFetch of
`https://www.icheme.org/media/9272/xxii-paper-60.pdf` on 2026-07-11. This paper states verbatim:
"The tolerability line passes through the R2P2 'anchor point', namely a frequency of 2×10⁻⁴ per
year ... for single events causing 50 or more deaths," and explicitly notes that **R2P2 does not
itself endorse a slope of −1** - it cites another HSE-commissioned report (Ball, D.J. and Floyd,
P.J. (1998). *Societal Risks*, report prepared for HSE) which also does not endorse a specific
slope. The −1 slope is a near-universal convention in later HSE guidance and independent literature,
not a value stated in R2P2 itself. The "broadly acceptable" line at two orders of magnitude below
the tolerability line is likewise a convention drawn from HSE ALARP guidance, not a fixed R2P2 number.

**Repurposing caveat:** R2P2's anchor point is a societal-risk criterion for a single major-hazard
*installation* (e.g. one chemical plant), evaluated against the population around that one site. It
was never intended as a national or building-stock-level benchmark. The −1 slope is itself
contested in the literature - Vince (2011) demonstrates that an apparently "scale-neutral" slope of
−1 still embeds substantial implicit scale aversion (a factor of ~137–2465 between the maximum
tolerable Potential Loss of Life for N=1000 vs N=1 events, versus the linear 1:1000 ratio a truly
scale-neutral criterion would imply). Using this line as a benchmark for an aggregated,
population-level PMD/e-bike fire FN-curve across the whole HDB stock requires treating "the hazard"
as if it were a single institutional actor bearing risk-reduction responsibility for the whole
population, which is not the R2P2 line's original framing (there, a single operator is responsible
for a single plant's risk to its neighbours). This mismatch must be stated explicitly if this line
is used as a benchmark rather than merely an illustrative reference.

---

## (b) Dutch VROM/Bevi societal risk "oriënterende waarde" (orientation value)

**Formula:** F = 10⁻³ / N² per installation-year, applicable for N ≥ 10, with a vertical cut-off
near N = 1000 (beyond which the orientation value framework is not applied/extended).

**Anchor points** (read directly off the published curve, and consistent with the formula):

| N (fatalities) | F (per year) |
|---|---|
| 10 | 1×10⁻⁵ |
| 100 | 1×10⁻⁷ |
| 1000 (cut-off) | 1×10⁻⁹ |

**Jurisdiction/scope:** The Netherlands. Originally developed for LPG-fuelling stations, later
extended to all Seveso-type hazardous establishments regulated under Dutch external-safety policy
(Besluit externe veiligheid inrichtingen, Bevi, given legal status via the 2004 External Safety
Decree). It is a per-installation, per-year "orientation value" (reference value), not a hard legal
limit - competent authorities must justify exceedances but are not automatically barred from
permitting them.

**Primary/secondary citations (both DOIs Crossref-verified 2026-07-11):**
- Jonkman, S.N., Jongejan, R., Maaskant, B. (2011). "The Use of Individual and Societal Risk
  Criteria Within the Dutch Flood Safety Policy-Nationwide Estimates of Societal Risk and Policy
  Applications." *Risk Analysis*, 31(2), 282–300. DOI: 10.1111/j.1539-6924.2010.01502.x. Crossref
  confirms this title/author/journal/volume/issue/page set exactly. Figure 2 of this paper plots the
  curve labelled "10⁻³/n² for establishments" with the anchor points reproduced above (read off the
  published figure axes: F=10⁻⁵ at n=10, decreasing to F≈10⁻⁹ near n=1000).
- Bottelberghs, P.H. (2000). "Risk analysis and safety policy developments in the Netherlands."
  *Journal of Hazardous Materials*, 71(1–3), 59–84. DOI: 10.1016/S0304-3894(99)00072-2. Crossref
  confirms this title/author/journal/volume/issue/page set exactly (note: pages are 59–84 per
  Crossref; one secondary listing elsewhere gives 117–123, which appears to be an indexing error -
  Crossref is treated as authoritative here). This is the earliest widely-cited English-language
  primary description of the Dutch major-hazards risk policy including the 10⁻³/N² criterion.

**Repurposing caveat:** Like R2P2, this is explicitly a *per-installation* criterion - the Bevi
Reference Manual and Jonkman et al. (2011) are both clear that it was designed to bound the risk
imposed on the public by one hazardous site at a time, not to bound aggregate societal risk across
a whole country or building type. Jonkman et al. (2011) is itself a directly relevant precedent for
this paper's repurposing move: they explicitly transplant the *same* quadratic Dutch major-hazards
societal-risk criterion from the single-installation context to a national, population-level context
(nationwide flood risk across all Dutch dike rings), summing FN-curves across ~53 dike-ring areas to
get a national curve, and comparing that aggregate national curve against the *unmodified*
per-installation orientation-value line as a reference/tangent benchmark (their Fig. 11 and Section
3.3.3–3.3.4). They flag the "cumulation" (correlated simultaneous failures) issue this creates, and
note the criterion's quadratic steepness (α=2) was retained "because the same slope is used within
the Dutch major hazards policy," calling it "a risk-averse criterion" without a rigorous derivation
for the aggregate-scale use. This is a useful, directly citable precedent for treating an
installation-level line as an aggregate/benchmark reference at building-stock scale, but Jonkman et
al. themselves treat it as a comparison device rather than a validated national-scale legal
criterion, and this paper should adopt the same framing (benchmark line, not a claim that the Dutch
government endorses this line for population-level PMD-fire risk).

---

## (c) Hong Kong societal-risk (FN) criteria for Potentially Hazardous Installations

**Formulae (read off the published log-log figure, consistent with the vertical cut-off near
N=1000):**
- Tolerability ("unacceptable" boundary) limit: F = 10⁻³ / N², for 1 ≤ N ≤ 1000, with a vertical
  cut-off at N ≈ 1000 (frequencies of exceeding 1000 fatalities are treated as intolerable regardless
  of frequency).
- De minimis ("acceptable"/negligible) limit: F = 10⁻⁵ / N², same N range.
- Individual risk limit (accompanying the FN criteria, not itself an FN line): 1×10⁻⁵ per year.

**Anchor points**

| N (fatalities) | F, tolerability/unacceptable line (per year) | F, de minimis/acceptable line (per year) |
|---|---|---|
| 1 | 1×10⁻³ | 1×10⁻⁵ |
| 10 | 1×10⁻⁵ | 1×10⁻⁷ |
| 100 | 1×10⁻⁷ | 1×10⁻⁹ |
| 1000 (cut-off) | 1×10⁻⁹ | 1×10⁻¹¹ |

**Jurisdiction/scope:** Hong Kong land-use planning around Potentially Hazardous Installations
(PHIs), administered jointly by the Hong Kong Planning Department and the Environmental Protection
Department (EPD) under the Hong Kong Risk Guidelines (HKRG) / Environmental Impact Assessment
Ordinance Technical Memorandum (EIAO-TM) Annex 4 framework. Dated to 1993 per the secondary source
below; direct confirmation from an HKSAR primary PDF (e.g. EIAO-TM Annex 4 itself) was not obtained
in this pass - see caveat.

**Citation (secondary, but the chain to the primary is documented and the journal-article link in
the chain is Crossref-verified):**
Van Coile, R., Hopkin, D., Lange, D., Jomaas, G., Bisby, L. (2019). "The Need for Hierarchies of
Acceptance Criteria for Probabilistic Risk Assessments in Fire Engineering." *Fire Technology*,
55(4), 1111–1146. DOI: 10.1007/s10694-018-0746-7. Crossref confirms this title/author/journal/
volume/issue/page set exactly (checked 2026-07-11). Figure 7 of this paper ("Illustrative FN-curves
applied for land-use planning") reproduces the Hong Kong (1993) and New South Wales (2007) curves,
citing as its own source: CCPS - Frank, W. and Farquharson, J. (2009). *Guidelines for Developing
Quantitative Safety Risk Criteria*. Center for Chemical Process Safety of the American Institute of
Chemical Engineers (AIChE), New York (their reference [8]). The same paper's Section 5.2 states the
Hong Kong individual risk limit as 1×10⁻⁵ per year, and gives the New South Wales individual-risk
range as 0.5×10⁻⁶ to 5×10⁻⁵ per year (not otherwise used in this file - recorded for completeness).

**Status of primary source:** I was not able to independently retrieve or read the original Hong
Kong Planning Department / EPD document (EIAO-TM Annex 4 or the "Hong Kong Risk Guidelines") in this
session - several `epd.gov.hk` EIA report pages that reference the criteria were fetched, but they
cite Annex 4 without reproducing its numeric content, and no clean PDF of Annex 4 itself was located.
**The exact anchor points above should therefore be treated as sourced via the Van Coile et al.
(2019) / CCPS (2009) secondary chain, not independently verified against the HKSAR primary document
- flag this explicitly if the manuscript leans heavily on the HK line**, and attempt to locate and
directly cite the HKSAR primary document before submission.

**Repurposing caveat:** The HK criteria are a land-use-planning tool for a single PHI's consultation
zone, evaluated against the specific population living/working near that one site - again a
per-installation framing, now being proposed as a benchmark for a territory-wide (Singapore-wide)
building-type-specific hazard. No paper found in this search performs or discusses this specific
repurposing (installation siting criterion → national building-stock fire-hazard criterion) for the
Hong Kong line; this is a gap this paper will need to justify on its own terms in the Discussion,
not one for which a published precedent exists.

---

## (d) Fire-specific societal-risk criteria proposals in the academic literature

Two Van Coile/Hopkin-coauthored papers were located and DOI-verified via Crossref
(`https://api.crossref.org/works/<DOI>`, checked 2026-07-11). **Neither proposes a single fixed
numerical FN anchor point/slope specific to building fire safety analogous to the HSE or Dutch
lines above** - this is stated explicitly so it is not misrepresented as a settled fire-specific
criterion in the manuscript.

### (d.1) Van Coile, Hopkin, Lange, Jomaas, Bisby (2019) - the hierarchy-of-acceptance-criteria paper

**Citation:** Van Coile, R., Hopkin, D., Lange, D., Jomaas, G., Bisby, L. (2019). "The Need for
Hierarchies of Acceptance Criteria for Probabilistic Risk Assessments in Fire Engineering." *Fire
Technology*, 55(4), 1111–1146. DOI: 10.1007/s10694-018-0746-7. Crossref-verified 2026-07-11
(title/authors/journal/volume/issue/pages all match exactly).

**What it actually contributes (read from the open-access post-print, University of Ghent
institutional repository, fetched 2026-07-11):** Not a new FN line. The paper (i) explicitly states
that "no generally accepted safety targets and semi-probabilistic design methodologies currently
exist" for fire safety engineering (Section 2.2); (ii) reviews existing FN-curve criteria from other
domains (HSE, Hong Kong, New South Wales - reproducing the HK/NSW figures used in Section (c) above)
as illustrative benchmarks; (iii) proposes a five-way hierarchy of acceptance concepts for fire PRA
(AC1 comparative-implicit-tolerability, AC2 de minimis, AC3 absolute/safety-target, AC4 ALARP sensu
stricto cost-benefit, AC5 comparative-with-explicit-tolerability), with a flowchart for choosing
between them and an explicit statement of which acceptance concept places more or less
responsibility on the designer; and (iv) argues that once tolerability is established via an
FN-curve, the ALARP region should be evaluated with a risk-neutral scalar cost-benefit analysis, not
by further risk-averse weighting. Its practical relevance to this paper is as *methodology for how
to use whichever benchmark FN line is chosen*, not as a source of new numeric anchor points.

### (d.2) Mohan, Van Coile, Hopkin, Jomaas, Caspeele (2021) - risk tolerability limits methodology

**Citation:** Mohan, A.T., Van Coile, R., Hopkin, D., Jomaas, G., Caspeele, R. (2021). "Risk
Tolerability Limits for Fire Engineering Design: Methodology and Reference Case Study." *Fire
Technology*, 57(5), 2235–2267. DOI: 10.1007/s10694-021-01118-w. Crossref-verified 2026-07-11
(title/authors/journal/volume/issue/pages all match exactly).

**What it contributes (from abstract and secondary summaries only - full text was not accessible in
this session, paywalled at Springer; flagged explicitly as unverified against primary text):**
Proposes "a simple framework for setting risk tolerability limits" via literature review plus a
feedback round with international fire safety professionals, addressing a gap the authors identify
in the UK guidance document PD 7974-7:2019, and demonstrates the framework on a reference case study
of a UK office building. One secondary source (not independently verified against the paper itself)
states an individual-risk tolerability limit of 1×10⁻⁴/year and a de minimis limit of 1×10⁻⁶/year
"based on PD 7974-7:2019" in connection with this paper - **this specific number pair could not be
confirmed against the primary text in this session and must not be cited in the manuscript until
verified directly against the paper (or PD 7974-7:2019 itself).**

**Conclusion for this axis of the search:** I could not find a DOI-verifiable Van Coile/Hopkin (or
other) paper that proposes one fixed, universal societal-risk FN anchor point/slope specifically for
building fire safety, analogous to the HSE 2×10⁻⁴-at-N=50 point or the Dutch 10⁻³/N² line. The
fire-safety literature's contribution to date is methodological (a hierarchy of acceptance concepts,
and a framework for deriving *project-specific* tolerability limits) rather than a single numeric
benchmark line. **This paper's own FN-curve for PMD/e-bike HDB fires will therefore need to be
benchmarked against the installation-level lines in (a)–(c) above with the repurposing caveats
stated, since no ready-made fire-specific population-level line exists in the literature** - this
absence is itself worth stating explicitly in the manuscript's framing of its contribution.

**Repurposing caveat (applies to using (a)/(b)/(c) as fire-safety benchmarks generally):** Van Coile
et al. (2019) is the one paper in this search that explicitly discusses reusing non-fire societal-risk
FN lines (HSE, HK, NSW) as reference/comparator lines inside a fire-engineering PRA framework - but
even they do so as illustrative examples of "what an FN-curve criterion looks like," not as an
endorsement that those specific numeric lines are the correct tolerability threshold for fire risk to
building occupants. No paper found here performs the further step this paper needs: aggregating fire
risk across an entire building-stock/population (all Singapore HDB blocks) rather than a single
building or installation. This should be stated as an explicit methodological choice and limitation
in the Discussion section, not implied to be supported by precedent.

---

## Summary table

| Line | Jurisdiction / scope | Formula / slope | Anchor point(s) | Primary citation |
|---|---|---|---|---|
| UK HSE R2P2 tolerability line | Single major-hazard installation, UK land-use planning | F = 10⁻²/N (slope −1, conventional, not itself stated in R2P2) | N=50, F=2×10⁻⁴/yr | HSE (2001) *Reducing Risks, Protecting People*, HSE Books, London; slope/point per Vince (2011) IChemE Hazards XXII |
| UK HSE R2P2 "broadly acceptable" line | Same | F = 10⁻⁴/N | N=50, F=2×10⁻⁶/yr | Same, via Vince (2011) |
| Dutch VROM/Bevi orientation value | Single Seveso-type installation, Netherlands | F = 10⁻³/N², N≥10, cut-off ~N=1000 | N=10, F=1×10⁻⁵/yr; N=100, F=1×10⁻⁷/yr; N=1000, F=1×10⁻⁹/yr | Jonkman, Jongejan & Maaskant (2011), *Risk Analysis* 31(2):282–300, DOI 10.1111/j.1539-6924.2010.01502.x; Bottelberghs (2000), *J. Hazard. Mater.* 71:59–84, DOI 10.1016/S0304-3894(99)00072-2 |
| Hong Kong tolerability (unacceptable) line | Single PHI, HK land-use planning | F = 10⁻³/N², cut-off ~N=1000 | N=1, F=1×10⁻³/yr; N=1000, F=1×10⁻⁹/yr | HK (1993) via CCPS (2009) *Guidelines for Developing Quantitative Safety Risk Criteria*, AIChE, as reproduced in Van Coile et al. (2019) *Fire Technology* 55(4):1111–1146, DOI 10.1007/s10694-018-0746-7 - primary HKSAR document not independently located |
| Hong Kong de minimis line | Same | F = 10⁻⁵/N², cut-off ~N=1000 | N=1, F=1×10⁻⁵/yr; N=1000, F=1×10⁻¹¹/yr | Same chain |
| Fire-safety-specific FN line | - | **None found**: literature offers a hierarchy-of-acceptance-concepts methodology and a project-specific tolerability-derivation framework, not a fixed numeric FN line | - | Van Coile et al. (2019) DOI 10.1007/s10694-018-0746-7; Mohan et al. (2021) DOI 10.1007/s10694-021-01118-w (full text of the latter not accessed - flagged) |

---

## Open items before this file can be treated as submission-ready

1. Locate and directly cite the primary Hong Kong Planning Department / EPD document (EIAO-TM
   Annex 4 or the Hong Kong Risk Guidelines) rather than relying on the CCPS(2009)/Van Coile(2019)
   secondary chain for section (c).
2. Obtain full-text access to Mohan et al. (2021) (Fire Technology, DOI 10.1007/s10694-021-01118-w)
   and either confirm or drop the individual-risk-limit figures (1×10⁻⁴/yr tolerability, 1×10⁻⁶/yr
   de minimis) currently flagged as unverified.
3. HSE's own hosted copy of R2P2 (`hse.gov.uk/risk/theory/r2p2.htm`) returned HTTP 404 when checked
   2026-07-11 - locate a current HSE-hosted or HSE-archived PDF of R2P2 itself for a first-party
   citation link, rather than relying solely on the Vince (2011) secondary quotation.


## HK Risk Guidelines primary source - 2026-07-11 (written by MIES after verifier-agent report; resolves Open Item 1; F161)

PRIMARY DOCUMENT (fetched and read in full by the verification agent, 2026-07-11):
Planning Department, HKSAR Government. Hong Kong Planning Standards and Guidelines,
Chapter 12 (Miscellaneous), Section 4: Potentially Hazardous Installations, paras 4.1-4.9
and Figure 4.2 (Societal Risk Guidelines for Acceptable Risk Levels).
URL: https://www.pland.gov.hk/file/tech_doc/hkpsg/full/pdf/ch12_en.pdf (online edition dated
December 2025 on title page; Section 4 unchanged from the long-standing CCPHI Risk Guidelines;
Figure 4.1 companion plan dated Nov 93). Accessed 2026-07-11.

WHAT THE PRIMARY ACTUALLY SPECIFIES (Figure 4.2, log-log, F 1e-2..1e-9 /yr, N 1..10,000):
- BOTH lines slope -1 (F proportional to 1/N), NOT -2:
  - Unacceptable line: (1, 1e-3), (10, 1e-4), (100, 1e-5), (1000, 1e-6)
  - Acceptable line:   (1, 1e-5), (10, 1e-6), (100, 1e-7), (1000, 1e-8)
- Vertical cutoff at N=1000 (para 4.4.3 verbatim: vertical cut-off line at the 1000 fatality
  level extending down to a frequency of 1 in a billion years) - N>1000 unacceptable regardless.
- Regions: UNACCEPTABLE / ALARP / ACCEPTABLE. Individual risk criterion 1e-5/yr (para 4.4.2).

VERDICT: the 1e-3/N^2 and 1e-5/N^2 forms carried in this file's Section (c) and in
fn_pmd_v2.py / the PMD manuscript are a SECONDARY-SOURCE CONFLATION (CCPS 2009 / Van Coile
2019 chain) with the genuinely quadratic Dutch Bevi criterion. N=1 anchors are correct; the
slope is wrong; error grows 10x per decade of N. The SG national FN paper corrected this in
fn_sg_national_v3.py + SG_FNcurve_v3.md (2026-07-11). The PMD paper MUST correct script,
figures, table and text before submission (warning posted in this pod's state.md).
