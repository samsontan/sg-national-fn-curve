# Supplementary Material

Supplementary material for:

> **A frequency-severity decomposition of national building-fire societal risk: evidence from seven jurisdictions**
> Samson Tan, Teoh Teik Toe, Paul Joseph, Khalid A. M. Moinuddin

Section numbers refer to the main article. Every model-based quantity here regenerates from `analysis/frequency_severity/run_all.py` in this repository; committed outputs are in `results/`.

## S1. Search log: national and building-stock F-N work from fire statistics

Supporting detail for Section 2.2. Each work below was checked against the question "does this construct a national or building-stock F-N curve from official fire statistics?". None does. Wang, Lu and Li, the sole located precedent, is described in the main text and is not repeated here.

Other national fire-statistics work is methodologically close but not F-N-framed. Manes and Rush re-derive structural fire-response statistics from England and from United States data against the guidance document PD 7974-7; that is a frequency-and-severity treatment rather than a fatality-exceedance curve. Sekizawa's international comparison of fire-death rates across countries is descriptive epidemiology with no exceedance construction, a reading re-verified against the full text. Matellini et al. model United Kingdom dwelling fires with a three-part Bayesian network calibrated on national 1990 to 2011 statistics; its outputs are single-dwelling fatality probabilities and cost utilities, not a frequency-of-N-fatalities curve. The full text was checked directly because a national F-N construction there would have been a second precedent; it is not one. Regional high-rise fire-risk models exist, notably a spatial Markov-chain and indicator-system model for Beijing; these characterise relative regional risk without an F-N curve or societal-risk criterion. F-N and Bayesian-network hybrids are appearing in adjacent infrastructure such as tunnels. Our own group has quantified structural-fire ignition frequency from Australian statistics and built a probabilistic fire-risk model for high-rise residential buildings; both are frequency-side or building-level rather than national F-N. Per-incident multiple-fatality data do exist outside the F-N literature: the NFPA publishes an annual catastrophic multiple-death fire series for the United States, and national fatality studies for Australia and New Zealand tabulate deaths per incident (Section 2.4); none constructs an F-N curve. No F-N curve built from Singapore fire statistics was located in the peer-reviewed record; SCDF's releases are operational statistics, not academic analysis.

One further lead remains open: an unattributed figure titled as an F-N curve for fire events in the Netherlands, circulating on a scholarly repository, could not be traced to a parent publication despite systematic searching. The precedent survey is therefore stated as bounded rather than exhaustive. The searches behind it were systematic: national F-N and societal-risk constructions were searched for Japan (including a Japanese-language pass over the national fire agency's annual statistical reports), Korea, Hong Kong, Taiwan, Macau, the Netherlands, the Nordic countries, the United States, Canada, the United Kingdom, Australia, New Zealand, and the ten Southeast Asian jurisdictions from Malaysia to Laos, with every promising candidate resolved by a full-text read where access allowed. The searches ran over English-language indexed databases; non-English grey literature and fire-service technical reports in the region's languages were not covered, and the novelty statements of Sections 1 and 2.5 are bounded accordingly.

The search also grounds the paper's gap statements G1 to G4 (Section 2.5): what exists at national scale is descriptive (rates, ratios and drawn curves), and no located work estimates a severity model, compares severity distributions across jurisdictions, treats the catastrophic tail as a separate regime, or supplies an operational layer.

## S2. Window sensitivity and classification scope tests

Supporting detail for Section 4.5.

### S2.1 Window sensitivity

The rate is sensitive to the observation window. The first half of the primary window (2012 to 2018, T = 7) gives F(N >= 1) = 3.0 per year and F(N >= 2) = 0.71 per year from 21 fatal events. The second half (2019 to 2025, T = 7) gives 3.7 and 1.0 per year from 26 events. The extended register (2005 to 2025, T = 21) gives 2.6 and 0.67 per year from 55 events. The higher second-half figures should not be read as a trend: the event counts are too small to separate a real change from sampling noise, and register completeness is itself time-varying because SCDF's per-incident itemisation improves over the period. The trend test reported in Section 4.14 (23 of 47 events in the last five of fourteen years, p of about 0.06) is the sharper statement of the same caution.

### S2.2 Classification scope tests

Two classification choices were scope-tested against the national curve. Excluding the 2021 Tuas dust-explosion fire from the building envelope, treating it like the open process-plant and on-vessel events, lowers F(N >= 1) from 3.36 to 3.29 per year and F(N >= 3) from 0.29 to 0.21 per year, and reduces the non-residential scope from 6 to 5 events; the national tail at N = 4 is unaffected because the four-death Geylang event is unrelated to this reclassification. The delayed-death reclassification of the 2022 Bedok North fire is reported in Section 4.2. Neither reclassification changes the paper's qualitative reading.

## S3. Full register-to-SCDF fatality reconciliation, 2005-2025

The complete year-by-year reconciliation is `data/scdf_fatality_reconciliation.csv` in this repository, with per-year notes distinguishing transcription-based matches, conditional matches, and the three genuine excess years (2019, 2021, 2025). Section 4.1 of the main text summarises; Section 4.16 carries the CTIF cross-reference. The reconciliation is on deaths; an incident-level reconciliation awaits the official event-resolved series recommended in Section 4.20.

## S4. Severity-model supplement

Supporting detail for Sections 3.1 to 3.6 and 4.9 to 4.15. All values regenerate from `run_all.py`; committed outputs in `results/`.

### S4.1 Per-jurisdiction fits at the observed support

`results/fits_per_jurisdiction.json`. Exponents fitted at each source's observed support with profile-likelihood intervals and parametric-bootstrap G-test p-values. Notes carried in the main text: Singapore's pass has low power at n = 47; England's failure is produced entirely by Grenfell Tower; the US distribution is percentage-derived, so its fit test is rounding-dominated and its interval is convolved with the rounding envelope; Sweden's reconstruction makes its fit test non-informative by construction; Australia's open ">4" bin moves the exponent by under 0.05 when placed at N = 5, 6 or spread to 8.

### S4.2 Common-basis (within-regime, untruncated) fits

`results/common_basis.json`. The comparative basis of Table 9: untruncated MLE on the within-regime record, England excluding Grenfell (n = 2,354), China excluding its three events of N >= 30 (n = 8,951).

### S4.3 Universality likelihood-ratio tests

`results/universality_lrt.json`. Common-exponent LRTs: all seven jurisdictions (b_common = 3.14, LRT = 675, df = 6, p ~ 1e-142); the five Western datasets (p ~ 1e-24); Singapore and China (p = 0.24, a low-power non-rejection given Singapore's interval). Reported here rather than in the main text because, at a combined n of about 25,000 with residual misspecification in the largest datasets, the p-values overstate the crispness of the conclusion; the main text's universality statement is the hierarchical spread sd(log b) = 0.23 (95% CI 0.13 to 0.48).

### S4.4 Hierarchical posterior

`results/hierarchical_posterior.json`. Random-walk Metropolis, 40,000 iterations, 10,000 burn-in; log b_j ~ Normal(mu, sigma^2), mu ~ Normal(log 3, 1), sigma ~ HalfNormal(0.5). Pooled values are a stated sensitivity only (Section 3.4).

### S4.5 Singapore tail scenarios and monitoring

`results/singapore_tail.json`. Per-Nmax profile refits (Table 10) and the marginal monitoring Bayes factors over Singapore's own untruncated posterior: 2.9, 6.6 and 13 over one, two and three event-free decades.

### S4.6 Grenfell dragon-king test and catastrophic pooling

`results/grenfell_dragon_king.json` and `results/catastrophic_pooling.json`. The within-regime England fit (b = 4.28 on N <= 6), the implied return period of a 70-or-more-death event, the look-elsewhere check, and the pooled any-cause catastrophic rate with leave-one-out return periods for Singapore.

### S4.7 Derived indices and the two-number estimator

`results/derived_indices.json`. Common-basis lambda x E[N] against observed deaths (the first-moment consistency check of Section 4.15, which has no tail power), truncation ratios, F(>= 10) upper bounds, estimator examples, and the delta-method standard error.

### S4.8 Split-sample check

Fitting 2012-2020 (24 events, b = 2.47 at the observed support) and predicting 2021-2025: observed events 23 against a 95% prediction interval [6, 24]; observed deaths 33 against [7, 33]; realised deaths-per-event 1.43 against 1.31 predicted (a 9% miss, the honest out-of-sample number quoted in Section 4.15). Severity exponents across halves are formally indistinguishable (equality LRT p = 0.38); the frequency side is marginal, not passed (Section 4.14).

## S5. International severity dataset: derivation notes

`data/international_severity_bins.csv` holds the per-incident severity bins for all seven jurisdictions with per-row scope, provenance grade and source. Derivations and caveats:

- **Singapore.** This repository's register, building-envelope scope, 47 events 2012-2025. Reconciled annually against SCDF totals (S3), which bounds single-fatality undercount, the bias that would most inflate apparent severity.
- **England.** Home Office in-depth review 2010/11-2018/19, Table 2, per-incident. Shoreham Airshow (N = 11, outdoor aircraft crash) excluded as out of building scope; Grenfell Tower (N = 71) included in the dataset and excluded from the within-regime fit.
- **United States.** USFA Topical Fire Report Series percentage distributions converted to pseudo-counts (about 1,900 fatal residential fires per year); rounding uncertainty is comparable to sampling uncertainty and is carried in the reported interval. Residential only; the 5+ bin is open.
- **Sweden.** Reconstructed from the published aggregate constraints (1,183 incidents, 1,292 deaths) and an earlier period's proportional split; the per-N split carries roughly +-15% uncertainty at N >= 2, and fit tests on this row are non-informative by construction.
- **Australia.** Coates et al. Table 8 exact counts; the ">4" bin is entered as an interval; the source states its final years under-report by 7-23%. Preventable residential scope excludes deliberate fires, which steepens the fitted exponent (Section 3.3).
- **New Zealand.** Lilley and McNoe per-incident counts; unintentional domestic scope (same steepening caveat). The source's stated death total (118) differs from its implied count (120); disclosed, not resolved.
- **China.** Wang, Lu and Li Table 4, cumulative counts converted to bins; the source merges 4-5 and 6-9, which the interval likelihood handles exactly. All-fires scope (not only buildings) and unverifiable single-fatality capture; both caveats carried in Sections 3.3 and 4.9.

Scope differences act on the fitted exponent with a known sign (excluding deliberate or non-building events steepens it), which is why cross-jurisdiction severity claims in the main text are anchored on the Singapore-England near-scope-matched pair.
