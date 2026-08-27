# sg-national-fn-curve

Data, model code and reproduction package for:

> **A frequency-severity decomposition of national building-fire societal risk: evidence from seven jurisdictions**
> Samson Tan, Teoh Teik Toe, Paul Joseph, Khalid A. M. Moinuddin
> Manuscript in preparation, August 2026. Citation to be updated on acceptance.

National building-fire F-N curves have been drawn descriptively but never modelled. This
repository holds everything behind a paper that changes that: a decomposition of a
jurisdiction's societal fire risk into a Poisson fatal-fire rate and a conditional severity
distribution (a discrete power law with upper limit), estimated by maximum likelihood with
goodness-of-fit testing across seven jurisdictions, about 25,000 fatal fires in all. The
Singapore per-incident register that anchors the study, its reconciliation against official
totals, the international severity dataset, the model code, and every committed result are
here, so any reported number can be checked or extended without asking us.

**Version 2.0 supersedes the v1.0 package** (tagged `v1.0-descriptive`), which accompanied an
earlier descriptive version of the paper. Everything in v1.0 is retained; v2.0 adds the
seven-jurisdiction severity dataset and the frequency-severity model layer.

## What the model finds

- **Rate and severity are uncoupled.** Fatal-fire rates per million population span a
  fourteen-fold range while the scalar deaths-per-fire ratio is nearly flat; conditional on a
  fatal fire, the probability of two or more deaths is 0.255 in Singapore against 0.060 in
  England (ratio 4.2, 95% CI 2.3-6.7). The lowest-rate jurisdiction carries the heaviest
  severity distribution.
- **A universal severity law is rejected** as an estimated spread, sd(log b) = 0.23
  (0.13-0.48), with exponents spanning about 2.7 to 4.3 on a common within-regime basis.
- **Catastrophic fires are regime failures, not tail draws.** Under England's own severity
  law fitted to its 2,354 ordinary fatal fires, Grenfell Tower was a 5,000-to-14,000-year
  event; pooling the four such regime failures across the seven jurisdictions' 9,379 million
  person-years gives an order-of-magnitude anchor for jurisdictions with no catastrophe in
  their record (Singapore: roughly 400 years, with all assumptions stated).
- **A two-number estimator**: the severity exponent is recoverable from the multi-fatality
  share alone (cap-insensitive above Nmax = 10), so any national table publishing fatal-fire
  and multi-fatality event frequencies can be converted into model parameters without
  per-incident microdata.

## Contents

| Path | What |
|---|---|
| `data/fatal_incident_register.csv` | The Singapore register. 58 incidents, 2005-2025, one row per fatal fire, each with a public source URL and access date |
| `data/scdf_fatality_reconciliation.csv` | Year-by-year register-to-SCDF reconciliation, including the years where the register exceeds the official total |
| `data/international_severity_bins.csv` | **New in v2.0.** Per-incident severity bins for seven jurisdictions (Singapore, China, England, USA, Sweden, Australia, New Zealand) with scope, provenance grade and source per row |
| `data/international_fn_curves_v2.json` | **New in v2.0.** The seven per-capita F-N curves with full derivation notes |
| `data/exposure_national.csv`, `data/frequency_series.csv` | Singapore dwelling-stock, population and fire-occurrence series |
| `data/register_provenance.md`, `data/frequency_provenance.md`, `data/criterion_lines.md` | Provenance and reference-line documentation |
| `analysis/fn_sg_national.py` | v1.0 descriptive analysis: stepped curves, bootstrap bands, leave-one-out envelope, normalisations |
| `analysis/frequency_severity/fs_model.py` | **New in v2.0.** The model: DPLDwUL severity (interval likelihoods, profile CIs, bootstrap GoF), common-basis rule, LRT machinery |
| `analysis/frequency_severity/run_all.py` | **New in v2.0.** One command reproduces every model-based quantity in the paper |
| `results/` | **New in v2.0.** Committed outputs: per-jurisdiction fits, common basis, universality tests, hierarchical posterior, Singapore tail scenarios and monitoring Bayes factors, Grenfell dragon-king test, catastrophic pooling, derived indices |
| `figures/` | All paper figures: the v1.0 descriptive set and the v2.0 model set (framework, worked example, risk map, fitted overlays, Grenfell, tail scenarios, monitoring, placement chart) |
| `supplementary/SUPPLEMENTARY.md` | Supplementary material S1-S5: search log, sensitivity tests, reconciliation notes, severity-model supplement, international dataset derivations |

## Reproduce

```
pip install numpy scipy matplotlib
python analysis/fn_sg_national.py                     # v1 descriptive layer
python analysis/frequency_severity/run_all.py         # v2 model layer -> results/
```

Both read only committed files and write only into `results/`. No number in the paper is
entered by hand.

## Use the estimator on your own jurisdiction

If your national statistics publish (a) the annual frequency of fatal fires and (b) the
annual frequency of fatal fires killing two or more people, the ratio (b)/(a) is the
multi-fatality index MFI, and the severity exponent solves
`1 - 1/sum(i^-b, i=1..Nmax) = MFI` (insensitive to Nmax above 10). See
`results/derived_indices.json` for worked values and the standard-error formula, and the
paper's placement chart for where your jurisdiction lands. Sources reporting only death
totals (the CTIF World Fire Statistics format) cannot feed the estimator; publishing the
two event frequencies would fix that, which is one of the paper's recommendations.

## Provenance and limits

The Singapore register was compiled entirely from public sources (SCDF annual releases,
data.gov.sg datasets, MHA records, press reports) with per-row provenance, and is reconciled
annually against SCDF's published totals; it meets or exceeds the official total in every
year, so residual uncertainty is classificatory, not missing incidents. International
severity data are graded per jurisdiction (per-incident, published counts, reconstructed, or
percentage-derived) and every model claim in the paper carries the grade's consequences.
The severity model is a working approximation: it fails formal goodness-of-fit in the two
largest datasets and the paper says so.

## License

MIT for code. Data compiled from public sources; per-row source URLs retained. If you use
the register, the international dataset or the model, please cite the paper above.
