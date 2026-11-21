# sg-national-fn-curve

Data, code and committed results for a study of societal fire risk in Singapore: a national F-N curve for building fires, 2012 to 2025, and a comparison with England on a matched basis.

> Samson Tan, Teoh Teik Toe, Paul Joseph, Khalid A. M. Moinuddin
> Manuscript in revision, November 2026. Title and citation to be added on acceptance.

## Version 3.0 replaces version 2.0

Version 2.0 (tag `v2.0-frequency-severity`) compared seven national datasets as if they covered the same incidents and populations. They do not. Version 3.0 withdraws that comparison and what was built on it:

- the seven-jurisdiction comparison matrix and the ranking drawn from it;
- the reading that Singapore has "the heaviest severity distribution";
- the test of a universal severity law and the hierarchical model;
- the pooled rate of catastrophic fires;
- the Swedish distribution by number of deaths, which no source publishes (the totals quoted in version 2.0 are from a report of the SP Technical Research Institute).

Versions 1.0 and 2.0 remain available under their tags.

## What version 3.0 does

- **Singapore register.** 62 fatal fires, 2005 to 2025, one row per fire, each with public source addresses and access dates. The primary window 2012 to 2025 holds 50 fatal building fires and 68 deaths. A death is counted with the fire that caused it, in the year of the fire.
- **Each published dataset is compared only with English record-level data re-cut to the same scope**, on matched denominators, with the English window cut to the years that overlap each dataset (in full for Singapore and the United States, in part for Australia and New Zealand). Singapore and three published datasets (United States, Australia, New Zealand) are compared with England. China and Sweden enter no comparison. Nothing is pooled and no ranking is drawn.
- **Rates are compared by number of deaths.** Over fourteen matched years Singapore's rate of single-death dwelling fires is 0.23 of England's per dwelling (95% interval 0.16 to 0.32). Its rate of fires with two or more deaths is not distinguishable from England's: 1.00 per dwelling (0.42 to 2.02).
- **The share of fatal fires with two or more deaths is a ratio of those two rates.** It is not a measure of severity. After allowances for age, household composition and death definition the single-death ratio is 0.6 to 0.7 and the two rates no longer differ significantly.
- **The curve is written two ways**: as a single power law of deaths per fire, and as two rates with a shape above one death. The data do not choose between them, and the results are reported for both.

Every number above is in `results/r2_comparison.json`.

## Contents

Data (`data/`)

- `fatal_incident_register.csv`: the Singapore register, released without the free-text description and notes fields. `fatal_incident_register_draft.csv` and `fatal_incident_register_redacted.csv` are identical copies kept so that every path cited in the manuscript resolves.
- `scdf_fatality_reconciliation.csv`: year-by-year reconciliation against the published official totals.
- `frequency_series.csv`, `exposure_national.csv`: national fire-occurrence and dwelling series for Singapore (`_draft` copies are identical).
- `home_office_fire_related_fatalities.csv`: the record-level fire-related fatalities dataset for England, year ending March 2026 edition, one row per fatality, converted from the published ODS file without change. Contains public sector information licensed under the Open Government Licence v3.0.
- `england_fatal_fires_official.csv`, `england_fire0502_cells.csv`: official numbers of fatal fires and deaths for England from table FIRE0502 (edition of 22 July 2026), summed by year and location, and at the published resolution. Contains public sector information licensed under the Open Government Licence v3.0.
- `england_fatal_dwelling_fires_by_household.csv`: fatal dwelling fires by type of household, 2010/11 to 2018/19, from Table 38 of the Home Office in-depth review.
- `exposure_comparison.csv`, `population_by_age.csv`: dwellings, households and population for the jurisdictions compared, each row with its official source.
- `revision_bases.csv`, `international_severity_bins.csv`: the published deaths-per-fire distributions of the audited datasets, one row per bin, with the audited scope and source of every row.
- `criterion_lines.md`, `register_provenance.md`, `frequency_provenance.md`, `claw_xlsx_verification.md`: sources and provenance notes.

Code (`analysis/`)

- `fn_sg_national.py`: stepped F-N curves, bootstrap bands, normalisations and sensitivity analyses for Singapore.
- `frequency_severity/fs_model.py`: the single power law (likelihood, profile intervals, fit test).
- `frequency_severity/run_all.py`, `run_revision.py`, `england_incident_level.py`: fits of the single law, Singapore tail scenarios, and England's distribution from the record-level data.
- `frequency_severity/run_r2.py`: the matched comparison by number of deaths, standardisations, capture bounds, the shape above one death, bounds on English fire sizes from official counts, and the test of the two forms where the number of large fires is known.
- `frequency_severity/r2_figs.py`: the model figures.
- `frequency_severity/sync_singapore_bins.py`: refreshes the Singapore rows of the two bins files from the register.

Results and figures: `analysis/fn_sg_results.json`, `results/*.json`, `figures/`, `figures_revision/`.

## Reproduce

Tested with Python 3.13, NumPy 2.4.3, SciPy 1.17.1 and Matplotlib 3.10.8. Install the packages in `requirements.txt`, then run from this directory:

```text
python analysis/fn_sg_national.py
python analysis/frequency_severity/sync_singapore_bins.py
python analysis/frequency_severity/run_all.py
python analysis/frequency_severity/run_revision.py
python analysis/frequency_severity/england_incident_level.py
python analysis/frequency_severity/run_r2.py
python analysis/frequency_severity/r2_figs.py
```

The scripts read only committed files. All numerical outputs are deterministic: the bootstrap bands, fit tests and simulations use fixed seeds. No number in the paper is entered by hand.

## Provenance and limits

The Singapore register was compiled from public sources (SCDF annual releases, data.gov.sg, court judgments, coronial findings and press reports) and is reconciled each year against SCDF's published totals. It meets or exceeds the official total in every year. The Singapore side of the matched comparison holds 34 single-death and 8 multiple-death fires, so the intervals are wide. The record ends at four deaths: return periods above that are properties of the assumed form of the law, not of the record.

## License

MIT for code. Data are compiled from public sources and keep their source addresses. If you use the register or the code, please cite the paper once it is published.
