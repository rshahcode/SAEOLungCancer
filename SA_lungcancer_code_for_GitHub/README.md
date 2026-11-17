# South Asia early onset lung cancer analysis code

This repository contains analysis code only. It deliberately excludes GBD downloads,
standard population values, generated numerical results, manuscripts, images, plotting
code, figure links, and personal directory paths. Inputs must be obtained separately
under the relevant IHME access and use terms. Nothing is downloaded automatically.

## Run

Install R packages `readxl`, `sandwich`, and `nlme`. The validated environment used
R 4.4.1. Run from any directory:

```sh
Rscript /path/to/code/run_analysis.R /path/to/external/inputs /path/to/external/results
Rscript /path/to/code/tests.R
```

Input and output directories must be outside this repository. Do not commit generated
results or downloads. The script refuses to write inside its code directory.

## Manual extraction specifications

Use GBD 2023, all population, cause **Tracheal, bronchus, and lung cancer** (cause 426).
Retain full precision and the original lower and upper limits. Use the exact filenames
below, renaming files without changing their contents. If a download is split, combine
parts with one header and do not duplicate records. The Results Tool is an interactive
source; this repository does not claim to provide a download API or an extraction
script that was not used.

1. `annual_rates.csv`: measures Incidence and Deaths; metric Rate; all years
   1990 through 2023 inclusive; seven separate age groups 15–19, 20–24, 25–29,
   30–34, 35–39, 40–44, and 45–49 years; Male, Female, and Both separately.
   Locations: Global; South Asia; Bangladesh; Bhutan; India; Nepal; Pakistan;
   Central Europe, Eastern Europe, and Central Asia; High-income; Latin America
   and Caribbean; North Africa and Middle East; Southeast Asia, East Asia, and
   Oceania; Sub-Saharan Africa. Expected records: 18,564. Raw `_name` columns
   and cleaned names are both accepted. Pooled 15–49 rates are not accepted.
2. `standard_population.xlsx`: original IHME GBD 2021 demographic Appendix Table
   S14 workbook. The first sheet must retain its original age labels in column A
   and standard-population percentages in column B. The script locates the seven
   age labels, rather than embedding their downloaded values.
3. `source_age_specific_changes.csv`: GBD-provided relative change for 1990–2023,
   Incidence and Deaths, Rate, same seven separate age groups, three sexes;
   all 13 locations listed for the annual rates. Expected records: 546.
   Keep `year_start`, `year_end`, `val`, `lower`, `upper`; fractional changes
   are multiplied by 100, not treated as percentage points.
4. `absolute_burden_1.csv` and `absolute_burden_2.csv`: the two non-overlapping
   parts of the Number export; pooled age 15–49 years; Incidence, Deaths,
   YLLs (Years of Life Lost), and DALYs (Disability-Adjusted Life Years);
   1990, 2019, and 2023; three sexes; Global, High-income, South Asia,
   Bangladesh, Bhutan, India, Nepal, Pakistan. Expected combined records: 288.
5. `risk_male.value.csv`, `risk_female.value.csv`, `risk_onlyboth.value.csv`:
   Deaths, Percent, age 15–49 years, 2023, South Asia and five countries;
   one file for each sex. Use cleaned columns `measure,location,sex,age,cause,
   rei,metric,year,val,upper,lower`. Values are fractions on the 0–1 scale.
   Include all 26 original exported nested categories, including All risk factors.
   Expected combined records: 468. The pipeline preserves source limits, including
   negative marginal lower limits, without silently clipping them.
6. `risk_both.sex.change.csv`: combined three-sex GBD relative change export for
   the same risk fractions, 1990–2023, with cleaned `location,sex,rei,val,lower,upper`
   columns and all other original fields. `risk_male.change.csv` and
   `risk_female.change.csv` are the dedicated exports used to cross-check it.
   There were 458 supplied change records. Missing asbestos changes remain
   unavailable; they are never replaced by zero.

## Methods and interpretation

Direct standardization uses fixed GBD 2021 weights restricted and renormalized within
each of ages 15–49, 15–39, and 40–49. Regional estimates are standardized from regional
age-specific rates, not averaged over country rates. Both-sex rates use the directly
exported both-sex estimates. Endpoint changes are ratios of standardized point estimates.

Weighted age-specific lower and upper limits are **approximate bounds**, not formal
95% uncertainty intervals. Joint draws were not available. No endpoint-change tests
or derived uncertainty intervals are produced. Source-provided age-specific, absolute
number and risk intervals retain their original interpretation. Risk outputs are
fractions, not risk-attributable death counts or standardized mortality rates; nested
categories overlap and are not summed.

EAPC is an exploratory log-linear regression summary of the annual standardized point
estimates. Reference confidence intervals use Newey–West covariance, Bartlett lag 3,
no prewhitening, finite-sample adjustment and a t reference with n−2 degrees of freedom.
These are conditional regression confidence intervals, not propagated GBD uncertainty.
Diagnostics and sensitivity outputs are retained for transparency, not model selection
for statistical significance. Alternative AR1 fits were added after examining residual
persistence; near-boundary fits are flagged as unstable. No inference is fitted to the
five observations in 2019–2023. The code contains no Joinpoint or age-period-cohort model.

`primary_eapc_1990_2023.csv` contains the 78 primary trend summaries, including lag
sensitivities and Benjamini–Hochberg adjusted p values. Regional age/period sensitivity
and AR1 results are separate outputs. Review these together; do not selectively report
only intervals that exclude zero. Curvature can make a single full-period slope an
incomplete description. Regression intervals do not establish significance of endpoint
percentage changes, country ranks, or sex differences.

## Checks

The pipeline rejects duplicate or incomplete age-rate keys, wrong causes or metrics,
invalid intervals and incompatible grids. It reconciles absolute mean counts across
sexes and countries without summing uncertainty limits. Input checksums and the R
session information are written to the external output directory. `tests.R` uses
synthetic in-memory values only. Publication requires separately confirming the GBD
risk-cause definitions and use/publication requirements; code validation does not
replace those checks.

