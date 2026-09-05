# Statistical analysis specification

## Study design

This ecological time-series study compares schizophrenia incidence, prevalence,
and disability-adjusted life years (DALYs) in China and the United States from
1990 through 2023 using GBD 2023 estimates. Analyses are stratified by sex.

## Descriptive estimates

For each country, sex, outcome, and year, the analysis reports age-standardized
rates and reconstructed counts. Endpoint summaries quantify absolute and
relative change between 1990 and 2023. Country and sex contrasts are calculated
on both absolute and relative scales.

## Temporal trends

Piecewise log-linear regression estimates annual percent change within segments
and average annual percent change over the complete series. Candidate
breakpoints are selected using the prespecified information-criterion procedure
implemented in `run_analysis.py`. Sensitivity analyses examine autoregressive
errors, alternative breakpoint specifications, practical-stability thresholds,
and exclusion of 2020--2023.

Incidence is the principal outcome for age-period-cohort analysis. Prevalence and
DALYs are descriptive disease-burden measures and are not interpreted as
incident-event risks.

## Age-period-cohort analysis

Five-year age and period groups define birth cohorts by period midpoint minus
age midpoint. The constrained log-linear model reports descriptive age curves,
age-specific temporal slopes, period curvature, cohort curvature, and a global
period slope. Population-weighted estimates are primary; equal-weight estimates
and a 1990--2019 time window assess sensitivity. These custom summaries are not
presented as conventional NCI estimable functions.

## Decomposition

Changes in counts are decomposed into population growth, population ageing, and
age-specific rate change. The primary analysis uses consecutive five-year age
groups available across both countries, both sexes, all years, and all three
outcomes. Annual and five-year chained decompositions, alternative factor order,
and broader age-bin definitions assess sensitivity. An incidence-only analysis
restricted to ages with positive source rates is supplementary.

## Uncertainty and interpretation

Uncertainty intervals supplied by GBD are propagated for reported endpoint
estimates. Trend and decomposition results are treated as model-based
descriptions of the GBD estimates. Country, sex, and outcome comparisons are
correlated views of the same modeled data and are not treated as independent
replications. Interpretation emphasizes direction, magnitude, consistency
across analyses, and sensitivity results rather than isolated significance
tests.

The analysis is descriptive and cannot establish causality. Differences may
reflect epidemiology, demography, health-system detection, data availability,
and GBD modeling assumptions.
