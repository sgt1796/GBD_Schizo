# Analysis specification

## Status

This document records the analysis decisions used for reproducible rebuilds of
the China--United States schizophrenia study. It is a versioned analysis
specification, not a claim of prospective registration. Any material change
should be documented here before the corresponding manuscript is submitted.

Current version: 0.2 (1 September 2026). Version 0.2 extends the demographic
decomposition from ages 15+ to all ages after confirming that the available
0--14, five-year, and 70+ age bins reproduce every reported all-age count.

## Research question

How did sex-specific schizophrenia incidence, prevalence, and disability-
adjusted life-year (DALY) burden differ between China and the United States
from 1990 through 2023, and how much of the change in counts is attributable to
population size, age composition, and age-specific rates?

## Primary scope

- Locations: China and United States of America.
- Sex strata: female and male, using the GBD sex variable.
- Period: 1990--2023.
- Outcomes: incidence, prevalence, and DALYs.
- Descriptive trend estimands: all-age counts and age-standardized rates per
  100,000.
- Decomposition estimands: all-age count changes, using the available 0--14
  group, 5-year age groups from 15--19 through 65--69 years, and the 70+
  terminal group.

YLDs are audited against DALYs and are not reported as a separate outcome when
the estimates are numerically identical. The GBD Percent metric,
probability-of-death extract, and risk-factor extract are outside the study
question and are excluded.

## Analysis principles

1. Native GBD 95% uncertainty intervals are reported only for native GBD
   estimates.
2. Ratios, endpoint changes, trend summaries, and decomposition components are
   point estimates unless posterior draws are available.
3. Annual GBD posterior means are correlated modeled estimates. They are not
   treated as independent observations for confirmatory hypothesis testing.
4. Piecewise log-linear models are descriptive summaries. Knot number and
   location may be selected reproducibly by an information criterion, but the
   resulting slopes do not receive inferential confidence intervals or
   significance labels in the primary manuscript.
5. Country and sex contrasts are ecological and descriptive. No health-system,
   policy, diagnostic, or causal mechanism is inferred from them.

## Demographic decomposition

The primary decomposition compares 1990 with 2023 and averages marginal
contributions over all six replacement orders for population size, age
composition, and age-specific rates. Components must close to the reconstructed
count change within numerical tolerance, and reconstructed counts must agree
with the source's reported all-age counts. Analyses beginning in 2000 and 2010,
plus annual- and five-year-chained decompositions, assess endpoint and path
sensitivity. The broad 0--14 and 70+ terminal groups limit interpretation of
within-bin age-pattern changes.

The submission build must use population estimates from the same official GBD
2023 release as the burden estimates. Population reconstructed from count/rate
pairs is permitted only for visibly provisional engineering builds.

## Secondary age-period-cohort analysis

The secondary analysis is restricted to incidence, ages 15--69 years, and six
five-year periods from 1994--1998 through 2019--2023. Only estimable drift and
nonlinear curvature summaries are interpreted. The entire 2019--2023 period is
removed in the endpoint-period sensitivity analysis. APC results are
descriptive and cannot identify independent causal age, period, or cohort
effects.

## Optional validation analyses

- Approximate weighting derived from marginal GBD UIs may be used only as a
  sensitivity analysis, not as sampling-variance weighting.
- Official NCI Joinpoint output may be compared with the open implementation as
  an optional validation analysis. It is not required for the reproducible
  primary analysis and must never be claimed unless generated in the official
  software.
- If GBD posterior draws become available, trend contrasts and decomposition
  uncertainty should be recomputed draw by draw and promoted to the primary
  analysis.

## Submission gate

A submission build requires the authenticated official GBD 2023 population
export and a successful full pipeline validation. Until then, all generated
manuscripts and decompositions remain explicitly provisional.
