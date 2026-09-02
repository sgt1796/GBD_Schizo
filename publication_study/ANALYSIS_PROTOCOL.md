# Analysis specification

## Status and versioning

This document records the analysis decisions used for reproducible rebuilds of
the China--United States schizophrenia study. It is a versioned analysis
specification written after data access, not a prospective registration. Any
material change must be recorded here before the corresponding manuscript is
designated final.

Current version: 0.3 (1 September 2026).

Version 0.3 separates the study into three analysis layers, replaces the former
embedded APC routine with an estimable-function module, defines independent
1994--2023 and 1990--2019 APC windows, requires fine-age production inputs, and
adds segmented-model and decomposition age-bin sensitivity analyses. The prior
descriptive trend redesign remains frozen unless a reproducible bug is found.

## Research question

How did sex-specific schizophrenia incidence, prevalence, and disability-
adjusted life-year (DALY) burden differ between China and the United States
from 1990 through 2023, and how do aggregate trajectories, age-specific temporal
patterns, and demographic accounting jointly characterize changes in burden?

## Three complementary analysis layers

### 1. Primary: descriptive burden and trend analysis

- Locations: China and United States of America.
- Sex strata: female and male, using the GBD sex variable.
- Period: 1990--2023.
- Outcomes: incidence, prevalence, and DALYs.
- Estimands: all-age counts, age-standardized rates per 100,000, endpoint
  changes, descriptive country/sex contrasts, and descriptive segmented-curve
  summaries.
- Curve specification: log-linear curves with zero to two breakpoints selected
  by BIC, with breakpoint location counted in the penalty and a four-year
  minimum segment. Alternative breakpoint counts, minimum segment lengths,
  calendar windows, and rate-scale curves are sensitivity analyses.
- Formal trajectory hypothesis tests, permutation tests, p values, q values,
  and conditional regression confidence intervals are excluded.

### 2. Secondary: age--period--cohort analysis

APC is retained intentionally because it addresses age-specific temporal
structure that cannot be recovered from all-age or age-standardized trends.
It is not a binary validation of segmented AAPC.

- Outcome: incidence.
- Preferred ages: 10--69 years in consecutive five-year groups when complete
  fine-age input is available.
- Provisional fallback: 15--69 years in consecutive five-year groups. The
  aggregated 0--14 group and open-ended older groups are not used in APC.
- Primary window: 1994--2023, six complete five-year periods.
- Sensitivity window: 1990--2019, six complete five-year periods.
- Estimands: net drift, local drift by age, longitudinal age curve, period
  relative risks, and cohort relative risks.
- Identification: the exact age = period - cohort dependency is handled by two
  linear trends plus orthogonal nonlinear age, period, and cohort curvatures.
  Three unconstrained linear effects are not reported as uniquely identifiable.
- Interpretation: period and cohort patterns are model-derived population
  patterns and are not causal exposure effects.

APC and segmented curves estimate different quantities. APC uses age-specific,
population-weighted cells and a global constrained model over a different age
range and calendar window; segmented trends summarize age-standardized annual
series with piecewise curves. Numerical agreement is neither required nor
expected. Every directional disagreement is retained in the cross-analysis
audit and investigated in terms of estimand, age coverage, window, age
standardization, population weighting, model form, and identification
constraints.

### 3. Secondary: demographic decomposition

The decomposition compares changes in counts attributable to:

1. population-size change;
2. age-structure change; and
3. age-specific-rate change.

For population vector P and age-specific rate vector R, burden is the sum of
P_a R_a across a complete age partition. Contributions are averaged over all
six replacement orders (Shapley averaging), so interactions are allocated
symmetrically and the three components close to the reconstructed count change
within numerical tolerance.

The production analysis requires the finest common five-year age groups from
0--4 through 90--94 years plus 95+ years in both burden and official population
exports. Provisional engineering builds may use the available 0--14, five-year
15--69, and 70+ partition, but that partition cannot clear the submission gate.

The main comparison is 1990 versus 2023. Analyses beginning in 2000 and 2010,
annual and five-year chained decompositions, and an intentionally collapsed
four-group decomposition assess path and age-bin sensitivity. A sensitivity
flag is raised when a component changes sign, its absolute-magnitude rank
changes, or its shift reaches 10% of the total count change.

Exact reconstruction of reported all-age counts is a QA check, not proof that
the available age bins are sufficiently granular for ageing attribution.
Decomposition components are accounting quantities and are not causal effects.

## Input requirements

The production build requires matching GBD 2023 burden and population exports
for China and the United States, both sexes, all years 1990--2023, and the same
fine age partition. Burden exports must contain Number and Rate metrics for
incidence, prevalence, and DALYs. Raw IHME exports must be preserved unchanged,
with query/export metadata, retrieval dates, and SHA-256 hashes.

Population reconstructed from matched count/rate pairs is permitted only in a
visibly provisional engineering build. When official population is supplied,
the pipeline records cell-level discrepancies from the reconstructed
denominator rather than silently substituting values.

When a YLD panel is included in an export, it is audited against DALYs and is
not reported as a separate outcome when the estimates are numerically
identical. YLD is not a required production input because it is not a study
outcome. The GBD Percent metric,
probability-of-death extract, and available risk-factor extract are outside the
research question and are excluded.

## Uncertainty and inference principles

1. Native GBD 95% uncertainty intervals are reported only for native GBD
   estimates.
2. Endpoint changes, ratios, AAPCs, APC estimands, trajectory contrasts, and
   decomposition components are point estimates unless posterior draws are
   available.
3. GBD point estimates are modeled and correlated across years, strata, and
   outcomes. Incidence, prevalence, and DALYs are not independent replications.
4. Residual diagnostics describe model fit and autocorrelation; they are not
   converted into GBD-level inference.
5. No claim such as a significant country, sex, period, cohort, or trajectory
   difference is made from marginal uncertainty intervals or model-conditional
   standard errors.
6. If posterior draws become available, uncertainty propagation is implemented
   as a separate upgrade without changing the descriptive estimands silently.

## Optional validation analyses

- Approximate weighting derived from marginal GBD uncertainty intervals may be
  used only as a sensitivity analysis, not as sampling-variance weighting.
- Official NCI Joinpoint output may be compared with the open segmented model
  as optional validation. It is not a primary dependency and must never be
  claimed unless generated in the official software.
- Posterior-draw analyses may propagate uncertainty through endpoints,
  selected trend summaries, APC estimands, and decomposition where feasible.

## Submission gate

A build may report `submission_ready=true` only when all of the following hold:

- authenticated official GBD 2023 population is used;
- matching fine-age burden Number and Rate panels are complete;
- both APC windows pass matrix, interval, cohort-index, reference, and
  synthetic-recovery tests;
- fine-age decomposition and age-bin sensitivity are complete;
- cross-method disagreements are recorded and explained;
- all authoritative tables, figures, documents, and provenance records are
  rebuilt from final code in a clean directory;
- the manuscript and supplements match the machine-readable outputs; and
- the complete reproducibility and visual-QA sequence passes.

Until those conditions are satisfied, generated documents and decomposition
results remain explicitly provisional.
