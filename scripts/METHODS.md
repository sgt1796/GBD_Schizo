# Statistical analysis specification

## Study design

This ecological time-series study compares schizophrenia incidence, prevalence,
and disability-adjusted life years (DALYs) in China and the United States from
1990 through 2023 using GBD 2023 estimates. Analyses are stratified by sex.

## Descriptive estimates

For each country, sex, outcome, and year, the analysis reports age-standardized
rates and reconstructed counts. Endpoint summaries quantify absolute and
relative change between 1990 and 2023. Country and sex contrasts are calculated
on both absolute and relative scales. Annual point-estimate rate ratios describe
the evolution of the China-to-United States and male-to-female disparities.
Age-specific endpoint changes describe the distribution of rate changes across
five-year age groups; percent changes are left undefined when either endpoint
rate is zero.

## Temporal trends

Piecewise log-linear regression estimates annual percent change within segments
and average annual percent change over the complete series. Candidate
breakpoints are selected using the prespecified information-criterion procedure
implemented in `run_analysis.py`. Sensitivity analyses examine autoregressive
errors, alternative breakpoint specifications, practical-stability thresholds,
and exclusion of 2020--2023.

Incidence alone is used for age-period-cohort analysis. Prevalence and DALYs
are disease-burden measures, not incident-event counts.

## Age-period-cohort analysis

The age-period-cohort incidence analysis uses the published `apc2()` R code of
Rosenberg et al. (*Cancer Epidemiol Biomarkers Prev*, 2014), distributed with
the National Cancer Institute's APC Web Tool (source pinned at commit
`9e86d92d3f95b76b610a37f6e3c539ee126b5efd`). The input is GBD modeled
incident counts summed within ages 10--14 through 65--69 and six five-year
periods from 1994--1998 through 2019--2023;
the offset is the corresponding summed GBD population. The implementation
fits its constrained age-period-cohort log-rate model and returns its standard
estimable functions: net drift, local drifts, longitudinal age rates, period
rate ratios, and cohort rate ratios. Reference categories are the
implementation's defaults. The equal-width 1990--1994 through 2015--2019
window is a sensitivity analysis.

The GBD counts are modeled posterior means, often fractional, not observed
independent Poisson event counts. Therefore only the APC point-estimate
functions are retained. The implementation's confidence intervals and Wald
tests are discarded; neither they nor GBD posterior uncertainty can be
interpreted without posterior draws and cross-cell dependence information.

## Decomposition

Changes in counts are decomposed into population growth, population ageing, and
age-specific rate change. The primary analysis uses consecutive five-year age
groups available across both countries, both sexes, all years, and all three
outcomes. Annual and five-year chained decompositions, alternative factor order,
and broader age-bin definitions assess sensitivity. An incidence-only analysis
restricted to ages with positive source rates is supplementary.

## Prevalence country gap, population ageing, and annual decomposition

Age-by-year heatmaps display the log2 ratio of age-specific prevalence rates in
China versus the United States, separately for females and males. The diverging
scale is centered at equality, and age-year cells with undefined ratios are
masked. The ageing analysis reports the percentage and number of prevalent
cases aged 65 years or older over time. The annual chained decomposition is
plotted cumulatively to show when population growth, population ageing, and
age-specific rate change contributed to the change in prevalent counts.

## Uncertainty and interpretation

Uncertainty intervals supplied by GBD are propagated for reported endpoint
estimates. Derived contrasts and decomposition components use posterior means
because posterior draws were not
available. Trend and decomposition results are treated as model-based
descriptions of the GBD estimates. Country, sex, and outcome comparisons are
correlated views of the same modeled data and are not treated as independent
replications. Interpretation emphasizes direction, magnitude, consistency
across analyses, and sensitivity results rather than isolated significance
tests.

The analysis is descriptive and cannot establish causality. Differences may
reflect epidemiology, demography, health-system detection, data availability,
and GBD modeling assumptions.
