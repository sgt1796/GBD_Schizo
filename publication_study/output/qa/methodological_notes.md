# Internal methodological notes

This file is generated from the current analysis outputs. Apparent disagreements are retained and investigated; they are not automatically described as validation failures.

## Cross-method disagreements

### China, Female, Incidence

- Directions: segmented 1990-2023 = increase; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; APC net drift = decrease.
- Likely factors: calendar window and piecewise-versus-endpoint estimand; opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Female, Prevalence

- Directions: segmented 1990-2023 = increase; ASR endpoint 1994-2023 = increase; selected-age crude endpoint 1994-2023 = increase; APC net drift = decrease.
- Likely factors: APC global drift versus crude endpoint change and APC constraints; opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Female, DALYs

- Directions: segmented 1990-2023 = increase; ASR endpoint 1994-2023 = increase; selected-age crude endpoint 1994-2023 = increase; APC net drift = decrease.
- Likely factors: APC global drift versus crude endpoint change and APC constraints; opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Male, Incidence

- Directions: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; APC net drift = decrease.
- Likely factors: opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Male, Prevalence

- Directions: segmented 1990-2023 = increase; ASR endpoint 1994-2023 = increase; selected-age crude endpoint 1994-2023 = increase; APC net drift = increase.
- Likely factors: opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Male, DALYs

- Directions: segmented 1990-2023 = increase; ASR endpoint 1994-2023 = increase; selected-age crude endpoint 1994-2023 = increase; APC net drift = decrease.
- Likely factors: APC global drift versus crude endpoint change and APC constraints; opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Female, Incidence

- Directions: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; APC net drift = decrease.
- Likely factors: opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Male, Incidence

- Directions: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; APC net drift = increase.
- Likely factors: APC global drift versus crude endpoint change and APC constraints; opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Male, Prevalence

- Directions: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; APC net drift = decrease.
- Likely factors: opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Male, DALYs

- Directions: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; APC net drift = decrease.
- Likely factors: opposing age-specific local drifts hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

## APC window sensitivity

- China, Female, DALYs: primary 1994-2023 net drift = -0.001645%/year (decrease); sensitivity 1990-2019 net drift = 0.001900%/year (increase). Both magnitudes should be inspected before interpreting the sign change.

## Decomposition age-bin sensitivity

- China, Female, Incidence: maximum component shift = 50.0% of total change; rank stable = False.
- China, Male, Incidence: maximum component shift = 49.0% of total change; rank stable = False.
- United States of America, Female, Incidence: maximum component shift = 17.6% of total change; rank stable = True.

These accounting components are descriptive, not causal effects. Exact all-age closure is a QA property and does not prove that the source age bins are sufficiently granular for ageing attribution.
