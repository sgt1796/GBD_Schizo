# Internal methodological notes

This file is generated from the current analysis outputs. Apparent disagreements are retained and investigated; they are not automatically described as validation failures.

## Cross-method disagreements

### China, Female, Incidence

- Practical labels: segmented 1990-2023 = practically stable; ASR endpoint 1994-2023 = practically stable; selected-age crude endpoint 1994-2023 = decrease; custom APC global period slope = practically stable.
- Likely factors: age coverage, age standardization, and population weighting; custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Female, Prevalence

- Practical labels: segmented 1990-2023 = practically stable; ASR endpoint 1994-2023 = practically stable; selected-age crude endpoint 1994-2023 = increase; custom APC global period slope = practically stable.
- Likely factors: age coverage, age standardization, and population weighting; custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Female, DALYs

- Practical labels: segmented 1990-2023 = practically stable; ASR endpoint 1994-2023 = practically stable; selected-age crude endpoint 1994-2023 = increase; custom APC global period slope = practically stable.
- Likely factors: age coverage, age standardization, and population weighting; custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Male, Incidence

- Practical labels: segmented 1990-2023 = practically stable; ASR endpoint 1994-2023 = practically stable; selected-age crude endpoint 1994-2023 = decrease; custom APC global period slope = practically stable.
- Likely factors: age coverage, age standardization, and population weighting; custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Male, Prevalence

- Practical labels: segmented 1990-2023 = practically stable; ASR endpoint 1994-2023 = practically stable; selected-age crude endpoint 1994-2023 = increase; custom APC global period slope = practically stable.
- Likely factors: age coverage, age standardization, and population weighting; custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### China, Male, DALYs

- Practical labels: segmented 1990-2023 = practically stable; ASR endpoint 1994-2023 = practically stable; selected-age crude endpoint 1994-2023 = increase; custom APC global period slope = practically stable.
- Likely factors: age coverage, age standardization, and population weighting; custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Female, Incidence

- Practical labels: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; custom APC global period slope = decrease.
- Likely factors: opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Male, Incidence

- Practical labels: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; custom APC global period slope = practically stable.
- Likely factors: custom APC global period slope versus crude endpoint change and model constraints; opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Male, Prevalence

- Practical labels: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; custom APC global period slope = decrease.
- Likely factors: opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

### United States of America, Male, DALYs

- Practical labels: segmented 1990-2023 = decrease; ASR endpoint 1994-2023 = decrease; selected-age crude endpoint 1994-2023 = decrease; custom APC global period slope = decrease.
- Likely factors: opposing age-specific slopes hidden by aggregate summaries.
- Interpretation: The compared methods target different age coverage, weighting, windows, and model functions. The discrepancy is retained as a substantive cross-method finding unless synthetic recovery or input QA fails.

## APC window sensitivity

All custom APC global-period-slope practical labels agreed across the two windows.

## Decomposition age-bin sensitivity

- China, Female, Incidence: maximum component shift = 50.0% of total change; rank stable = False.
- China, Male, Incidence: maximum component shift = 49.0% of total change; rank stable = False.
- United States of America, Female, Incidence: maximum component shift = 17.6% of total change; rank stable = True.

These accounting components are descriptive, not causal effects. Exact all-age closure is a QA property and does not prove that the source age bins are sufficiently granular for ageing attribution.
