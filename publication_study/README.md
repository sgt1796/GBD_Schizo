# China-US schizophrenia GBD 2023 publication study

This directory contains the reproducible redesign of the original exploratory
report. The legacy files in the repository root are preserved.

## Production build

```powershell
python publication_study\publication_analysis.py `
  --population-csv path\to\GBD_2023_population.csv `
  --population-release "GBD 2023" `
  --nci-results-csv path\to\nci_results_template.csv `
  --output-dir publication_study\output

python publication_study\build_documents.py `
  --analysis-dir publication_study\output
```

## Provisional build

When the official population export is unavailable, the pipeline can be tested
with populations reverse-engineered from burden count/rate pairs:

```powershell
python publication_study\publication_analysis.py `
  --allow-proxy-population `
  --output-dir publication_study\output
```

This mode is for engineering and manuscript drafting only. It sets
`submission_ready=false`, labels every decomposition output as provisional, and
does not satisfy the study protocol.

## Statistical safeguards

- Primary outcomes are incidence, prevalence, and DALYs. Percent metrics and
  probability-of-death outputs are excluded.
- DALYs and YLDs are audited for identity; YLDs are not reported twice.
- Native endpoint 95% GBD UIs remain distinct from model confidence intervals.
- Derived ratios and changes are point estimates because posterior draws are absent.
- The included segmented trend model uses residual permutation selection and is
  not described as NCI Joinpoint.
- APC inference is limited to incidence, equal five-year periods, net/local
  drift, and nonlinear curvature contrasts.

## Tests

```powershell
python -m pytest publication_study\tests -q
```

See `EXTERNAL_INPUTS_REQUIRED.md` for the two submission gates.
