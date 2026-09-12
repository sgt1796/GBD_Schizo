# Schizophrenia burden analysis

This directory contains the analysis comparing schizophrenia burden in China
and the United States from 1990 to 2023.

## Run

From the repository root:

```powershell
python scripts/prepare_burden.py
python scripts/prepare_population.py
python scripts/run_analysis.py
```

The APC step requires R and `Rscript` on `PATH` (or `RSCRIPT_PATH` set to its
executable). The age-period-cohort step uses pinned R code published by
Rosenberg et al. and distributed by the National Cancer Institute.

The preparation scripts extract and combine the retained IHME GBD 2023 source
archives in `GBD_data/`. The analysis script reads the two canonical CSV files
in that directory and writes tables and figures under `output/` at the
repository root. The main figures include the age-specific China-US prevalence
gap, ageing of the prevalent population, and the endpoint demographic
decomposition. Segmented trends and cumulative annual decomposition are
supplementary figures.

## Contents

- `prepare_burden.py`: prepares the age-, sex-, year-, and measure-specific
  burden dataset.
- `prepare_population.py`: prepares the matching population dataset.
- `nci_apc_analysis.py`: validates and groups incidence counts and population.
- `nci_apc_runner.R`: runs the published APC estimable functions and exports point estimates.
- `nci_apc_reference.R`: pinned source code from the National Cancer Institute.
- `run_analysis.py`: produces the reported analyses, sensitivity analyses,
  tables, and figures.
- `METHODS.md`: statistical analysis specification.
- `DATA_INPUTS.md`: source archive and canonical dataset description.

Install the pinned runtime packages with:

```powershell
python -m pip install -r scripts/requirements.txt
```

The input data are IHME estimates rather than individual-level observations.
Results should therefore be interpreted as comparative population-level trends,
not causal effects or individual risks.
