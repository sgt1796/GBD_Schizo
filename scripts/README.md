# Schizophrenia burden analysis

This directory contains the analysis used for the journal manuscript comparing
schizophrenia burden in China and the United States from 1990 to 2023.

## Run

From the repository root:

```powershell
python scripts/prepare_burden.py
python scripts/prepare_population.py
python scripts/run_analysis.py
```

The preparation scripts extract and combine the retained IHME GBD 2023 source
archives in `GBD_data/`. The analysis script reads the two canonical CSV files
in that directory and writes manuscript tables, figures, and documents under
`output/` at the repository root.

## Contents

- `prepare_burden.py`: prepares the age-, sex-, year-, and measure-specific
  burden dataset.
- `prepare_population.py`: prepares the matching population dataset.
- `apc_analysis.py`: age-period-cohort calculations used by the main analysis.
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
