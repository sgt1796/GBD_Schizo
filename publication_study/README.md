# China-US schizophrenia GBD 2023 analysis

This directory contains the reproducible analysis of schizophrenia incidence,
prevalence, and DALYs in China and the United States from 1990 through 2023.
The root-level `schizo_gbd_analysis_results` files are the reader-facing
report.

## Environment

```powershell
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe `
  --requirement analysis\requirements.txt
```

## Inputs

The canonical inputs are:

- `data/GBD_2023_schizophrenia_fine_age_China_US.csv`
- `data/GBD_2023_population_China_US.csv`
- `data/metadata/burden_export.json`
- `data/metadata/population_export.json`
- `data/metadata/structural_zero_provenance.json`

The retained source ZIP files are under `schizo/` and `addition_inputs/`.
See `DATA_INPUTS.md` for dimensions, provenance, and exclusions.

## Run the analysis

Use a new or empty output directory:

```powershell
.\.venv\Scripts\python.exe analysis\gbd_analysis.py `
  --burden-csv data\GBD_2023_schizophrenia_fine_age_China_US.csv `
  --population-csv data\GBD_2023_population_China_US.csv `
  --burden-metadata-json data\metadata\burden_export.json `
  --population-metadata-json data\metadata\population_export.json `
  --output-dir analysis\results_new
```

The checked analysis results are under `analysis/results/`:

- `tables/`: CSV tables and `analysis_tables.xlsx`
- `figures/`: main and supporting figures
- `qa/`: numerical validation and methodological notes
- `build_metadata.json`: input fingerprints and analysis settings
- `data_dictionary.csv`: generated-field definitions

## Tests and integrity checks

```powershell
.\.venv\Scripts\pytest.exe -q analysis\tests
.\.venv\Scripts\python.exe analysis\verify_analysis.py
```

## Interpretation limits

Native GBD uncertainty intervals apply to original GBD estimates. Derived
changes, ratios, segmented trends, APC summaries, and decomposition components
are point estimates because posterior draws and cross-estimate covariance were
not available. Country and sex comparisons are descriptive and do not establish
causal effects of health systems or policies.
