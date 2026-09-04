# Data inputs

## Required dimensions

- Locations: China and United States of America
- Sexes: Female and Male
- Years: every integer from 1990 through 2023
- Outcomes: incidence, prevalence, and DALYs
- Metrics: Number and Rate
- Burden ages: five-year groups from 0-4 through 95+, plus All ages and
  Age-standardized
- Population ages: the 20 five-year groups from 0-4 through 95+
- Population measure: persons

## Canonical files

The normalized burden and population files are stored under `data/`, with
retrieval metadata and source-file hashes under `data/metadata/`. Raw source
archives remain unchanged under `schizo/` and `addition_inputs/`.

The source export contains 2,720 exact zero Number/Rate cells. The canonical
input preserves them unchanged. The verified pattern is incidence at ages 0-9
and 80+, and prevalence and DALYs at ages 0-9, across both locations, sexes,
metrics, and all years. The export does not establish whether these are
biological zeros or a model-support convention.

## Exclusions

- The GBD Percent metric is excluded because its denominator and age basis
  differ by outcome.
- Probability of death is outside this nonfatal-burden analysis.
- The available risk-factor archive covers only one limited risk branch and is
  not a complete comparative risk assessment.
- YLDs are not a required outcome. If present, they are audited against DALYs.

## Preparation

```powershell
.\.venv\Scripts\python.exe prepare_schizo_inputs.py
.\.venv\Scripts\python.exe analysis\prepare_production_burden.py
.\.venv\Scripts\python.exe analysis\prepare_production_population.py
```

These commands normalize labels, validate dimensional completeness and
uncertainty bounds, retain provenance, and write the two canonical CSV files.
