# China-US schizophrenia GBD 2023 publication study

This directory contains the reproducible analysis and manuscript build for the
China-US comparison. The root-level `prepare_schizo_inputs.py` script and the
raw archives under `schizo/` are upstream inputs to this pipeline. Other
root-level analysis scripts and reports are earlier exploratory work and are
not publication-pipeline inputs.

All commands below are PowerShell commands run from the repository root
(`GBD_Schizo`). The workflow is tested with CPython 3.12, and direct Python
dependencies are pinned in `publication_study/requirements.txt`.

## 1. Create the environment

```powershell
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe `
  --requirement publication_study\requirements.txt
```

## 2. Prepare the burden inputs

```powershell
.\.venv\Scripts\python.exe prepare_schizo_inputs.py
.\.venv\Scripts\python.exe publication_study\prepare_production_burden.py
.\.venv\Scripts\python.exe publication_study\prepare_production_population.py
.\.venv\Scripts\python.exe publication_study\verify_reproducibility.py
```

The preparation step reads three GBD 2023 archives that contribute to prepared
outputs:

- `schizo/IHME-GBD_2023_DATA-8a231d74-1.zip` and
  `schizo/IHME-GBD_2023_DATA-d3e13215-1.zip` are combined into
  `prepared_inputs/cause_all.csv`;
- `schizo/IHME-GBD_2023_DATA-8d0bd0df-1.zip` produces the separately retained
  probability-of-death table, which the publication analysis excludes.

When present, `schizo/IHME-GBD_2023_DATA-5ef7a575-1.zip` is inspected solely to
document why its single risk branch is excluded; it is not an analysis input.
The preparation report records source row counts, hashes, dimensions,
duplicates, exclusions, and validation warnings in
`prepared_inputs/GBD_1990_2023_schizophrenia_preparation_report.txt`. The
verifier checks the Python environment and these source/prepared artifacts; it
does not run the analysis.

The second preparation command builds the production fine-age burden file from
the preserved split downloads under `addition_inputs/`. The first download
contained an accidental `65-74 years` aggregate and the second supplied the
required `70-74 years` group. The script removes the aggregate and unused
Percent rows, normalizes `<5 years` to `0-4 years`, validates every required
dimensional cell, and writes
`data/GBD_2023_schizophrenia_fine_age_China_US.csv` deterministically. It does
not create or infer the mandatory official population file.

The third preparation command reads the preserved GBD 2023 Population/Number
ZIP, normalizes the IHME `<5 years` label to `0-4 years`, validates exactly
2,720 positive location-sex-age-year cells and their source uncertainty bounds,
and writes `data/GBD_2023_population_China_US.csv` with an explicit `GBD 2023`
release marker.

### IHME data terms

The raw ZIP archives and derived CSV files remain subject to IHME's data terms;
any license applied to repository code does not override those terms. IHME/GHDx
states that downloadable data may be used, shared, modified, or built upon by
non-commercial users under the
[IHME Free-of-Charge Non-Commercial User Agreement](https://www.healthdata.org/sites/default/files/files/free-of-charge_non-commercial_user_agreement.pdf).
Anyone using or redistributing the ZIP or CSV data must review the current IHME
terms, remain within the permitted use, and cite the applicable GBD release
using the IHME citation supplied with the source archives. This repository does
not represent the data as unrestricted or Creative Commons data; commercial-use
questions must be directed to IHME.

## 3. Run the automated tests

```powershell
.\.venv\Scripts\python.exe -m pytest publication_study\tests -q
```

## 4. Build a provisional package

Use this mode for code verification and manuscript development while the
official population export is unavailable:

```powershell
.\.venv\Scripts\python.exe publication_study\publication_analysis.py `
  --allow-proxy-population `
  --output-dir publication_study\build\provisional

.\.venv\Scripts\python.exe publication_study\build_documents.py `
  --analysis-dir publication_study\build\provisional

.\.venv\Scripts\python.exe publication_study\verify_reproducibility.py `
  --analysis-dir publication_study\build\provisional
```

This mode reconstructs population from matched burden count/rate pairs. It is
deterministic, labels the population and decomposition as provisional, and sets
`submission_ready` to `false`. A successful verifier exit means that the
package is internally complete; it does **not** mean that the study is ready for
submission.

The fine-age schizophrenia export has exact zero Number and Rate values in both
`0-4 years` and `5-9 years`. Their separate populations therefore cannot be
reconstructed from burden ratios. Do not run the fine-age input with proxy mode
or invent a split of the combined residual; use the official population export
for the fine-age production analysis.

Trend selection is deterministic: the primary open model chooses zero to two
breakpoints by BIC and does not use a resampling-based significance test.

## 5. Build the production package

Production requires matching fine-age schizophrenia burden and official
population exports from GBD 2023, plus completed provenance sidecars for both.
Their exact contracts are in `EXTERNAL_INPUTS_REQUIRED.md`; use
`gbd_export_metadata_template.json` for burden and
`gbd_population_export_metadata_template.json` for population.

```powershell
.\.venv\Scripts\python.exe publication_study\publication_analysis.py `
  --burden-csv data\GBD_2023_schizophrenia_fine_age_China_US.csv `
  --population-csv data\GBD_2023_population_China_US.csv `
  --population-release "GBD 2023" `
  --burden-metadata-json data\metadata\burden_export.json `
  --population-metadata-json data\metadata\population_export.json `
  --output-dir publication_study\build\production

.\.venv\Scripts\python.exe publication_study\build_documents.py `
  --analysis-dir publication_study\build\production

.\.venv\Scripts\python.exe publication_study\verify_reproducibility.py `
  --analysis-dir publication_study\build\production `
  --require-submission-ready
```

Do not build production over an older provisional directory. A fresh output
directory prevents stale tables or documents from being mistaken for results
from the current inputs. Both analysis and document commands reject nonempty
output directories.

### Optional official NCI validation

Each run exports 12 series plus settings under
`<analysis-dir>/nci_joinpoint_inputs/`. If registered NCI software is available,
run those inputs, normalize the official output as described in
`EXTERNAL_INPUTS_REQUIRED.md`, and rebuild with:

```powershell
--nci-results-csv data\nci_results_normalized.csv
```

This adds separate `nci_validation_*` tables alongside the primary independent
segmented-trend tables; it is not a submission gate. The file
`nci_results_template.csv` is a header-only schema template; it is not an
analysis result and must never be passed unchanged.

As checked on 2026-09-01, NCI distributes version 6.1.0. The importer accepts a
single consistently reported, parseable version and records it in metadata;
old generated settings that name 6.0.1 must be regenerated with the current
code. Record the version actually used and validate any version transition.
Never relabel one version as another or describe the independent Python model
as NCI output.

Posterior GBD draws are also optional. If they become available, they would
materially strengthen uncertainty propagation for contrasts, changes, and
decomposition, but the current point-estimate analysis remains reproducible
without them and explicitly reports that limitation.

## Input contracts

The default burden input is `prepared_inputs/cause_all.csv`. A replacement may
be supplied with `--burden-csv`. It must contain at least these columns:

```text
location_name, sex_name, age_name, measure_name, metric_name,
year, val, lower, upper
```

The analysis retains China and United States of America; Female and Male;
1990-2023; required Incidence, Prevalence, and DALY outcomes; optional YLDs for
the duplication audit; and Number and Rate metrics.
The provisional prepared file contains the all-age count, age-standardized
rate, and a broad endpoint age partition. Production requires Number and Rate
panels in consecutive five-year ages from 0-4 through 90-94 plus 95+ for all
three outcomes, as well as All ages Number and Age-standardized Rate panels for
the primary endpoint and trend analyses.

The canonical population schema is:

```text
location_name,sex_name,age_name,year,population,gbd_release
```

It requires exactly 2,720 unique rows: 2 locations x 2 sexes x 20 age groups x
34 years. The groups run from `0-4 years` through `90-94 years` plus `95+ years`.
Values must be finite and positive. Exact
accepted labels and the
optional NCI normalized-result schema are listed in
`EXTERNAL_INPUTS_REQUIRED.md`.

## Generated package

Each analysis directory contains machine-readable CSV/XLSX tables, a data
dictionary, main and supplementary figures, QA results, optional NCI-validation
inputs, and `build_metadata.json`. The document builder adds the manuscript,
supplement, statistical methods appendix, GATHER checklist, and a document
manifest under `documents/`.

The current primary trend schema is descriptive. Segmented, trajectory, APC,
and cross-analysis tables contain no p values, q values, regression confidence
intervals, or parallelism-test decisions. APC outputs separately report net
drift, local drift, longitudinal age curves, period RRs, and cohort RRs for the
1994-2023 primary and 1990-2019 sensitivity windows. The decomposition is
consistently identified as all-age: audit fields use `all_age_groups` and
`all_age_year_cells`, while decomposition tables use
`all_age_count_start_reconstructed` and
`all_age_count_end_reconstructed`. The separate
`all_age_count_reconstruction.csv` table verifies that age-specific
population-rate cells sum to reported all-age counts. The production partition
has 20 groups; the provisional input has 13 and is visibly barred from
submission. `decomposition_age_bin_sensitivity.csv`,
`cross_analysis_consistency.csv`, and `qa/methodological_notes.md` retain
sensitivity flags and cross-method disagreements.
Run into a fresh directory so removed permutation-era tables cannot survive
beside the current schema.

The authoritative final checks are:

- `qa/validation_summary.json` for data and output invariants;
- `build_metadata.json` for input provenance and submission status; and
- `verify_reproducibility.py --require-submission-ready` for a single failing
  command if the official-population, fine-age, or provenance gate remains
  unresolved.

## Document rendering and visual QA

`build_documents.py` generates DOCX files and `document_manifest.json`; it does
not generate PDFs. The current DOCX files are the authoritative editable
documents. Any pre-existing files under `rendered/`, or PDFs elsewhere in the
repository, may be stale and must not be treated as current deliverables.

For final submission, render new PDFs from the exact final DOCX files using the
target word-processing environment, then inspect every page for pagination,
table overflow, figure resolution, fonts, captions, cross-references, and
provisional labels. A PDF becomes an authoritative submission artifact only
after that fresh render and visual QA. The reproducibility verifier checks the
DOCX package but does not perform or certify visual QA.

## Statistical scope

- Primary outcomes are incidence, prevalence, and DALYs. Percent metrics and
  probability-of-death outputs are excluded.
- When YLDs are supplied, DALYs and YLDs are audited for numerical identity and
  YLDs are not reported twice. A wholly absent YLD panel is recorded as an
  optional audit that was not available; a partial or non-identical panel fails
  validation.
- Native 95% GBD uncertainty intervals are reported only for original GBD
  estimates. Descriptive segmented trends and AAPCs have no model confidence
  intervals or p values.
- Derived ratios, changes, and decomposition components are point estimates
  because posterior draws are absent.
- The built-in BIC-selected segmented model is an independent implementation
  and must not be described as NCI Joinpoint.
- APC analysis is secondary and restricted to incidence, equal five-year ages
  (preferably 10-69; provisional fallback 15-69), two six-period windows, and
  identifiable drifts/nonlinear relative-risk functions. The production
  all-age decomposition uses the complete fine-age partition.

See `EXTERNAL_INPUTS_REQUIRED.md` before treating any build as submit-ready.
