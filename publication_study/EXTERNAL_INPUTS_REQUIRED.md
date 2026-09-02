# External inputs and optional validation

Two matching authenticated exports remain mandatory before submission: a
fine-age schizophrenia burden export and official population denominators from
the same GBD 2023 release. Each needs a provenance sidecar recording the exact
query, retrieval date, export identifier, and hashes of the preserved raw
files. Official NCI output and posterior GBD draws are optional
validation/strengthening inputs. None of these inputs may be inferred,
relabeled, or manufactured by the Python pipeline.

The burden/population data are governed by IHME terms independently of any code
license. Users of redistributed source ZIPs, prepared CSVs, or a newly obtained
population export must comply with the current
[IHME Free-of-Charge Non-Commercial User Agreement](https://www.healthdata.org/sites/default/files/files/free-of-charge_non-commercial_user_agreement.pdf)
and cite the applicable GBD release. Nothing in this repository grants
unrestricted or Creative Commons rights to IHME data.

## 1. Mandatory: matching fine-age GBD 2023 exports

Download both schizophrenia burden estimates and population estimates from the
same GBD 2023 results release. Preserve the original files unchanged under a
non-generated raw-data directory.

### Required dimensions

- Locations: `China`; `United States of America`
- Sexes: `Female`; `Male`
- Years: every integer from `1990` through `2023`
- Burden ages: `0-4 years`, `5-9 years`, `10-14 years`, every consecutive
  five-year group through `90-94 years`, and `95+ years`, plus `All ages` and
  `Age-standardized`
- Population ages: the 20 fine groups from `0-4 years` through `95+ years`
- Burden measures: schizophrenia incidence, prevalence, and DALYs
- Burden metrics: Number and Rate
- Population measure/metric: population number in persons

Do not substitute `0-14 years` or `70+ years` for their component groups in the
production inputs. The former 13-bin files remain suitable only for provisional
engineering and sensitivity work.

For each burden outcome, the 20 fine groups must contain both Number and Rate.
The frozen primary analysis additionally requires `All ages` Number and
`Age-standardized` Rate for every location-sex-year cell. Other exported
age-metric combinations may be retained, but they do not replace those two
summary panels.

### Canonical CSV schema

Use the header in `population_input_template.csv`:

```text
location_name,sex_name,age_name,year,population,gbd_release
```

There must be one and only one population row per location-sex-age-year
combination: 2,720 rows in total (2 locations x 2 sexes x 20 ages x 34 years).
`population` must be numeric, finite, positive, and in
persons (not thousands). Set `gbd_release` to `GBD 2023` on every row. The
loader also recognizes the column aliases `location`, `sex`, `age`, and one of
`pop`, `value`, or `val` for population, but canonical names are preferred for
an auditable archive.

Run the file with:

```powershell
--population-csv data\GBD_2023_population_China_US.csv `
--population-release "GBD 2023" `
--burden-csv data\GBD_2023_schizophrenia_fine_age_China_US.csv `
--burden-metadata-json data\metadata\burden_export.json `
--population-metadata-json data\metadata\population_export.json
```

`--population-release` is a provenance assertion supplied by the analyst; it
does not authenticate the CSV. Use `gbd_export_metadata_template.json` for the
burden sidecar and `gbd_population_export_metadata_template.json` for the
population sidecar. Replace every placeholder and list the preserved raw files
using repository-relative paths and actual lowercase SHA-256 hashes. The
pipeline verifies the query dimensions and hashes. Do not use population
reconstructed from burden count/rate pairs in a submitted analysis.

The 20 fine population groups cover all ages and are used in decomposition. As a
cross-source QA check, summing each age-specific reconstruction
(`population x rate / 100,000`) must agree with the corresponding reported
all-age count within the pipeline's documented numerical/rounding tolerance.
The secondary APC analysis uses the consecutive groups from 10--14 through
65--69 years when complete. It excludes open-ended ages and retains a visible
15--69 fallback only for provisional builds.

## 2. Optional: official NCI Joinpoint validation

NCI requires the end user to accept its Terms of Use and register before
downloading the application. The Python analysis exports 12 input series plus
`input_manifest.csv` and `analysis_settings.json` under
`<analysis-dir>/nci_joinpoint_inputs/`. Those files identify the recorded
series and settings; they are not evidence that NCI was run.

Official NCI output is useful for validating or replacing the repository's
independent segmented-trend implementation, but it is not required to open the
submission gate.

### Version compatibility issue

As checked on 2026-09-01, the NCI website lists Joinpoint **6.1.0** as the
current desktop and command-line release. The current importer does not
hard-code that version: it requires every normalized row to carry the same
parseable version and records the value in build metadata. Older generated
settings in saved output directories may still name **6.0.1**; regenerate those
inputs with the current code rather than relabeling old artifacts.

Before importing output from a different version, verify the analysis settings
and normalization fields, record any changed defaults, and rerun the complete
package. Do not change a version label, and do not describe output from the
independent Python segmented model as NCI output. No official NCI output is
included or claimed by this repository at present. Confirm the current release
at <https://surveillance.cancer.gov/joinpoint/> and follow the registered
download terms at <https://surveillance.cancer.gov/joinpoint/download>.

### Normalized CSV contract

Copy `nci_results_template.csv` to a new file and populate it from official NCI
exports. The template contains only a header and is not itself valid input.
Keep all listed columns, using blanks for fields that do not apply to a row.
Every row must state the software version actually used.

Accepted identifiers are:

- `location_name`: `China` or `United States of America`
- `sex_name`: `Female` or `Male`
- `measure_name`: `Incidence`, `Prevalence`, or `DALYs`
- `analysis_type`: `trend`, `segment`, `comparison`, or `fitted`

The file must contain:

- 12 `trend` series (one per location-sex-measure combination), with
  `joinpoint_count`, `joinpoint_years`, `aapc`, `aapc_lower_ci`, and
  `aapc_upper_ci`;
- one or more `segment` rows for every one of the 12 series, with
  `segment_index`, `start_year`, `end_year`, `apc`, `apc_lower_ci`, and
  `apc_upper_ci`;
- either no `comparison` rows or exactly 12 optional rows, with
  `comparison_family`, `stratum`, `group_a`, `group_b`, and
  `parallelism_p_value`: six China-versus-US comparisons within sex and outcome,
  and six Female-versus-Male comparisons within location and outcome; and
- exactly 408 unique `fitted` rows (12 series x 34 years), with `year` and
  `fitted` for every year from 1990 through 2023.

For `comparison_family=country`, use sex as `stratum`, `China` as `group_a`,
and `United States of America` as `group_b`. For
`comparison_family=sex`, use location as `stratum`, `Female` as `group_a`, and
`Male` as `group_b`.

Normalization is a transcription step, not a re-estimation step. Archive the
native NCI session/project and output files alongside the normalized CSV so the
mapping can be audited. Software confidence limits and p values are preserved
as reference-only fields; they are not promoted to primary inference because
GBD posterior draws and cross-estimate covariance are unavailable.

## 3. Optional: posterior GBD draws

Posterior draws would allow uncertainty to be propagated through changes,
ratios, trajectory contrasts, and decomposition. Their absence is explicitly
reported, and they are not a submission gate for the current descriptive
point-estimate study. If draws are obtained, preserve draw identifiers and the
same location-sex-age-year definitions, then revise and validate the statistical
workflow before making draw-based interval claims.

## Submission gate

After running the final analysis and document builds in a fresh directory, run:

```powershell
.\.venv\Scripts\python.exe publication_study\verify_reproducibility.py `
  --analysis-dir publication_study\build\production `
  --require-submission-ready
```

The command must exit successfully. `build_metadata.json` must report
`population_status` as `official_GBD_2023`, `fine_age_burden_validated` and
`source_metadata_complete` as true, and `submission_ready` as true. The QA
summary must agree. NCI-import status may be false because it is optional.
Passing this mechanical gate does not replace scientific review of the source
export, model settings, or manuscript claims.

It also does not certify rendered documents. The automated build produces DOCX
files, not authoritative PDFs. Render the final DOCX files afresh and complete
page-by-page visual QA before designating any PDF as a submission artifact;
older files under `rendered/` may be stale.
