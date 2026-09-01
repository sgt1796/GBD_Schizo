# External inputs required before submission

The code and provisional publication package are complete, but two externally
licensed/authenticated inputs cannot be truthfully manufactured by the pipeline.

## 1. Official GBD 2023 population export

Download population estimates from the same GBD 2023 results release as the
burden data. Required dimensions:

- Locations: China; United States of America
- Years: 1990-2023
- Sexes: Female; Male
- Ages: 15-19, 20-24, ..., 65-69, and 70+ years
- Metric: population number

Save a CSV with the columns shown in `population_input_template.csv`. Either add
`gbd_release=GBD 2023` to every row or run with
`--population-release "GBD 2023"`. The strict production command refuses a
population file without that release marker.

## 2. Official NCI Joinpoint 6.0.1 output

NCI requires the end user to accept its Terms of Use and register with personal
details. The inputs and exact settings are generated in
`output/nci_joinpoint_inputs/`. Run those files in the official application,
then normalize the exported results to `nci_results_template.csv`.

The independent segmented-regression output included in the provisional build
is deliberately labeled **not NCI Joinpoint**. It must not be represented as an
NCI result.

## Submission gate

`output/build_metadata.json` must contain:

```json
{"population_status": "official_GBD_2023", "submission_ready": true}
```

before the manuscript is submitted.
