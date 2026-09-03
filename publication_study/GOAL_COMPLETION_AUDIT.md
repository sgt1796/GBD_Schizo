# Publication-readiness completion audit

Audit date: 3 September 2026

This audit maps the 27-item completion plan to the authoritative package in
`publication_study/output`. The reproducible computational gate passes with
`submission_ready=true`. The four DOCX files were rendered with LibreOffice and
all 30 PDF pages passed page-by-page layout inspection.

## Current conclusion

The required fine-age schizophrenia burden and official GBD 2023 population
exports are present, preserved, documented, and validated. All population-
dependent analyses were rerun from an empty output directory. APC now covers
incidence, prevalence, and DALYs in every country-sex panel. The tracked legacy
package was replaced, so retired permutation-era tables are absent from the
authoritative output.

The analysis, tables, figures, DOCX files, and rendered QA PDFs pass automated
reproducibility, numerical, schema, structural, and visual checks. A clean Python
3.12 environment reproduced the canonical inputs, passed all 41 tests, rebuilt
the full package, and produced 61 core files byte-identical to the authoritative
output. Authors still need to supply repository, funding, and contributor details
where applicable; those administrative metadata are not analytic findings.

## Requirement status

| Item | Status | Evidence |
|---:|---|---|
| 1 | Complete | Descriptive BIC-selected segmented regression is retained; inferential permutation/F-test output is absent; residual diagnostics and descriptive contrasts remain. |
| 2 | Complete | Provisional outputs were not polished before the final inputs; the final package was rebuilt only after the official population was integrated. |
| 3 | Complete | `ANALYSIS_PROTOCOL.md` defines the three analysis layers, distinct estimands, post-data-access analysis specification, and point-estimate limitations. |
| 4 | Complete | Preserved raw burden and population ZIPs, exact query sidecars, and canonical fine-age CSVs are present and hash-validated. |
| 5 | Complete | Official population is loaded and key-validated; comparison against reconstructible count-rate denominators shows numerical agreement, with unavailable zero/zero cells explicitly recorded. The production command does not use proxy population. |
| 6 | Complete | APC is isolated in `apc_analysis.py`; core ages are consecutive five-year groups from 10-14 through 65-69. |
| 7 | Complete | Primary 1994-2023 and sensitivity 1990-2019 windows each contain six complete five-year periods. |
| 8 | Complete | Net drift, local drift, longitudinal age curve, period RR, and cohort RR use an identifiable estimable-function parameterization. |
| 9 | Complete | APC is interpreted as complementary age-specific structure. China-male incidence now agrees in direction; US-male incidence remains discordant and is retained in the contradiction audit. |
| 10 | Complete | APC is based on central estimates and contains no nominal confidence intervals or significance claims. |
| 11 | Complete | Synthetic recovery, complete matrix, interval-width, cohort-index, reference, missing-cell, duplicate-cell, and multi-outcome isolation tests pass. |
| 12 | Complete | Fine 20-group Shapley decomposition uses official population and reports population-size, age-structure, and age-specific-rate components. |
| 13 | Complete | All-age closure is reported as QA only; the generated methodological note explicitly rejects granularity validation from closure. |
| 14 | Complete | Finest-versus-collapsed decomposition sensitivity reports sign, magnitude, and component-rank changes; three panels are flagged. |
| 15 | Complete | Descriptive and population-dependent analyses were rerun after replacement; count reconstruction and trend results pass numerical checks. |
| 16 | Complete | Breakpoint-count, minimum-length, window, rate-scale, and UI-width-weighted sensitivities are generated. |
| 17 | Complete | `cross_analysis_consistency.csv` has 12 unique country-sex-outcome rows and finite APC, trend, and decomposition fields for every row. |
| 18 | Complete | Ten panels with directional disagreement and/or opposing local drifts are explained in `qa/methodological_notes.md`; no implementation failure is inferred. |
| 19 | Complete | `publication_study/output` was regenerated from an empty directory; verifier confirms that legacy permutation-era files are absent. |
| 20 | Complete | Manuscript Methods distinguish descriptive trends, APC estimable functions, and demographic decomposition, including each method's uncertainty limits. |
| 21 | Complete | Results lead with magnitude and direction, report APC heterogeneity and decomposition components, and avoid significance-driven or causal accounting language. |
| 22 | Complete | Discussion integrates count, standardized-rate, sex, age, APC, and demographic patterns without causal policy claims. |
| 23 | Complete | Limitations state modeled-estimate, covariance, posterior-draw, method-comparison, and outcome-dependence constraints. |
| 24 | Not triggered | Posterior draws were unavailable and are documented as an optional future enhancement rather than a submission gate. |
| 25 | Complete | Export parameters, retrieval dates, hashes, portable paths, pinned versions, preparation commands, and deterministic canonical inputs are recorded. |
| 26 | Complete | A fresh Python 3.12 environment reproduced both canonical inputs, passed 41 tests, rebuilt analyses and documents, passed `py_compile` and the verifier, and reproduced 61 core files byte-for-byte. Numerical, terminology, schema, structural, six-figure, and 30-page PDF visual checks pass. |
| 27 | Complete | Official population, fine-age data, APC validation, decomposition, inconsistency notes, clean outputs, matching manuscript/supplement, provenance, and rendered-document QA are present; `submission_ready=true`. |

## Validated production evidence

- Authoritative output: `publication_study/output`.
- Full test suite: 41 passing tests.
- Reproducibility verifier: 0 failures and 0 warnings.
- Burden canonical SHA-256:
  `2311f85973d75572bc1c812fef4f995e1da5bb600b736d2961aceecdc9b57689`.
- Population canonical SHA-256:
  `5b67a430816174e048f71a7a1442061801253558fca40866626df7b464020c1e`.
- Population input: 2,720 rows, 20 ages, 34 years, two locations, and two sexes.
- APC: 12 complete primary panels and 12 complete sensitivity panels; directions
  agree in 11 of 12 panels. The sole change is China-female DALYs at magnitudes
  close to zero (-0.001645%/year versus +0.001900%/year).
- All-age reconstruction: 408 rows; maximum absolute relative error
  `2.41e-11%`.
- Decomposition closure: maximum absolute error `1.16e-10`.
- DOCX structural QA: all four files are valid OOXML archives, contain all
  expected sections, have no empty tables, and contain no obsolete APC wording.
- PDF QA: LibreOffice rendered 30 pages (14 manuscript, 12 supplement, two
  statistical appendix, and two GATHER checklist pages); no clipping, overflow,
  out-of-page text, blank pages, or unreadable audit tables remained after review.
- Fresh-environment reproducibility: 41 tests and the submission verifier pass;
  61 CSV, JSON, Markdown, and PNG outputs match the authoritative package by
  SHA-256.

## Author-supplied journal-upload metadata

1. Add the persistent code/data repository link after deposit.
2. Confirm funding and CRediT contributor statements in the unblinded submission.
