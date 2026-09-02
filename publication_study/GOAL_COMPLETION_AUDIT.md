# Publication-readiness completion audit

Audit date: 2 September 2026

This audit maps the 27-item research completion plan to current repository
evidence. A provisional engineering result is not counted as proof of a
production-data requirement. The authoritative mechanical gate remains
`verify_reproducibility.py --require-submission-ready`.

## Current conclusion

The fine-age GBD 2023 schizophrenia export is now present, prepared, and
validated, but the study is not submission-ready. The repository still lacks
the matching official GBD 2023 population export and its completed provenance
sidecar. Consequently, final APC and decomposition results, integrated
manuscript interpretation, production documents, and PDF visual QA remain
pending.

## Requirement status

| Item | Status | Authoritative evidence or remaining condition |
|---:|---|---|
| 1 | Complete | Descriptive BIC-selected segmented regression remains unchanged; inferential permutation/F-test fields are excluded; residual diagnostics and descriptive contrasts are retained; 39 tests pass. |
| 2 | Complete for current phase | Existing checked-in outputs remain explicitly provisional/obsolete. Clean audit builds are under `publication_study/build/`, not promoted as final outputs. |
| 3 | Complete | `ANALYSIS_PROTOCOL.md` defines the three layers, intentional APC retention, distinct estimands, post-data-access analysis specification, and point-estimate uncertainty limits. |
| 4 | Partly complete; blocked by population input | The preserved split burden ZIPs, deterministic preparation script, 17,544-row canonical burden CSV, and completed burden sidecar are present. All required fine-age and summary cells validate. The matching official population export remains absent. |
| 5 | Implemented; production execution pending | Official population loader, complete-key validation, proxy comparison, discrepancy table, and production command without `--allow-proxy-population` are implemented and pass a synthetic end-to-end test. Real-data rerun is pending item 4. |
| 6 | Complete | APC is isolated in `apc_analysis.py`; equal consecutive age intervals are enforced; aggregated/open-ended ages are excluded. |
| 7 | Complete | Primary 1994-2023 and sensitivity 1990-2019 six-period windows are implemented and tested. |
| 8 | Complete | Net drift, local drift, longitudinal age curve, period RR, and cohort RR are generated using an identifiable estimable-function parameterization. |
| 9 | Implemented; official-population confirmation pending | A fine-burden APC-only audit using exactly reconstructed 10-69 denominators found that the prior China-male incidence sign discrepancy resolved (segmented -0.0014%, APC -0.0170%), while the US-male discrepancy persisted (segmented -0.0756%, APC +0.0383%) and China-female estimates also straddled zero (segmented +0.0125%, APC -0.0172%). Both APC windows retained the same signs. These are provisional point-estimate checks, not final population-validated conclusions. |
| 10 | Complete | APC outputs are point estimates; no nominal Poisson/Wald inference is promoted as GBD posterior uncertainty. |
| 11 | Complete | Tests cover complete matrices, interval width, cohort indices, references, synthetic recovery, missing cells, and duplicate cells. |
| 12 | Implemented; production execution pending | Fine 20-group Shapley decomposition and age-specific-rate terminology pass the synthetic production test. Real estimates require item 4. |
| 13 | Complete | All-age reconstruction is a separately labeled QA check and the specification warns that closure does not validate age-bin granularity. |
| 14 | Complete in code; final result pending | Finest-versus-collapsed sensitivity reports sign, rank, and magnitude changes. Production conclusions require item 4. |
| 15 | Descriptive portion complete; population-dependent portion pending | The replacement summary panels match the frozen input to at most 4.66e-10 absolute (2.50e-16 relative). Segmented summaries and sensitivities are identical, and fitted values differ by at most 5.69e-14. Population-dependent reruns still require the official population. |
| 16 | Complete | Breakpoint-count, minimum-length, calendar-window, rate-scale, and UI-width-weighted sensitivities are generated. |
| 17 | Complete in code; final table pending | `cross_analysis_consistency.csv` contains the required trend, APC, and decomposition fields. Final values require item 4. |
| 18 | Complete provisionally; final audit pending | The fine-burden APC check explicitly preserves the US-male and China-female near-zero sign disagreements rather than reconciling them. Their different age coverage, age standardization/population weighting, global APC drift, and piecewise aggregate estimands remain the explanatory candidates defined in the analysis specification. Final tables and notes must be regenerated and reviewed with official population. |
| 19 | Pending item 4 | Final authoritative outputs must be generated in a new empty directory. Both build entry points now reject nonempty output directories. |
| 20 | Implemented provisionally; final rewrite pending | Methods distinguish trend, APC, and decomposition estimands and uncertainty. Final wording must be checked against production outputs. |
| 21 | Implemented provisionally; final rewrite pending | Current Results avoid significance language and causal decomposition claims. Effect-specific APC/decomposition prose must be finalized from production results. |
| 22 | Implemented provisionally; final strengthening pending | Discussion distinguishes aggregate, age-specific, and demographic interpretations and treats cohort curves as noncausal. Production-specific integration remains pending. |
| 23 | Complete in framework; final document check pending | Point-estimate, covariance, posterior-draw, non-independence, APC, trend, and decomposition limitations are stated. |
| 24 | Not triggered | Posterior draws are absent and are explicitly optional rather than a blocker. |
| 25 | Partly complete | Packages are pinned, paths are portable, both burden preparations are byte-reproducible, commands are documented, and burden export IDs, retrieval date, dimensions, raw hashes, and transformations are recorded and validated. Population provenance still requires the official export. |
| 26 | Partly complete | Current environment, preparation, 37 tests, analysis, APC, documents, `py_compile`, verifier, cross-table checks, terminology scan, `git diff --check`, and figure QA pass provisionally. Final fresh production run and rendered-document QA require items 4 and 19. |
| 27 | Not satisfied | `submission_ready` correctly remains false. Fine-age burden is validated; the official population and its sidecar remain the decisive data gate before final regeneration and review. |

## Validated engineering evidence

- `publication_study/tests`: 39 passing tests.
- Synthetic production-path test: official-population schema, fine-age burden,
  no optional YLD panel, APC ages 10-69, 20-group decomposition, provenance,
  production manuscript wording, and `submission_ready=true` all execute.
- Clean available-data build:
  `publication_study/build/audit_provisional_v4`.
- Available-data verifier: zero failures and one expected submission-gate
  warning.
- `cause_all.csv` preparation: byte-identical SHA-256 before and after rebuild.
- Fine-age burden preparation: deterministic SHA-256
  `2311f85973d75572bc1c812fef4f995e1da5bb600b736d2961aceecdc9b57689`;
  all three outcomes have complete 2,720-cell Number and Rate panels.
- Replacement descriptive comparison: all segmented summaries and
  specification sensitivities are unchanged; maximum fitted-value difference
  is numerical roundoff (`5.69e-14`).
- Proxy limitation discovered and retained as a hard error: both `0-4 years`
  and `5-9 years` have zero Number and Rate, so their separate populations
  cannot be inferred without fabricating an allocation.
- Provisional fine-burden APC audit: all 1,632 population cells for ages 10-69
  are exactly reconstructible; primary and sensitivity APC net-drift signs
  agree for all four location-sex incidence panels. US-male and China-female
  APC-versus-segmented signs differ at very small magnitudes and remain flagged.
- No obsolete post-data-access or rate-component wording in current source or
  newly generated documents.

## Exact unblock inputs

Place the following authenticated files under the documented `data/` paths:

1. matching official GBD 2023 population CSV;
2. completed population sidecar based on
   `gbd_population_export_metadata_template.json`; and
3. every unchanged raw population file referenced by that sidecar.

Afterward run the production commands in `README.md` into a new directory and
require `--require-submission-ready` to pass before final manuscript and PDF QA.
