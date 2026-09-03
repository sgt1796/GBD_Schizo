"""Read-only checks for the publication environment and generated package."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent

ANALYSIS_FILES = (
    "build_metadata.json",
    "data_dictionary.csv",
    "qa/validation_summary.json",
    "tables/publication_tables.xlsx",
    "tables/endpoint_summary.csv",
    "tables/all_age_count_reconstruction.csv",
    "tables/segmented_summary.csv",
    "tables/segmented_segments.csv",
    "tables/segmented_specification_sensitivity.csv",
    "tables/trajectory_contrasts.csv",
    "tables/trend_excluding_2020_2023.csv",
    "tables/decomposition.csv",
    "tables/decomposition_age_bin_sensitivity.csv",
    "tables/annual_chained_decomposition.csv",
    "tables/fiveyear_chained_decomposition.csv",
    "tables/apc_summary.csv",
    "tables/apc_local_drift.csv",
    "tables/apc_age_curve.csv",
    "tables/apc_period_rr.csv",
    "tables/apc_cohort_rr.csv",
    "tables/apc_sensitivity_summary_1990_2019.csv",
    "tables/apc_window_sensitivity.csv",
    "tables/apc_primary_direction_agreement.csv",
    "tables/cross_analysis_consistency.csv",
    "tables/cross_method_contradictions.csv",
    "qa/methodological_notes.md",
    "figures/main/figure_1_asr_trends.png",
    "figures/main/figure_2_segmented_trends.png",
    "figures/main/figure_3_age_patterns.png",
    "figures/main/figure_4_decomposition.png",
    "figures/supplement/figure_s1_counts.png",
    "figures/supplement/figure_s2_apc_estimable_functions.png",
    "nci_joinpoint_inputs/input_manifest.csv",
    "nci_joinpoint_inputs/analysis_settings.json",
)

DOCUMENT_FILES = (
    "documents/manuscript_BMC_Public_Health.docx",
    "documents/supplementary_material.docx",
    "documents/statistical_methods_appendix.docx",
    "documents/GATHER_checklist.docx",
    "documents/document_manifest.json",
)

STALE_ANALYSIS_FILES = (
    "tables/pairwise_parallelism.csv",
    "tables/prepandemic_trend_sensitivity.csv",
    "tables/apc_excluding_2020_2023.csv",
    "tables/apc_excluding_2019_2023.csv",
    "tables/apc_period_curvature.csv",
    "tables/apc_cohort_curvature.csv",
)

ALL_AGE_SCHEMA = {
    "tables/data_audit.csv": {
        "all_age_groups",
        "all_age_year_cells",
    },
    "tables/all_age_count_reconstruction.csv": {
        "reconstructed_all_age_count",
        "reported_all_age_count",
        "relative_error_pct",
        "within_tolerance",
        "tolerance_pct",
    },
    "tables/decomposition.csv": {
        "all_age_count_start_reconstructed",
        "all_age_count_end_reconstructed",
        "population_size_change",
        "age_structure_change",
        "age_specific_rate_change",
    },
    "tables/annual_chained_decomposition.csv": {
        "all_age_count_start_reconstructed",
        "all_age_count_end_reconstructed",
    },
    "tables/fiveyear_chained_decomposition.csv": {
        "all_age_count_start_reconstructed",
        "all_age_count_end_reconstructed",
    },
}

PRIMARY_DESCRIPTIVE_TABLES = (
    "tables/segmented_summary.csv",
    "tables/segmented_segments.csv",
    "tables/trajectory_contrasts.csv",
    "tables/apc_summary.csv",
    "tables/apc_local_drift.csv",
    "tables/apc_age_curve.csv",
    "tables/apc_period_rr.csv",
    "tables/apc_cohort_rr.csv",
)

LEGACY_AGE_FIELDS = {
    "adult_age_groups",
    "adult_age_year_cells",
    "adult_count_start_reconstructed",
    "adult_count_end_reconstructed",
    "age_15_plus_groups",
    "age_15_plus_year_cells",
    "age_15_plus_count_start_reconstructed",
    "age_15_plus_count_end_reconstructed",
}

LEGACY_METADATA_FIELDS = {
    "permutations",
    "seed",
    "legacy_permutations_argument_ignored",
    "legacy_seed_argument_ignored",
}


class Results:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.warnings: list[str] = []

    def passed(self, message: str) -> None:
        print(f"PASS  {message}")

    def info(self, message: str) -> None:
        print(f"INFO  {message}")

    def warn(self, message: str) -> None:
        self.warnings.append(message)
        print(f"WARN  {message}")

    def fail(self, message: str) -> None:
        self.failures.append(message)
        print(f"FAIL  {message}")


def load_json(path: Path, results: Results) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        results.fail(f"cannot read valid JSON from {path}: {exc}")
        return {}


def pinned_requirements(path: Path, results: Results) -> list[tuple[str, str]]:
    pins: list[tuple[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.count("==") != 1:
            results.fail(f"requirement is not exactly pinned: {line}")
            continue
        name, expected = (part.strip() for part in line.split("==", 1))
        pins.append((name, expected))
    return pins


def check_environment(results: Results) -> None:
    if sys.version_info[:2] == (3, 12):
        results.passed(f"CPython {sys.version.split()[0]}")
    else:
        results.fail(
            f"CPython 3.12 is required for the recorded environment; found {sys.version.split()[0]}"
        )

    requirements = HERE / "requirements.txt"
    if not requirements.is_file():
        results.fail(f"missing {requirements}")
        return

    for name, expected in pinned_requirements(requirements, results):
        try:
            actual = version(name)
        except PackageNotFoundError:
            results.fail(f"missing package {name}=={expected}")
            continue
        if actual == expected:
            results.passed(f"{name}=={actual}")
        else:
            results.fail(f"{name}: expected {expected}, found {actual}")

    for relative in (
        "prepare_schizo_inputs.py",
        "schizo/IHME-GBD_2023_DATA-8a231d74-1.zip",
        "schizo/IHME-GBD_2023_DATA-d3e13215-1.zip",
        "schizo/IHME-GBD_2023_DATA-8d0bd0df-1.zip",
        "prepared_inputs/cause_all.csv",
        "prepared_inputs/GBD_1990_2023_ProbabilityOfDeath_ChinaUS_Schizophrenia.csv",
        "prepared_inputs/GBD_1990_2023_schizophrenia_preparation_report.txt",
        "publication_study/publication_analysis.py",
        "publication_study/apc_analysis.py",
        "publication_study/build_documents.py",
        "publication_study/prepare_production_burden.py",
        "publication_study/prepare_production_population.py",
        "publication_study/population_input_template.csv",
        "publication_study/gbd_export_metadata_template.json",
        "publication_study/gbd_population_export_metadata_template.json",
        "publication_study/nci_results_template.csv",
        "addition_inputs/IHME-GBD_2023_DATA-774041bd-1.zip",
        "addition_inputs/IHME-GBD_2023_DATA-22ef74c2-1.zip",
        "addition_inputs/IHME-GBD_2023_DATA-a9a792bb-1.zip",
        "data/GBD_2023_schizophrenia_fine_age_China_US.csv",
        "data/GBD_2023_population_China_US.csv",
        "data/metadata/burden_export.json",
        "data/metadata/population_export.json",
    ):
        path = REPO_ROOT / relative
        if path.is_file() and path.stat().st_size > 0:
            results.passed(f"source file {relative}")
        else:
            results.fail(f"missing or empty source file {relative}")

    excluded_risk = REPO_ROOT / "schizo" / "IHME-GBD_2023_DATA-5ef7a575-1.zip"
    if excluded_risk.is_file() and excluded_risk.stat().st_size > 0:
        results.info("optional excluded risk archive is present for provenance reporting")
    else:
        results.info("optional excluded risk archive is absent (not an analysis input)")


def require_files(base: Path, names: tuple[str, ...], results: Results) -> None:
    missing = [
        name
        for name in names
        if not (base / name).is_file() or (base / name).stat().st_size == 0
    ]
    if missing:
        for name in missing:
            results.fail(f"missing generated file {base / name}")
    else:
        results.passed(f"{len(names)} required generated files")


def csv_columns(path: Path, results: Results) -> set[str]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return set(next(csv.reader(handle)))
    except (OSError, StopIteration, csv.Error) as exc:
        results.fail(f"cannot read CSV header from {path}: {exc}")
        return set()


def is_inferential_result_field(column: str) -> bool:
    name = column.lower()
    return (
        "p_value" in name
        or "q_value" in name
        or "model_ci" in name
        or name.endswith("_lower_ci")
        or name.endswith("_upper_ci")
        or name in {"f_statistic", "df1", "df2", "parallel_trends_at_q_0_05"}
    )


def check_current_table_schema(base: Path, results: Results) -> None:
    stale = [name for name in STALE_ANALYSIS_FILES if (base / name).exists()]
    if stale:
        for name in stale:
            results.fail(f"legacy generated file is present; rebuild in a fresh directory: {base / name}")
    else:
        results.passed("no legacy permutation-era table files")

    for name, required in ALL_AGE_SCHEMA.items():
        columns = csv_columns(base / name, results)
        missing = required - columns
        legacy = LEGACY_AGE_FIELDS & columns
        if missing:
            results.fail(f"{name} is missing all-age columns: {sorted(missing)}")
        if legacy:
            results.fail(f"{name} contains legacy adult columns: {sorted(legacy)}")
        if columns and not missing and not legacy:
            results.passed(f"{name} uses all_age field names")

    for name in PRIMARY_DESCRIPTIVE_TABLES:
        columns = csv_columns(base / name, results)
        inferential = sorted(column for column in columns if is_inferential_result_field(column))
        if inferential:
            results.fail(f"{name} contains inferential result fields: {inferential}")
        elif columns:
            results.passed(f"{name} contains descriptive fields only")


def check_apc_completeness(base: Path, results: Results) -> None:
    expected_outcomes = {"Incidence", "Prevalence", "DALYs"}
    expected_panels = 12
    paths = {
        "APC summary": base / "tables" / "apc_summary.csv",
        "APC window sensitivity": base / "tables" / "apc_window_sensitivity.csv",
        "cross-analysis consistency": base / "tables" / "cross_analysis_consistency.csv",
    }
    loaded: dict[str, list[dict[str, str]]] = {}
    for label, path in paths.items():
        try:
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                loaded[label] = list(csv.DictReader(handle))
        except (OSError, csv.Error) as exc:
            results.fail(f"cannot inspect {label}: {exc}")
            return
    for label, rows in loaded.items():
        panels = {
            (row.get("location_name"), row.get("sex_name"), row.get("measure_name"))
            for row in rows
        }
        outcomes = {row.get("measure_name") for row in rows}
        if len(rows) == expected_panels and len(panels) == expected_panels and outcomes == expected_outcomes:
            results.passed(f"{label} covers all 12 country-sex-outcome panels")
        else:
            results.fail(
                f"{label} is incomplete: rows={len(rows)}, unique panels={len(panels)}, outcomes={sorted(str(value) for value in outcomes)}"
            )
    consistency = loaded["cross-analysis consistency"]
    drifts = []
    for row in consistency:
        try:
            drifts.append(float(row.get("apc_net_drift_1994_2023", "")))
        except (TypeError, ValueError):
            drifts.append(float("nan"))
    if len(drifts) == expected_panels and all(math.isfinite(value) for value in drifts):
        results.passed("cross-analysis consistency has finite APC net drift for every outcome")
    else:
        results.fail("cross-analysis consistency has missing or non-finite APC net drift")


def check_qa(qa: dict, results: Results) -> None:
    expected_true = (
        "all_primary_panels_complete_34_years",
        "all_age_count_reconstruction_within_tolerance",
        "yld_daly_audit_passed",
        "primary_outputs_exclude_percent_metric",
        "internal_validation_passed",
    )
    expected_zero = (
        "duplicate_dimensional_keys",
        "invalid_ui_rows",
        "negative_rows",
    )
    for key in expected_true:
        if qa.get(key) is True:
            results.passed(f"QA {key}=true")
        else:
            results.fail(f"QA {key} is not true")
    for key in expected_zero:
        if qa.get(key) == 0:
            results.passed(f"QA {key}=0")
        else:
            results.fail(f"QA {key} is not zero")

    yld_status = qa.get("yld_daly_audit_status")
    if yld_status in {"verified_identical", "not_available"}:
        results.passed(f"QA yld_daly_audit_status={yld_status}")
    else:
        results.fail(
            "QA yld_daly_audit_status must be verified_identical or not_available; "
            f"found {yld_status}"
        )

    if qa.get("formal_trend_inference_performed") is False:
        results.passed("QA formal_trend_inference_performed=false")
    else:
        results.fail("QA formal_trend_inference_performed is not false")

    closure = qa.get("maximum_absolute_decomposition_closure_error")
    if isinstance(closure, (int, float)) and math.isfinite(closure) and closure <= 1e-6:
        results.passed(f"QA decomposition closure error={closure:.3g}")
    else:
        results.fail(f"QA decomposition closure error is missing, non-finite, or >1e-6: {closure}")

    reconstruction_rows = qa.get("all_age_count_reconstruction_rows")
    if reconstruction_rows == 408:
        results.passed("QA all-age count reconstruction has 408 panels")
    else:
        results.fail(f"QA all-age count reconstruction rows: expected 408, found {reconstruction_rows}")

    reconstruction_error = qa.get("all_age_count_reconstruction_max_absolute_relative_error_pct")
    reconstruction_tolerance = qa.get("all_age_count_reconstruction_tolerance_pct")
    reconstruction_numbers = (reconstruction_error, reconstruction_tolerance)
    if (
        all(isinstance(value, (int, float)) and math.isfinite(value) for value in reconstruction_numbers)
        and reconstruction_error <= reconstruction_tolerance
    ):
        results.passed(
            "QA all-age reconstruction error "
            f"{reconstruction_error:.3g}% <= {reconstruction_tolerance:.3g}% tolerance"
        )
    else:
        results.fail(
            "QA all-age reconstruction error is missing, non-finite, or exceeds tolerance: "
            f"error={reconstruction_error}, tolerance={reconstruction_tolerance}"
        )


def check_analysis(base: Path, require_ready: bool, results: Results) -> None:
    require_files(base, ANALYSIS_FILES, results)
    require_files(base, DOCUMENT_FILES, results)
    check_current_table_schema(base, results)
    check_apc_completeness(base, results)

    nci_series = list((base / "nci_joinpoint_inputs").glob("*.csv"))
    nci_series = [path for path in nci_series if path.name != "input_manifest.csv"]
    if len(nci_series) == 12:
        results.passed("12 generated NCI input series")
    else:
        results.fail(f"expected 12 generated NCI input series, found {len(nci_series)}")

    metadata = load_json(base / "build_metadata.json", results)
    qa = load_json(base / "qa" / "validation_summary.json", results)
    if not metadata or not qa:
        return
    check_qa(qa, results)

    official_population = metadata.get("population_status") == "official_GBD_2023"
    fine_age_burden = metadata.get("fine_age_burden_validated") is True
    source_metadata = metadata.get("source_metadata_complete") is True
    official_nci = qa.get("official_nci_results_imported") is True
    internal_validation = qa.get("internal_validation_passed") is True
    ready = metadata.get("submission_ready") is True
    legacy_metadata = sorted(LEGACY_METADATA_FIELDS & metadata.keys())
    if legacy_metadata:
        results.fail(f"metadata contains removed CLI fields: {legacy_metadata}")
    if metadata.get("formal_trend_inference_performed") is not False:
        results.fail("metadata formal_trend_inference_performed is not false")
    if metadata.get("official_nci_results_optional") is not True:
        results.fail("metadata does not identify official NCI output as optional")
    if qa.get("population_is_official_gbd_2023") is not official_population:
        results.fail("population provenance disagrees between metadata and QA")
    if qa.get("submission_ready") is not ready:
        results.fail("submission status disagrees between metadata and QA")
    if qa.get("fine_age_burden_validated") is not fine_age_burden:
        results.fail("fine-age burden status disagrees between metadata and QA")
    if qa.get("source_metadata_complete") is not source_metadata:
        results.fail("source-metadata status disagrees between metadata and QA")
    if ready != (
        official_population and fine_age_burden and source_metadata and internal_validation
    ):
        results.fail(
            "submission_ready is inconsistent with population provenance, fine-age input, or internal QA"
        )
    elif ready:
        results.passed("official-population gate supplied; submission_ready=true")
    else:
        unresolved = []
        if not official_population:
            unresolved.append("official GBD 2023 population is absent")
        if not fine_age_burden:
            unresolved.append("matching fine-age burden panels are absent")
        if not source_metadata:
            unresolved.append("complete export metadata sidecars are absent")
        if not internal_validation:
            unresolved.append("internal QA did not pass")
        message = "submission gate remains closed: " + "; ".join(unresolved)
        if require_ready:
            results.fail(message)
        else:
            results.warn(message)

    if official_nci:
        results.passed("optional official NCI validation imported")
    else:
        results.info("optional official NCI validation not imported (not a submission gate)")

    rendered_pdfs = sorted((base / "rendered").rglob("*.pdf")) if (base / "rendered").is_dir() else []
    if rendered_pdfs:
        results.warn(
            "rendered PDFs are present but are not validated by this command; regenerate them "
            "from the current DOCX files and complete page-by-page visual QA"
        )
    else:
        results.info(
            "PDF rendering is external to this automated verifier; consult "
            "GOAL_COMPLETION_AUDIT.md for the recorded page-by-page QA."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the pinned environment and a generated publication package."
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        help="Generated analysis directory. Omit to check only the environment and source files.",
    )
    parser.add_argument(
        "--require-submission-ready",
        action="store_true",
        help="Fail unless the official matching GBD 2023 population is recorded.",
    )
    args = parser.parse_args()
    if args.require_submission_ready and args.analysis_dir is None:
        parser.error("--require-submission-ready requires --analysis-dir")
    return args


def main() -> int:
    args = parse_args()
    results = Results()
    check_environment(results)
    if args.analysis_dir is not None:
        check_analysis(args.analysis_dir.resolve(), args.require_submission_ready, results)

    print()
    print(f"Summary: {len(results.failures)} failure(s), {len(results.warnings)} warning(s)")
    return 1 if results.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
