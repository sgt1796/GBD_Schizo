"""Read-only integrity checks for the schizophrenia analysis results."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_RESULTS = HERE / "results"

REQUIRED_FILES = (
    "build_metadata.json",
    "data_dictionary.csv",
    "qa/validation_summary.json",
    "qa/methodological_notes.md",
    "tables/analysis_tables.xlsx",
    "tables/endpoint_summary.csv",
    "tables/segmented_summary.csv",
    "tables/segmented_ar1_sensitivity.csv",
    "tables/decomposition.csv",
    "tables/decomposition_path_sensitivity.csv",
    "tables/apc_descriptive_summary.csv",
    "tables/source_export_zero_audit.csv",
    "figures/main/figure_1_asr_trends.png",
    "figures/main/figure_2_segmented_trends.png",
    "figures/main/figure_3_age_patterns.png",
    "figures/main/figure_4_decomposition.png",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    failures: list[str] = []
    for relative in REQUIRED_FILES:
        if not (DEFAULT_RESULTS / relative).is_file():
            failures.append(f"missing result file: {relative}")

    metadata_path = DEFAULT_RESULTS / "build_metadata.json"
    qa_path = DEFAULT_RESULTS / "qa" / "validation_summary.json"
    if metadata_path.is_file() and qa_path.is_file():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
        if metadata.get("analysis_ready") is not True:
            failures.append("build metadata does not mark analysis_ready=true")
        if qa.get("analysis_ready") is not True:
            failures.append("QA summary does not mark analysis_ready=true")
        if qa.get("internal_validation_passed") is not True:
            failures.append("internal validation did not pass")
        if qa.get("all_primary_panels_complete_34_years") is not True:
            failures.append("annual panels are incomplete")
        if qa.get("all_age_count_reconstruction_within_tolerance") is not True:
            failures.append("all-age reconstruction is outside tolerance")
        if qa.get("maximum_absolute_decomposition_closure_error", 1) >= 1e-8:
            failures.append("decomposition closure error exceeds tolerance")
        for key in ("submission_ready", "scientific_review_gate_passed"):
            if key in metadata or key in qa:
                failures.append(f"obsolete workflow field remains: {key}")

        for field in ("burden_csv", "population_csv"):
            relative = metadata.get(field)
            expected = metadata.get(f"{field}_sha256")
            path = ROOT / relative if relative else None
            if not path or not path.is_file():
                failures.append(f"missing input referenced by {field}")
            elif expected != sha256(path):
                failures.append(f"SHA-256 mismatch for {field}")

    if failures:
        for failure in failures:
            print(f"FAIL  {failure}")
        print(f"\nSummary: {len(failures)} failure(s)")
        return 1

    print(f"PASS  {len(REQUIRED_FILES)} required result files")
    print("PASS  analysis metadata and numerical QA")
    print("PASS  canonical input SHA-256 fingerprints")
    print("\nSummary: 0 failure(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
