from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile


BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent
SOURCE_DIR = PROJECT_DIR / "schizo"
OUTPUT_DIR = BASE_DIR / "prepared_inputs"

CAUSE_OUTPUT = OUTPUT_DIR / "cause_all.csv"
POD_OUTPUT = OUTPUT_DIR / "GBD_1990_2023_ProbabilityOfDeath_ChinaUS_Schizophrenia.csv"
REPORT_OUTPUT = OUTPUT_DIR / "GBD_1990_2023_schizophrenia_preparation_report.txt"

CAUSE_SOURCES = [
    "IHME-GBD_2023_DATA-8a231d74-1.zip",
    "IHME-GBD_2023_DATA-d3e13215-1.zip",
]
POD_SOURCES = [
    "IHME-GBD_2023_DATA-8d0bd0df-1.zip",
]
EXCLUDED_RISK_SOURCES = [
    "IHME-GBD_2023_DATA-5ef7a575-1.zip",
]

CAUSE_COLUMNS = [
    "population_group_id",
    "population_group_name",
    "measure_id",
    "measure_name",
    "location_id",
    "location_name",
    "sex_id",
    "sex_name",
    "age_id",
    "age_name",
    "cause_id",
    "cause_name",
    "metric_id",
    "metric_name",
    "year",
    "val",
    "upper",
    "lower",
]

RISK_COLUMNS = CAUSE_COLUMNS[:12] + ["rei_id", "rei_name"] + CAUSE_COLUMNS[12:]

DIMENSION_COLUMNS = [
    "population_group_id",
    "measure_id",
    "location_id",
    "sex_id",
    "age_id",
    "cause_id",
    "metric_id",
    "year",
]


def numeric(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return 10**9


def row_sort_key(row: dict[str, str]) -> tuple[int, ...]:
    return (
        numeric(row["population_group_id"]),
        numeric(row["measure_id"]),
        numeric(row["location_id"]),
        numeric(row["sex_id"]),
        numeric(row["age_id"]),
        numeric(row["cause_id"]),
        numeric(row["metric_id"]),
        numeric(row["year"]),
    )


def read_zip_csv(zip_name: str) -> tuple[list[dict[str, str]], list[str], str]:
    zip_path = SOURCE_DIR / zip_name
    if not zip_path.exists():
        raise FileNotFoundError(f"Missing source zip: {zip_path}")

    with ZipFile(zip_path) as archive:
        csv_names = [name for name in archive.namelist() if name.lower().endswith(".csv")]
        if len(csv_names) != 1:
            raise ValueError(f"{zip_name} should contain exactly one CSV, found {csv_names}")

        csv_name = csv_names[0]
        with archive.open(csv_name) as raw:
            reader = csv.DictReader(line.decode("utf-8-sig") for line in raw)
            rows = list(reader)
            fieldnames = reader.fieldnames or []

    return rows, fieldnames, csv_name


def collect_sources(zip_names: list[str], expected_columns: list[str]) -> tuple[list[dict[str, str]], dict[str, int]]:
    all_rows: list[dict[str, str]] = []
    source_counts: dict[str, int] = {}

    for zip_name in zip_names:
        rows, fieldnames, _ = read_zip_csv(zip_name)
        if fieldnames != expected_columns:
            raise ValueError(
                f"Unexpected columns in {zip_name}: {fieldnames}. Expected: {expected_columns}"
            )
        source_counts[zip_name] = len(rows)
        all_rows.extend(rows)

    all_rows.sort(key=row_sort_key)
    return all_rows, source_counts


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def duplicate_count(rows: list[dict[str, str]], dimensions: list[str]) -> int:
    counts = Counter(tuple(row[column] for column in dimensions) for row in rows)
    return sum(count - 1 for count in counts.values() if count > 1)


def summarize(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    keys = [
        "measure_name",
        "location_name",
        "sex_name",
        "age_name",
        "cause_name",
        "metric_name",
        "year",
    ]
    result: dict[str, list[str]] = {}
    for key in keys:
        values = {row[key] for row in rows}
        if key == "year":
            result[key] = sorted(values, key=numeric)
        else:
            result[key] = sorted(values)
    return result


def summarize_excluded_risk(zip_name: str) -> dict[str, object]:
    rows, fieldnames, _ = read_zip_csv(zip_name)
    if fieldnames != RISK_COLUMNS:
        raise ValueError(f"Unexpected risk columns in {zip_name}: {fieldnames}")

    values: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        for column in [
            "measure_name",
            "location_name",
            "sex_name",
            "age_name",
            "cause_name",
            "rei_name",
            "metric_name",
            "year",
        ]:
            values[column].add(row[column])

    return {
        "rows": len(rows),
        "measures": sorted(values["measure_name"]),
        "locations": sorted(values["location_name"]),
        "sexes": sorted(values["sex_name"]),
        "ages": sorted(values["age_name"]),
        "cause": sorted(values["cause_name"]),
        "risks": sorted(values["rei_name"]),
        "metrics": sorted(values["metric_name"]),
        "years": sorted(values["year"], key=numeric),
    }


def format_list(values: list[str]) -> str:
    return ", ".join(values)


def write_report(
    cause_rows: list[dict[str, str]],
    pod_rows: list[dict[str, str]],
    cause_counts: dict[str, int],
    pod_counts: dict[str, int],
    excluded_risk: dict[str, object],
    cause_duplicates: int,
    pod_duplicates: int,
) -> None:
    cause_summary = summarize(cause_rows)
    pod_summary = summarize(pod_rows)
    years = cause_summary["year"]
    pod_years = pod_summary["year"]
    warnings: list[str] = []

    if years != pod_years:
        warnings.append("Cause and probability-of-death year ranges differ.")
    if cause_duplicates:
        warnings.append(f"Cause table has {cause_duplicates} duplicate dimensional keys.")
    if pod_duplicates:
        warnings.append(
            f"Probability-of-death table has {pod_duplicates} duplicate dimensional keys."
        )
    if set(cause_summary["sex_name"]) != {"Female", "Male"}:
        warnings.append("Cause table does not contain exactly Female and Male sex strata.")
    if set(pod_summary["sex_name"]) != {"Female", "Male"}:
        warnings.append(
            "Probability-of-death table does not contain exactly Female and Male sex strata."
        )

    lines = [
        "GBD 2023 schizophrenia input preparation report",
        "================================================",
        "",
        "Requested locations: China, United States of America",
        f"Requested years: {years[0]}-{years[-1]}",
        "Preparation method: combine only full-period GBD 2023 downloads; no historical/new-vintage splice.",
        "Sex strata: Female, Male",
        "Upper age bin in source data: 70+ years",
        "",
        "Cause of death or injury sources:",
    ]

    for name, count in cause_counts.items():
        lines.append(f"- {name}: {count:,} rows")
    lines.extend(
        [
            f"- combined: {len(cause_rows):,} rows; "
            f"{'no duplicate dimensional keys' if cause_duplicates == 0 else f'{cause_duplicates:,} duplicate dimensional keys'}; "
            f"years {years[0]}-{years[-1]}",
            f"- measures: {format_list(cause_summary['measure_name'])}",
            f"- ages: {format_list(cause_summary['age_name'])}",
            "",
            "Probability of death sources:",
        ]
    )

    for name, count in pod_counts.items():
        lines.append(f"- {name}: {count:,} rows")
    lines.extend(
        [
            f"- combined: {len(pod_rows):,} rows; "
            f"{'no duplicate dimensional keys' if pod_duplicates == 0 else f'{pod_duplicates:,} duplicate dimensional keys'}; "
            f"years {pod_years[0]}-{pod_years[-1]}",
            f"- ages: {format_list(pod_summary['age_name'])}",
            "",
            "Risk-factor sources intentionally omitted:",
        ]
    )

    for zip_name in EXCLUDED_RISK_SOURCES:
        lines.append(
            f"- {zip_name}: {excluded_risk['rows']:,} rows; "
            f"rei labels in source: {format_list(excluded_risk['risks'])}; "
            "not written because these labels represent a single sexual-violence risk branch "
            "rather than a multi-risk analytic table."
        )

    lines.extend(
        [
            "",
            "Prepared tables:",
            f"- {CAUSE_OUTPUT.name}: {len(cause_rows):,} rows",
            f"- {POD_OUTPUT.name}: {len(pod_rows):,} rows",
            "",
            "Consumer compatibility:",
            "- Output CSV columns match the corresponding breast cancer prepared input schemas.",
            "- The schizophrenia cause table keeps both Female and Male strata.",
            "- No schizophrenia risk-factor table was generated.",
            "",
            "Validation warnings: " + ("; ".join(warnings) if warnings else "none"),
        ]
    )

    REPORT_OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    cause_rows, cause_counts = collect_sources(CAUSE_SOURCES, CAUSE_COLUMNS)
    pod_rows, pod_counts = collect_sources(POD_SOURCES, CAUSE_COLUMNS)
    excluded_risk = summarize_excluded_risk(EXCLUDED_RISK_SOURCES[0])

    cause_duplicates = duplicate_count(cause_rows, DIMENSION_COLUMNS)
    pod_duplicates = duplicate_count(pod_rows, DIMENSION_COLUMNS)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(CAUSE_OUTPUT, cause_rows, CAUSE_COLUMNS)
    write_csv(POD_OUTPUT, pod_rows, CAUSE_COLUMNS)
    write_report(
        cause_rows,
        pod_rows,
        cause_counts,
        pod_counts,
        excluded_risk,
        cause_duplicates,
        pod_duplicates,
    )

    print(f"Wrote {CAUSE_OUTPUT} ({len(cause_rows):,} rows)")
    print(f"Wrote {POD_OUTPUT} ({len(pod_rows):,} rows)")
    print(f"Wrote {REPORT_OUTPUT}")


if __name__ == "__main__":
    main()
