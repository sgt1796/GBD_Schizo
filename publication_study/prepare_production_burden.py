from __future__ import annotations

import hashlib
from itertools import product
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "addition_inputs"
OUTPUT = ROOT / "data" / "GBD_2023_schizophrenia_fine_age_China_US.csv"

BASE_EXPORT = SOURCE_DIR / "IHME-GBD_2023_DATA-774041bd-1.zip"
CORRECTION_EXPORT = SOURCE_DIR / "IHME-GBD_2023_DATA-22ef74c2-1.zip"

LOCATIONS = ("China", "United States of America")
SEXES = ("Female", "Male")
YEARS = tuple(range(1990, 2024))
MEASURES = ("Incidence", "Prevalence", "DALYs")
METRICS = ("Number", "Rate")
FINE_AGES = (
    "0-4 years",
    "5-9 years",
    "10-14 years",
    "15-19 years",
    "20-24 years",
    "25-29 years",
    "30-34 years",
    "35-39 years",
    "40-44 years",
    "45-49 years",
    "50-54 years",
    "55-59 years",
    "60-64 years",
    "65-69 years",
    "70-74 years",
    "75-79 years",
    "80-84 years",
    "85-89 years",
    "90-94 years",
    "95+ years",
)

KEYS = (
    "location_name",
    "sex_name",
    "age_name",
    "measure_name",
    "metric_name",
    "year",
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_single_csv(zip_path: Path) -> pd.DataFrame:
    if not zip_path.is_file():
        raise FileNotFoundError(f"Missing raw export: {zip_path}")
    with ZipFile(zip_path) as archive:
        members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
        if len(members) != 1:
            raise ValueError(f"{zip_path.name} must contain exactly one CSV; found {members}")
        return pd.read_csv(archive.open(members[0]), low_memory=False)


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    required_columns = {
        *KEYS,
        "cause_name",
        "val",
        "lower",
        "upper",
    }
    missing = required_columns - set(frame.columns)
    if missing:
        raise ValueError(f"Raw burden export is missing columns: {sorted(missing)}")

    out = frame.copy()
    out["measure_name"] = out["measure_name"].replace(
        {"DALYs (Disability-Adjusted Life Years)": "DALYs"}
    )
    out["age_name"] = out["age_name"].replace({"<5 years": "0-4 years"})
    for column in ("year", "val", "lower", "upper"):
        out[column] = pd.to_numeric(out[column], errors="raise")
    out["year"] = out["year"].astype(int)
    return out


def required_index() -> pd.MultiIndex:
    records = list(product(LOCATIONS, SEXES, FINE_AGES, MEASURES, METRICS, YEARS))
    # Retain the exported all-age crude Rate as well as Number. It is not a
    # production requirement, but it permits an auditable provisional total-
    # population reconstruction when a zero-burden age has Number=Rate=0.
    records.extend(product(LOCATIONS, SEXES, ("All ages",), MEASURES, METRICS, YEARS))
    records.extend(
        product(LOCATIONS, SEXES, ("Age-standardized",), MEASURES, ("Rate",), YEARS)
    )
    return pd.MultiIndex.from_tuples(records, names=KEYS)


def validate(frame: pd.DataFrame) -> None:
    causes = set(frame["cause_name"].dropna().astype(str).str.strip())
    if causes != {"Schizophrenia"}:
        raise ValueError(f"Expected only Schizophrenia; found {sorted(causes)}")
    if frame.duplicated(list(KEYS)).any():
        raise ValueError("Canonical burden data contain duplicate dimensional keys.")

    expected = required_index()
    actual = pd.MultiIndex.from_frame(frame[list(KEYS)])
    missing = expected.difference(actual)
    unexpected = actual.difference(expected)
    if len(missing) or len(unexpected):
        raise ValueError(
            "Canonical burden dimensions are incomplete: "
            f"missing={len(missing)}, unexpected={len(unexpected)}"
        )

    estimates = frame[["lower", "val", "upper"]]
    if not np.isfinite(estimates.to_numpy(dtype=float)).all():
        raise ValueError("Canonical burden estimates must be finite.")
    if (estimates < 0).any().any():
        raise ValueError("Canonical burden estimates must be nonnegative.")
    if ((frame["lower"] > frame["val"]) | (frame["val"] > frame["upper"])).any():
        raise ValueError("Canonical burden uncertainty bounds do not contain val.")


def main() -> None:
    base = canonicalize(read_single_csv(BASE_EXPORT))
    correction = canonicalize(read_single_csv(CORRECTION_EXPORT))

    correction_ages = set(correction["age_name"])
    if correction_ages != {"70-74 years"}:
        raise ValueError(
            f"Correction export must contain only 70-74 years; found {sorted(correction_ages)}"
        )

    # The base query accidentally selected the overlapping 65-74 aggregate.
    # Remove it, append the corrected 70-74 export, and retain only analysis inputs.
    merged = pd.concat(
        [base.loc[base["age_name"].ne("65-74 years")], correction],
        ignore_index=True,
    )
    required_cell = (
        merged["location_name"].isin(LOCATIONS)
        & merged["sex_name"].isin(SEXES)
        & merged["year"].isin(YEARS)
        & merged["measure_name"].isin(MEASURES)
        & merged["metric_name"].isin(METRICS)
        & (
            merged["age_name"].isin(FINE_AGES)
            | merged["age_name"].eq("All ages")
            | (
                merged["age_name"].eq("Age-standardized")
                & merged["metric_name"].eq("Rate")
            )
        )
    )
    canonical = merged.loc[required_cell].copy()
    canonical = canonical.sort_values(list(KEYS), kind="stable").reset_index(drop=True)
    validate(canonical)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(OUTPUT, index=False, lineterminator="\n")
    print(f"Wrote {len(canonical):,} validated rows to {OUTPUT.relative_to(ROOT)}")
    print(f"Output SHA-256: {file_sha256(OUTPUT)}")
    print(f"Base ZIP SHA-256: {file_sha256(BASE_EXPORT)}")
    print(f"Correction ZIP SHA-256: {file_sha256(CORRECTION_EXPORT)}")


if __name__ == "__main__":
    main()
