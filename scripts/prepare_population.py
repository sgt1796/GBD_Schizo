from __future__ import annotations

from itertools import product
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "GBD_data"
SOURCE = DATA_DIR / "IHME-GBD_2023_DATA-a9a792bb-1.zip"
OUTPUT = DATA_DIR / "GBD_2023_population_China_US.csv"

LOCATIONS = ("China", "United States of America")
SEXES = ("Female", "Male")
YEARS = tuple(range(1990, 2024))
AGES = (
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
KEYS = ("location_name", "sex_name", "age_name", "year")


def read_single_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing raw population export: {path}")
    with ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
        if len(members) != 1:
            raise ValueError(f"{path.name} must contain exactly one CSV; found {members}")
        return pd.read_csv(archive.open(members[0]), low_memory=False)


def validate_raw(frame: pd.DataFrame) -> pd.DataFrame:
    required = {*KEYS, "measure_name", "metric_name", "val", "lower", "upper"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Raw population export is missing columns: {sorted(missing)}")

    out = frame.copy()
    out["age_name"] = out["age_name"].replace({"<5 years": "0-4 years"})
    for column in ("year", "val", "lower", "upper"):
        out[column] = pd.to_numeric(out[column], errors="raise")
    out["year"] = out["year"].astype(int)

    if set(out.measure_name) != {"Population"}:
        raise ValueError(f"Expected only Population; found {sorted(set(out.measure_name))}")
    if set(out.metric_name) != {"Number"}:
        raise ValueError(f"Expected only Number; found {sorted(set(out.metric_name))}")
    if out.duplicated(list(KEYS)).any():
        raise ValueError("Raw population export contains duplicate dimensional keys.")

    expected = pd.MultiIndex.from_tuples(
        list(product(LOCATIONS, SEXES, AGES, YEARS)), names=KEYS
    )
    actual = pd.MultiIndex.from_frame(out[list(KEYS)])
    missing_cells = expected.difference(actual)
    unexpected_cells = actual.difference(expected)
    if len(missing_cells) or len(unexpected_cells):
        raise ValueError(
            "Population dimensions are incomplete: "
            f"missing={len(missing_cells)}, unexpected={len(unexpected_cells)}"
        )

    estimates = out[["lower", "val", "upper"]].to_numpy(dtype=float)
    if not np.isfinite(estimates).all() or (estimates <= 0).any():
        raise ValueError("Population estimates and bounds must be finite and positive.")
    if ((out.lower > out.val) | (out.val > out.upper)).any():
        raise ValueError("Population uncertainty bounds do not contain val.")
    return out


def main() -> None:
    raw = validate_raw(read_single_csv(SOURCE))
    canonical = raw[[*KEYS, "val"]].rename(columns={"val": "population"})
    canonical["gbd_release"] = "GBD 2023"
    canonical = canonical.sort_values(list(KEYS), kind="stable").reset_index(drop=True)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(OUTPUT, index=False, lineterminator="\n")
    print(f"Wrote {len(canonical):,} validated rows to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
