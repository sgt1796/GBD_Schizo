"""Prepare GBD incidence cells for the pinned NCI APC R implementation."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


AGE_LABELS = tuple(f"{age}-{age + 4} years" for age in range(10, 70, 5))
WINDOWS = {"1994_2023": 1994, "1990_2019": 1990}
RESULT_NAMES = ("summary", "age_curve", "local_drift", "period_rr", "cohort_rr")


def _rscript() -> str:
    explicit = os.environ.get("RSCRIPT_PATH")
    if explicit:
        if not Path(explicit).is_file():
            raise FileNotFoundError(f"RSCRIPT_PATH does not exist: {explicit}")
        return explicit
    found = shutil.which("Rscript")
    if found:
        return found
    raise FileNotFoundError("Rscript is required for APC; set RSCRIPT_PATH or add it to PATH")


def run_nci_apc(
    burden: pd.DataFrame, population: pd.DataFrame, tables_dir: Path,
    locations: tuple[str, ...], sexes: tuple[str, ...],
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    keys = ["location_name", "sex_name", "age_name", "year"]
    cases = burden[
        burden.measure_name.eq("Incidence")
        & burden.metric_name.eq("Number")
        & burden.age_name.isin(AGE_LABELS)
    ][keys + ["val"]].rename(columns={"val": "events"})
    pop = population[population.age_name.isin(AGE_LABELS)][keys + ["population"]]
    cells = cases.merge(pop, on=keys, how="outer", validate="one_to_one", indicator=True)
    if not cells._merge.eq("both").all():
        raise ValueError("APC incidence counts and population are not matched")
    if not np.isfinite(cells[["events", "population"]].to_numpy()).all():
        raise ValueError("APC cells contain non-finite values")
    if (cells[["events", "population"]] <= 0).any().any():
        raise ValueError("APC requires positive incidence counts and population")
    cells = cells.drop(columns="_merge")
    cells["age_lower"] = cells.age_name.str.extract(r"^(\d+)-").astype(int)
    grouped = []
    for window, start in WINDOWS.items():
        part = cells[cells.year.between(start, start + 29)].copy()
        part["window"] = window
        part["period_start"] = start + 5 * ((part.year - start) // 5)
        part = part.groupby(
            ["window", "location_name", "sex_name", "age_lower", "period_start"],
            as_index=False,
        )[["events", "population"]].sum()
        expected = len(locations) * len(sexes) * len(AGE_LABELS) * 6
        if len(part) != expected:
            raise ValueError(f"APC {window} has {len(part)} cells; expected {expected}")
        grouped.append(part)
    inputs = pd.concat(grouped, ignore_index=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    input_path = tables_dir / "apc_input_cells.csv"
    inputs.to_csv(input_path, index=False)
    script_dir = Path(__file__).resolve().parent
    completed = subprocess.run(
        [_rscript(), str(script_dir / "nci_apc_runner.R"), str(input_path),
         str(tables_dir), str(script_dir / "nci_apc_reference.R")],
        capture_output=True, text=True,
    )
    if completed.returncode:
        raise RuntimeError(f"APC R analysis failed:\n{completed.stderr}")
    raw = {name: pd.read_csv(tables_dir / f"apc_{name}.csv") for name in RESULT_NAMES}
    primary = {name: frame[frame.window.eq("1994_2023")].drop(columns="window").copy()
               for name, frame in raw.items()}
    sensitivity = {name: frame[frame.window.eq("1990_2019")].drop(columns="window").copy()
                   for name, frame in raw.items()}
    return primary, sensitivity, raw
