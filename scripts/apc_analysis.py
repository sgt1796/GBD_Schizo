"""Descriptive age-period-cohort estimable-function analysis.

This module deliberately operates on GBD point estimates.  It reports no
sampling tests or confidence intervals because posterior draws and their
cross-cell covariance are unavailable.  The linear age-period-cohort identity
is handled by separating two linear trends from orthogonal nonlinear
curvatures; unconstrained age, period, and cohort slopes are never reported as
three independently identifiable effects.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re

import numpy as np
import pandas as pd
from scipy.linalg import null_space


INFERENCE_NOTE = (
    "Descriptive model of GBD posterior-mean point estimates; no GBD-level "
    "hypothesis test or confidence interval is available without posterior draws."
)

BASE_APC_AGES = (
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
)
PREFERRED_APC_AGES = ("10-14 years", *BASE_APC_AGES)


@dataclass(frozen=True)
class APCWindow:
    """A complete sequence of equal-width calendar periods."""

    name: str
    start_year: int
    end_year: int
    period_width: int = 5

    @property
    def period_starts(self) -> tuple[int, ...]:
        return tuple(range(self.start_year, self.end_year + 1, self.period_width))

    @property
    def period_labels(self) -> tuple[str, ...]:
        return tuple(
            f"{start}-{start + self.period_width - 1}"
            for start in self.period_starts
        )

    def period_for_year(self, year: int) -> str | None:
        if not self.start_year <= int(year) <= self.end_year:
            return None
        offset = int(year) - self.start_year
        start = self.start_year + (offset // self.period_width) * self.period_width
        return f"{start}-{start + self.period_width - 1}"


PRIMARY_WINDOW = APCWindow("primary_1994_2023", 1994, 2023)
SENSITIVITY_WINDOW = APCWindow("sensitivity_1990_2019", 1990, 2019)


def age_bounds(label: str) -> tuple[int, int]:
    match = re.fullmatch(r"(\d+)-(\d+) years", str(label).strip())
    if not match:
        raise ValueError(f"APC requires a closed age interval; found {label!r}.")
    lower, upper = (int(value) for value in match.groups())
    if upper < lower:
        raise ValueError(f"Invalid age interval {label!r}.")
    return lower, upper


def age_midpoint(label: str) -> float:
    lower, upper = age_bounds(label)
    return (lower + upper) / 2.0


def validate_equal_age_width(ages: tuple[str, ...], expected_width: int = 5) -> None:
    if len(set(ages)) != len(ages):
        raise ValueError("APC age groups contain duplicates.")
    intervals = [age_bounds(age) for age in ages]
    widths = {upper - lower + 1 for lower, upper in intervals}
    if widths != {expected_width}:
        raise ValueError(
            f"APC age groups must all span {expected_width} years; found widths {sorted(widths)}."
        )
    for previous, current in zip(intervals, intervals[1:]):
        if current[0] != previous[1] + 1:
            raise ValueError("APC age groups must be consecutive and non-overlapping.")


def select_apc_ages(available_ages: set[str]) -> tuple[str, ...]:
    """Prefer 10--69 years, falling back visibly to 15--69 when necessary."""
    if set(PREFERRED_APC_AGES) <= available_ages:
        return PREFERRED_APC_AGES
    if set(BASE_APC_AGES) <= available_ages:
        return BASE_APC_AGES
    missing = sorted(set(BASE_APC_AGES) - available_ages)
    raise ValueError(f"APC input is missing required five-year age groups: {missing}")


def curvature_basis(level_count: int) -> np.ndarray:
    """Return the nonlinear complement to an intercept and a linear trend."""
    if level_count < 3:
        raise ValueError("At least three ordered levels are required for curvature.")
    linear = np.column_stack(
        [np.ones(level_count), np.arange(level_count, dtype=float)]
    )
    return null_space(linear.T)


def _fit_weighted_log_rate(
    rate: np.ndarray,
    design: np.ndarray,
    population: np.ndarray,
    weighting: str = "population",
) -> dict[str, np.ndarray | float | int]:
    rate = np.asarray(rate, dtype=float)
    design = np.asarray(design, dtype=float)
    population = np.asarray(population, dtype=float)
    if (rate <= 0).any() or not np.isfinite(rate).all():
        raise ValueError("APC rates must be finite and positive.")
    if (population <= 0).any() or not np.isfinite(population).all():
        raise ValueError("APC populations must be finite and positive.")
    if weighting == "population":
        root_weight = np.sqrt(population / np.mean(population))
    elif weighting == "equal":
        root_weight = np.ones_like(population)
    else:
        raise ValueError("APC weighting must be 'population' or 'equal'.")
    weighted_design = design * root_weight[:, None]
    weighted_response = np.log(rate) * root_weight
    beta = np.linalg.pinv(weighted_design) @ weighted_response
    fitted = design @ beta
    residual = np.log(rate) - fitted
    return {
        "beta": beta,
        "fitted_log_rate": fitted,
        "residual": residual,
        "rank": int(np.linalg.matrix_rank(weighted_design)),
        "columns": int(design.shape[1]),
        "weighted_rss": float(np.sum((root_weight * residual) ** 2)),
    }


def _validate_source_panel(
    annual: pd.DataFrame,
    ages: tuple[str, ...],
    window: APCWindow,
    locations: tuple[str, ...],
    sexes: tuple[str, ...],
    measures: tuple[str, ...],
) -> None:
    validate_equal_age_width(ages, window.period_width)
    keys = ["location_name", "sex_name", "measure_name", "age_name", "year"]
    duplicates = int(annual.duplicated(keys).sum())
    if duplicates:
        raise ValueError(f"APC input contains {duplicates} duplicated age-year cells.")
    expected_years = window.end_year - window.start_year + 1
    expected_per_panel = len(ages) * expected_years
    expected_total = len(locations) * len(sexes) * len(measures) * expected_per_panel
    if len(annual) != expected_total:
        raise ValueError(
            f"APC input has {len(annual)} annual cells; expected {expected_total} "
            f"({expected_per_panel} per location-sex-outcome panel)."
        )
    sizes = annual.groupby(["location_name", "sex_name", "measure_name"]).size()
    expected_panels = len(locations) * len(sexes) * len(measures)
    if len(sizes) != expected_panels or not sizes.eq(expected_per_panel).all():
        raise ValueError(
            "APC input does not contain complete location-sex-outcome age-year panels."
        )


def build_apc_cells(
    burden: pd.DataFrame,
    population: pd.DataFrame,
    window: APCWindow,
    locations: tuple[str, ...],
    sexes: tuple[str, ...],
    ages: tuple[str, ...] | None = None,
    measures: tuple[str, ...] = ("Incidence",),
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Pool annual outcome estimates into equal five-year APC cells."""
    if len(window.period_starts) != 6:
        raise ValueError("The analysis specification requires exactly six APC periods.")
    if window.period_starts[-1] + window.period_width - 1 != window.end_year:
        raise ValueError("APC window does not end on a complete period boundary.")

    if not measures:
        raise ValueError("At least one APC outcome is required.")
    available_by_measure = {
        measure: set(burden.loc[burden.measure_name.eq(measure), "age_name"])
        for measure in measures
    }
    missing_measures = [
        measure for measure, available in available_by_measure.items() if not available
    ]
    if missing_measures:
        raise ValueError(f"APC input is missing outcomes: {missing_measures}.")
    available = set.intersection(*available_by_measure.values())
    ages = tuple(ages or select_apc_ages(available))
    validate_equal_age_width(ages, window.period_width)
    keys = ["location_name", "sex_name", "age_name", "year"]
    rates = burden[
        burden.measure_name.isin(measures)
        & burden.metric_name.eq("Rate")
        & burden.age_name.isin(ages)
        & burden.location_name.isin(locations)
        & burden.sex_name.isin(sexes)
        & burden.year.between(window.start_year, window.end_year)
    ][["measure_name", *keys, "val"]].rename(columns={"val": "rate"})
    pop = population[
        population.age_name.isin(ages)
        & population.location_name.isin(locations)
        & population.sex_name.isin(sexes)
        & population.year.between(window.start_year, window.end_year)
    ][keys + ["population"]]
    annual = rates.merge(pop, on=keys, how="inner", validate="many_to_one")
    _validate_source_panel(annual, ages, window, locations, sexes, measures)

    annual = annual.copy()
    annual["period"] = annual.year.map(window.period_for_year)
    annual["age_midpoint"] = annual.age_name.map(age_midpoint)
    period_midpoints = {
        label: start + (window.period_width - 1) / 2.0
        for label, start in zip(window.period_labels, window.period_starts)
    }
    annual["period_midpoint"] = annual.period.map(period_midpoints)
    annual["cohort_midpoint"] = annual.period_midpoint - annual.age_midpoint
    annual["events"] = annual.population * annual.rate / 100_000.0
    cells = annual.groupby(
        [
            "location_name",
            "sex_name",
            "measure_name",
            "age_name",
            "age_midpoint",
            "period",
            "period_midpoint",
            "cohort_midpoint",
        ],
        as_index=False,
        observed=True,
    ).agg(events=("events", "sum"), population=("population", "sum"))
    cells["rate"] = cells.events / cells.population * 100_000.0
    cells["window"] = window.name
    cells["age_coverage"] = f"{ages[0].replace(' years', '')} to {ages[-1].replace(' years', '')}"
    return cells, ages


def _validate_cell_matrix(
    cells: pd.DataFrame, ages: tuple[str, ...], window: APCWindow
) -> None:
    keys = ["age_name", "period"]
    if cells.duplicated(keys).any():
        raise ValueError("APC pooled matrix contains duplicated age-period cells.")
    expected = pd.MultiIndex.from_product(
        [ages, window.period_labels], names=keys
    )
    observed = pd.MultiIndex.from_frame(cells[keys])
    missing = expected.difference(observed)
    extra = observed.difference(expected)
    if len(missing) or len(extra):
        raise ValueError(
            f"APC pooled matrix is incomplete: {len(missing)} missing and {len(extra)} extra cells."
        )


def _fit_panel(
    panel: pd.DataFrame,
    ages: tuple[str, ...],
    window: APCWindow,
    weighting: str = "population",
) -> dict[str, pd.DataFrame]:
    _validate_cell_matrix(panel, ages, window)
    panel = panel.copy()
    measures = panel.measure_name.drop_duplicates().tolist()
    if len(measures) != 1:
        raise ValueError("Each APC model panel must contain exactly one outcome.")
    measure = measures[0]
    age_lookup = {value: index for index, value in enumerate(ages)}
    period_lookup = {
        value: index for index, value in enumerate(window.period_labels)
    }
    panel["age_index"] = panel.age_name.map(age_lookup).astype(int)
    panel["period_index"] = panel.period.map(period_lookup).astype(int)
    panel["cohort_index"] = panel.period_index - panel.age_index
    cohorts = tuple(sorted(panel.cohort_index.unique()))
    cohort_lookup = {value: index for index, value in enumerate(cohorts)}
    panel["cohort_level"] = panel.cohort_index.map(cohort_lookup).astype(int)

    age_basis = curvature_basis(len(ages))
    period_basis = curvature_basis(len(window.period_labels))
    cohort_basis = curvature_basis(len(cohorts))
    age_center = panel.age_index - (len(ages) - 1) / 2.0
    period_center = panel.period_index - (len(window.period_labels) - 1) / 2.0
    design = np.column_stack(
        [
            np.ones(len(panel)),
            age_center,
            period_center,
            age_basis[panel.age_index],
            period_basis[panel.period_index],
            cohort_basis[panel.cohort_level],
        ]
    )
    fit = _fit_weighted_log_rate(
        panel.rate.to_numpy(), design, panel.population.to_numpy(), weighting
    )
    if fit["rank"] != fit["columns"]:
        raise ValueError(
            f"APC estimable-function design is rank deficient ({fit['rank']} of {fit['columns']})."
        )
    beta = np.asarray(fit["beta"])
    age_start = 3
    period_start = age_start + age_basis.shape[1]
    cohort_start = period_start + period_basis.shape[1]
    age_nonlinear = age_basis @ beta[age_start:period_start]
    period_nonlinear = period_basis @ beta[period_start:cohort_start]
    cohort_nonlinear = cohort_basis @ beta[cohort_start:]

    reference_age = "40-44 years" if "40-44 years" in ages else ages[len(ages) // 2]
    reference_period = window.period_labels[len(window.period_labels) // 2]
    reference_cohort_level = len(cohorts) // 2
    reference_cohort_index = cohorts[reference_cohort_level]
    cohort_midpoint_lookup = (
        panel[["cohort_index", "cohort_midpoint"]]
        .drop_duplicates()
        .set_index("cohort_index")
        .cohort_midpoint
        .to_dict()
    )

    age_linear = beta[1] * (np.arange(len(ages)) - (len(ages) - 1) / 2.0)
    age_curve = age_linear + age_nonlinear
    reference_age_index = ages.index(reference_age)
    age_rows = pd.DataFrame(
        {
            "age_name": ages,
            "age_midpoint": [age_midpoint(age) for age in ages],
            "longitudinal_age_rr": np.exp(age_curve - age_curve[reference_age_index]),
            "reference_age": reference_age,
        }
    )
    period_rows = pd.DataFrame(
        {
            "period": window.period_labels,
            "period_midpoint": [
                start + (window.period_width - 1) / 2.0
                for start in window.period_starts
            ],
            "period_rr": np.exp(
                period_nonlinear
                - period_nonlinear[window.period_labels.index(reference_period)]
            ),
            "reference_period": reference_period,
        }
    )
    cohort_rows = pd.DataFrame(
        {
            "cohort_index": cohorts,
            "cohort_midpoint": [cohort_midpoint_lookup[value] for value in cohorts],
            "cohort_rr": np.exp(
                cohort_nonlinear - cohort_nonlinear[reference_cohort_level]
            ),
            "reference_cohort_midpoint": cohort_midpoint_lookup[
                reference_cohort_index
            ],
        }
    )

    local_rows = []
    for age in ages:
        age_panel = panel[panel.age_name.eq(age)].sort_values("period_midpoint")
        x = age_panel.period_midpoint.to_numpy(float)
        x = x - x.mean()
        local_design = np.column_stack([np.ones(len(x)), x])
        local_fit = _fit_weighted_log_rate(
            age_panel.rate.to_numpy(),
            local_design,
            age_panel.population.to_numpy(),
            weighting,
        )
        local_rows.append(
            {
                "age_name": age,
                "age_midpoint": age_midpoint(age),
                "local_drift": 100.0 * (math.exp(float(local_fit["beta"][1])) - 1.0),
                "formal_inference_performed": False,
                "inference_note": INFERENCE_NOTE,
            }
        )

    net_drift = 100.0 * (
        math.exp(float(beta[2]) / window.period_width) - 1.0
    )
    summary = pd.DataFrame(
        [
            {
                "measure_name": measure,
                "window": window.name,
                "start_year": window.start_year,
                "end_year": window.end_year,
                "period_count": len(window.period_labels),
                "period_width_years": window.period_width,
                "age_group_count": len(ages),
                "age_coverage": f"{ages[0].replace(' years', '')} to {ages[-1].replace(' years', '')}",
                "net_drift": net_drift,
                "reference_age": reference_age,
                "reference_period": reference_period,
                "reference_cohort_midpoint": cohort_midpoint_lookup[
                    reference_cohort_index
                ],
                "design_rank": fit["rank"],
                "design_columns": fit["columns"],
                "weighted_log_rate_rss": fit["weighted_rss"],
                "weighting": weighting,
                "formal_inference_performed": False,
                "identifiability_note": (
                    "The exact age = period - cohort dependency is handled by two "
                    "linear trends plus orthogonal nonlinear curvatures; three "
                    "unconstrained linear effects are not estimated."
                ),
                "interpretation": INFERENCE_NOTE,
            }
        ]
    )
    panel["fitted_rate"] = np.exp(np.asarray(fit["fitted_log_rate"]))
    return {
        "summary": summary,
        "local_drift": pd.DataFrame(local_rows),
        "age_curve": age_rows,
        "period_rr": period_rows,
        "cohort_rr": cohort_rows,
        "cells": panel,
    }


def run_apc(
    burden: pd.DataFrame,
    population: pd.DataFrame,
    window: APCWindow,
    locations: tuple[str, ...],
    sexes: tuple[str, ...],
    ages: tuple[str, ...] | None = None,
    measures: tuple[str, ...] = ("Incidence",),
    weighting: str = "population",
) -> dict[str, pd.DataFrame]:
    """Run APC point-estimate analyses for every requested outcome panel."""
    cells, selected_ages = build_apc_cells(
        burden, population, window, locations, sexes, ages, measures
    )
    outputs = {
        "summary": [],
        "local_drift": [],
        "age_curve": [],
        "period_rr": [],
        "cohort_rr": [],
        "cells": [],
    }
    for (location, sex, measure), panel in cells.groupby(
        ["location_name", "sex_name", "measure_name"], sort=True
    ):
        fitted = _fit_panel(panel, selected_ages, window, weighting)
        for name, frame in fitted.items():
            frame = frame.copy()
            if "sex_name" not in frame:
                frame.insert(0, "sex_name", sex)
            if "location_name" not in frame:
                frame.insert(0, "location_name", location)
            if "measure_name" not in frame:
                frame.insert(2, "measure_name", measure)
            if "weighting" not in frame:
                frame.insert(3, "weighting", weighting)
            outputs[name].append(frame)
    return {
        name: pd.concat(frames, ignore_index=True)
        for name, frames in outputs.items()
    }
