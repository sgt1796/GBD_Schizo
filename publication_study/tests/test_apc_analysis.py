from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import apc_analysis as aa


def synthetic_apc_inputs():
    ages = aa.BASE_APC_AGES
    window = aa.PRIMARY_WINDOW
    age_basis = aa.curvature_basis(len(ages))
    period_basis = aa.curvature_basis(len(window.period_labels))
    cohort_values = tuple(range(-(len(ages) - 1), len(window.period_labels)))
    cohort_basis = aa.curvature_basis(len(cohort_values))
    cohort_lookup = {value: index for index, value in enumerate(cohort_values)}
    age_coefficients = np.linspace(-0.03, 0.03, age_basis.shape[1])
    period_coefficients = np.linspace(0.02, -0.02, period_basis.shape[1])
    cohort_coefficients = np.linspace(-0.015, 0.015, cohort_basis.shape[1])
    age_linear = 0.08
    period_linear = 0.025

    burden_rows = []
    population_rows = []
    for age_index, age in enumerate(ages):
        for year in range(window.start_year, window.end_year + 1):
            period_index = (year - window.start_year) // window.period_width
            cohort_index = period_index - age_index
            eta = (
                4.0
                + age_linear * (age_index - (len(ages) - 1) / 2.0)
                + period_linear
                * (period_index - (len(window.period_labels) - 1) / 2.0)
                + age_basis[age_index] @ age_coefficients
                + period_basis[period_index] @ period_coefficients
                + cohort_basis[cohort_lookup[cohort_index]] @ cohort_coefficients
            )
            burden_rows.append(
                {
                    "location_name": "Test location",
                    "sex_name": "Test sex",
                    "age_name": age,
                    "year": year,
                    "measure_name": "Incidence",
                    "metric_name": "Rate",
                    "val": np.exp(eta),
                }
            )
            population_rows.append(
                {
                    "location_name": "Test location",
                    "sex_name": "Test sex",
                    "age_name": age,
                    "year": year,
                    "population": 100_000.0 + 500.0 * age_index,
                }
            )
    expected = {
        "net_drift": 100.0 * (np.exp(period_linear / window.period_width) - 1.0),
        "age_curve": age_linear
        * (np.arange(len(ages)) - (len(ages) - 1) / 2.0)
        + age_basis @ age_coefficients,
        "period_curve": period_basis @ period_coefficients,
        "cohort_curve": cohort_basis @ cohort_coefficients,
    }
    return pd.DataFrame(burden_rows), pd.DataFrame(population_rows), expected


def test_apc_windows_are_six_complete_five_year_periods():
    assert aa.PRIMARY_WINDOW.period_labels == (
        "1994-1998",
        "1999-2003",
        "2004-2008",
        "2009-2013",
        "2014-2018",
        "2019-2023",
    )
    assert aa.SENSITIVITY_WINDOW.period_labels == (
        "1990-1994",
        "1995-1999",
        "2000-2004",
        "2005-2009",
        "2010-2014",
        "2015-2019",
    )


def test_apc_age_intervals_must_be_equal_and_consecutive():
    aa.validate_equal_age_width(aa.BASE_APC_AGES)
    with pytest.raises(ValueError, match="must all span 5 years"):
        aa.validate_equal_age_width(("0-14 years", *aa.BASE_APC_AGES))
    with pytest.raises(ValueError, match="consecutive"):
        aa.validate_equal_age_width(("15-19 years", "25-29 years", "30-34 years"))


def test_cohort_index_and_complete_matrix_construction():
    burden, population, _ = synthetic_apc_inputs()
    cells, ages = aa.build_apc_cells(
        burden,
        population,
        aa.PRIMARY_WINDOW,
        ("Test location",),
        ("Test sex",),
        aa.BASE_APC_AGES,
    )
    assert ages == aa.BASE_APC_AGES
    assert len(cells) == len(ages) * 6
    assert cells[["age_name", "period"]].drop_duplicates().shape[0] == len(cells)
    assert np.allclose(
        cells.cohort_midpoint,
        cells.period_midpoint - cells.age_midpoint,
    )


def test_apc_synthetic_estimands_are_recovered():
    burden, population, expected = synthetic_apc_inputs()
    result = aa.run_apc(
        burden,
        population,
        aa.PRIMARY_WINDOW,
        ("Test location",),
        ("Test sex",),
        aa.BASE_APC_AGES,
    )
    summary = result["summary"].iloc[0]
    assert summary.design_rank == summary.design_columns
    assert summary.net_drift == pytest.approx(expected["net_drift"], abs=1e-9)

    age = result["age_curve"]
    reference_age = list(aa.BASE_APC_AGES).index("40-44 years")
    expected_age_rr = np.exp(
        expected["age_curve"] - expected["age_curve"][reference_age]
    )
    assert age.longitudinal_age_rr.to_numpy() == pytest.approx(
        expected_age_rr, abs=1e-9
    )

    period = result["period_rr"]
    reference_period = len(aa.PRIMARY_WINDOW.period_labels) // 2
    expected_period_rr = np.exp(
        expected["period_curve"] - expected["period_curve"][reference_period]
    )
    assert period.period_rr.to_numpy() == pytest.approx(
        expected_period_rr, abs=1e-9
    )

    cohort = result["cohort_rr"]
    reference_cohort = len(expected["cohort_curve"]) // 2
    expected_cohort_rr = np.exp(
        expected["cohort_curve"] - expected["cohort_curve"][reference_cohort]
    )
    assert cohort.cohort_rr.to_numpy() == pytest.approx(
        expected_cohort_rr, abs=1e-9
    )
    assert age.loc[age.age_name.eq("40-44 years"), "longitudinal_age_rr"].iloc[0] == pytest.approx(1.0)
    assert period.loc[period.period.eq(summary.reference_period), "period_rr"].iloc[0] == pytest.approx(1.0)
    assert cohort.loc[
        cohort.cohort_midpoint.eq(summary.reference_cohort_midpoint), "cohort_rr"
    ].iloc[0] == pytest.approx(1.0)


def test_apc_runs_every_requested_outcome_without_cross_panel_mixing():
    burden, population, expected = synthetic_apc_inputs()
    measures = ("Incidence", "Prevalence", "DALYs")
    burden = pd.concat(
        [burden.assign(measure_name=measure) for measure in measures],
        ignore_index=True,
    )
    result = aa.run_apc(
        burden,
        population,
        aa.PRIMARY_WINDOW,
        ("Test location",),
        ("Test sex",),
        aa.BASE_APC_AGES,
        measures,
    )
    assert set(result["summary"].measure_name) == set(measures)
    assert len(result["summary"]) == len(measures)
    assert result["summary"].net_drift.to_numpy() == pytest.approx(
        np.repeat(expected["net_drift"], len(measures)), abs=1e-9
    )
    for name, frame in result.items():
        assert set(frame.measure_name) == set(measures), name


def test_apc_missing_and_duplicated_cells_fail_loudly():
    burden, population, _ = synthetic_apc_inputs()
    with pytest.raises(ValueError):
        aa.build_apc_cells(
            burden.iloc[:-1],
            population,
            aa.PRIMARY_WINDOW,
            ("Test location",),
            ("Test sex",),
            aa.BASE_APC_AGES,
        )
    duplicated = pd.concat([burden, burden.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError):
        aa.build_apc_cells(
            duplicated,
            population,
            aa.PRIMARY_WINDOW,
            ("Test location",),
            ("Test sex",),
            aa.BASE_APC_AGES,
        )
