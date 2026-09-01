from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import publication_analysis as pa


@pytest.fixture(scope="module")
def burden():
    return pa.load_burden(pa.DEFAULT_BURDEN)


@pytest.fixture(scope="module")
def proxy_population(burden):
    return pa.infer_proxy_population(burden)


def test_primary_outcomes_and_years_complete(burden):
    primary = burden[
        burden.measure_name.isin(pa.OUTCOMES)
        & (burden.metric_name == "Rate")
        & (burden.age_name == pa.ASR)
    ]
    sizes = primary.groupby(["location_name", "sex_name", "measure_name"]).year.nunique()
    assert len(sizes) == 12
    assert (sizes == 34).all()


def test_yld_daly_are_numerically_identical(burden):
    result = pa.verify_yld_daly_identity(burden).iloc[0]
    assert bool(result.numerically_identical)
    assert result.max_relative_difference_val < 1e-7


def test_proxy_population_is_complete_and_reconstructs(burden, proxy_population):
    pa.validate_population(proxy_population)
    _, _, reconstruction = pa.audit_burden(burden, proxy_population)
    assert reconstruction.relative_error_pct.abs().quantile(0.99) < 1e-6


def test_official_population_requires_release_marker(tmp_path, proxy_population):
    path = tmp_path / "population.csv"
    proxy_population.drop(columns="population_source").to_csv(path, index=False)
    with pytest.raises(ValueError, match="GBD 2023"):
        pa.load_official_population(path, "unknown")


def test_apc_periods_are_equal_width_and_prespecified():
    assert pa.apc_period(1993) is None
    assert pa.apc_period(1994) == "1994-1998"
    assert pa.apc_period(1998) == "1994-1998"
    assert pa.apc_period(1999) == "1999-2003"
    assert pa.apc_period(2023) == "2019-2023"


def test_decomposition_closes_exactly():
    p0 = np.array([100.0, 200.0, 300.0])
    p1 = np.array([120.0, 250.0, 360.0])
    r0 = np.array([0.01, 0.02, 0.03])
    r1 = np.array([0.015, 0.018, 0.035])
    out = pa.decompose_change(p0, p1, r0, r1)
    assert abs(out["closure_error"]) < 1e-10
    expected = np.sum(p1 * r1) - np.sum(p0 * r0)
    assert out["total_change"] == pytest.approx(expected)


def test_bh_adjustment_is_monotone_in_rank():
    p = pd.Series([0.01, 0.04, 0.03, 0.20])
    q = pa.bh_adjust(p)
    ordered = q[np.argsort(p.to_numpy())]
    assert np.all(np.diff(ordered) >= -1e-12)
    assert np.all((0 <= q) & (q <= 1))


def test_segmented_model_detects_strong_change():
    years = np.arange(1990, 2024)
    x = years - 1990
    log_rate = 2.0 + 0.01 * x + 0.09 * np.maximum(0, years - 2005)
    panel = pd.DataFrame({
        "location_name": "China", "sex_name": "Female", "measure_name": "Incidence",
        "year": years, "val": np.exp(log_rate), "lower": np.exp(log_rate) * 0.9,
        "upper": np.exp(log_rate) * 1.1,
    })
    summary, _, _ = pa.segmented_summary(panel, 199, np.random.default_rng(12))
    assert summary["joinpoint_count"] >= 1
    knots = [int(x) for x in summary["joinpoint_years"].split(",") if x]
    assert min(abs(k - 2005) for k in knots) <= 1


def test_endpoint_outputs_exclude_percent_and_duplicate_yld(burden):
    table = pa.endpoint_table(burden)
    assert set(table.measure_name) == set(pa.OUTCOMES)
    assert not table.metric_name.str.contains("Percent", case=False).any()


def test_secondary_apc_dimensions(burden, proxy_population):
    apc = pa.run_secondary_apc(burden, proxy_population)
    assert set(apc["summary"].measure_name) == {"Incidence"}
    assert apc["cells"].period.nunique() == 6
    assert apc["cells"].age_name.nunique() == 11

