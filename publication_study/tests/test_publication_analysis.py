import json
from pathlib import Path
import sys
from types import SimpleNamespace

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


@pytest.fixture(scope="module")
def descriptive_trends(burden):
    return pa.run_segmented(burden)


@pytest.fixture(scope="module")
def apc_results(burden, proxy_population):
    return pa.run_secondary_apc(burden, proxy_population)


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
    assert bool(result.complete_expected_panel)
    assert result.matched_cells == result.expected_cells == 3944
    assert result.duplicate_cells == 0
    for field in ("val", "lower", "upper"):
        assert result[f"max_relative_difference_{field}"] < pa.YLD_DALY_IDENTITY_TOLERANCE


def test_burden_cause_is_verified(tmp_path):
    bad = pd.read_csv(pa.DEFAULT_BURDEN, nrows=1)
    bad["cause_name"] = "Not schizophrenia"
    path = tmp_path / "wrong_cause.csv"
    bad.to_csv(path, index=False)
    with pytest.raises(ValueError, match="only cause_name='Schizophrenia'"):
        pa.load_burden(path)


def test_proxy_population_is_complete_and_reconstructs(burden, proxy_population):
    pa.validate_population(proxy_population)
    _, _, reconstruction = pa.audit_burden(burden, proxy_population)
    assert reconstruction.relative_error_pct.abs().quantile(0.99) < 1e-6
    assert len(proxy_population) == 1768
    all_age = pa.all_age_count_reconstruction(burden, reconstruction)
    assert len(all_age) == 408
    assert all_age.within_tolerance.all()
    assert all_age.relative_error_pct.abs().max() < 1e-10


def test_all_age_decomposition_schema_is_explicit(burden, proxy_population):
    audit, _, _ = pa.audit_burden(burden, proxy_population)
    decomposition = pa.run_decomposition(
        burden, proxy_population, windows=((1990, 2023),)
    )
    assert {"all_age_groups", "all_age_year_cells"} <= set(audit.columns)
    assert (audit.all_age_groups == 13).all()
    assert (audit.all_age_year_cells == 13 * 34).all()
    assert not any(column.startswith(("adult_", "age_15_plus_")) for column in audit.columns)
    assert {
        "all_age_count_start_reconstructed",
        "all_age_count_end_reconstructed",
    } <= set(decomposition.columns)
    assert not any(column.startswith(("adult_", "age_15_plus_")) for column in decomposition.columns)
    assert pa.DECOMPOSITION_AGES[0] == "0-14 years"
    assert "0-14 years" not in pa.APC_AGES
    assert "70+ years" not in pa.APC_AGES


def test_official_population_requires_release_marker(tmp_path, proxy_population):
    path = tmp_path / "population.csv"
    population = proxy_population.drop(columns="population_source").copy()
    population["gbd_release"] = "GBD 2023"
    population.to_csv(path, index=False)
    with pytest.raises(ValueError, match="GBD 2023"):
        pa.load_official_population(path, "unknown")

    population["gbd_release"] = "GBD 2021"
    population.to_csv(path, index=False)
    with pytest.raises(ValueError, match="consistently identify GBD 2023"):
        pa.load_official_population(path, "GBD 2023")

    population.drop(columns="gbd_release").to_csv(path, index=False)
    with pytest.raises(ValueError, match="gbd_release"):
        pa.load_official_population(path, "GBD 2023")


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


def test_trend_direction_is_explicit_and_handles_missing_values():
    assert pa.trend_direction(0.2) == "increase"
    assert pa.trend_direction(-0.2) == "decrease"
    assert pa.trend_direction(0.0) == "stable"
    assert pa.trend_direction(np.nan) == "not available"


def test_descriptive_bic_curve_detects_strong_change_without_inference():
    years = np.arange(1990, 2024)
    x = years - 1990
    noise = np.random.default_rng(12).normal(0.0, 0.002, len(years))
    log_rate = 2.0 + 0.01 * x + 0.09 * np.maximum(0, years - 2005) + noise
    panel = pd.DataFrame({
        "location_name": "China", "sex_name": "Female", "measure_name": "Incidence",
        "year": years, "val": np.exp(log_rate), "lower": np.exp(log_rate) * 0.9,
        "upper": np.exp(log_rate) * 1.1,
    })
    summary, _, _ = pa.segmented_summary(panel)
    assert summary["joinpoint_count"] >= 1
    knots = [int(x) for x in summary["joinpoint_years"].split(",") if x]
    assert min(abs(k - 2005) for k in knots) <= 1
    assert summary["selection_method"].startswith("minimum BIC")
    assert not summary["formal_inference_performed"]
    assert "permutation_p_0_vs_1" not in summary
    assert "aapc_lower_model_ci" not in summary


def test_residual_autocorrelation_is_reported_not_used_for_inference():
    diagnostics = pa.residual_diagnostics(np.arange(1.0, 10.0))
    assert diagnostics["lag1_residual_autocorrelation"] == pytest.approx(1.0)
    assert diagnostics["material_residual_autocorrelation"]
    assert diagnostics["durbin_watson"] > 0


def test_trajectory_contrasts_are_descriptive(burden):
    contrasts = pa.build_trajectory_contrasts(burden)
    assert len(contrasts) == 12
    assert "p_value" not in contrasts
    assert "q_value_bh_within_family" not in contrasts
    assert not contrasts.formal_inference_performed.any()
    assert np.isfinite(
        contrasts.annualized_endpoint_change_difference_b_minus_a_pct_points
    ).all()
    assert np.isfinite(contrasts.rms_annual_log_change_difference).all()


def test_endpoint_outputs_exclude_percent_and_duplicate_yld(burden):
    table = pa.endpoint_table(burden)
    assert set(table.measure_name) == set(pa.OUTCOMES)
    assert not table.metric_name.str.contains("Percent", case=False).any()


def test_secondary_apc_dimensions(apc_results):
    apc = apc_results
    assert set(apc["summary"].measure_name) == {"Incidence"}
    assert apc["cells"].period.nunique() == 6
    assert apc["cells"].age_name.nunique() == 11
    assert "net_drift_lower_model_ci" not in apc["summary"]
    assert not apc["summary"].formal_inference_performed.any()


def test_apc_primary_direction_agreement_is_reported(
    burden, descriptive_trends, apc_results
):
    primary, _, _ = descriptive_trends
    comparison = pa.compare_primary_apc_directions(
        burden, primary, apc_results["summary"]
    )
    assert len(comparison) == 4
    expected = comparison.apc_net_drift_direction == comparison.primary_segmented_direction
    assert comparison.apc_vs_segmented_direction_agreement.equals(expected)
    assert set(comparison.primary_segmented_direction) <= {"increase", "decrease", "stable"}
    assert set(comparison.apc_net_drift_direction) <= {"increase", "decrease", "stable"}


def test_pandemic_sensitivity_windows_are_unambiguous(
    burden, proxy_population
):
    trend, _, _ = pa.run_segmented(burden, end_year=2019)
    assert (trend.start_year == 1990).all()
    assert (trend.end_year == 2019).all()
    apc = pa.run_secondary_apc(burden, proxy_population, include_last_period=False)
    assert (apc["summary"].end_year == 2018).all()
    assert "2019-2023" not in set(apc["cells"].period)


def test_nci_version_is_recorded_without_hard_coding_one_release():
    versions = pd.Series(["NCI Joinpoint 6.1.0", "6.1.0"])
    assert pa._consistent_software_version(versions) == "6.1.0"
    with pytest.raises(ValueError, match="mix software versions"):
        pa._consistent_software_version(pd.Series(["6.0.1", "6.1.0"]))


def test_portable_workspace_paths_use_forward_slashes():
    assert pa.portable_path(pa.DEFAULT_BURDEN) == "prepared_inputs/cause_all.csv"


def test_sha256_fingerprint_is_content_based(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("abc", encoding="utf-8")
    assert pa.file_sha256(source) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_table_writer_removes_only_known_stale_outputs(tmp_path):
    for filename in pa.OBSOLETE_TABLE_CSVS:
        (tmp_path / filename).write_text("obsolete", encoding="utf-8")
    for name in pa.OPTIONAL_NCI_TABLES:
        (tmp_path / f"{name}.csv").write_text("stale optional", encoding="utf-8")
    unrelated = tmp_path / "user_notes.csv"
    unrelated.write_text("preserve me", encoding="utf-8")

    pa.write_tables({"current_table": pd.DataFrame({"value": [1]})}, tmp_path)

    assert all(not (tmp_path / filename).exists() for filename in pa.OBSOLETE_TABLE_CSVS)
    assert all(not (tmp_path / f"{name}.csv").exists() for name in pa.OPTIONAL_NCI_TABLES)
    assert unrelated.read_text(encoding="utf-8") == "preserve me"
    assert (tmp_path / "current_table.csv").exists()
    assert (tmp_path / "publication_tables.xlsx").exists()


def test_full_descriptive_pipeline_smoke(tmp_path):
    args = SimpleNamespace(
        output_dir=tmp_path,
        burden_csv=pa.DEFAULT_BURDEN,
        population_csv=None,
        population_release="GBD 2023",
        allow_proxy_population=True,
        nci_results_csv=None,
    )
    result = pa.run(args)
    assert not result["metadata"]["submission_ready"]
    assert not result["metadata"]["formal_trend_inference_performed"]
    assert "trajectory_contrasts" in result["tables"]
    assert "all_age_count_reconstruction" in result["tables"]
    assert "trend_excluding_2020_2023" in result["tables"]
    assert "apc_excluding_2019_2023" in result["tables"]
    assert "apc_primary_direction_agreement" in result["tables"]
    assert len(result["metadata"]["burden_csv_sha256"]) == 64
    validation = json.loads((tmp_path / "qa" / "validation_summary.json").read_text(encoding="utf-8"))
    assert validation["all_age_count_reconstruction_rows"] == 408
    assert validation["all_age_count_reconstruction_within_tolerance"]
    assert (tmp_path / "qa" / "validation_summary.json").exists()

