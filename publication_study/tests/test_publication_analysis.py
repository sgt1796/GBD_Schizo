import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import publication_analysis as pa
import build_documents as bd


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
    pa.validate_required_burden_summaries(burden)
    primary = burden[
        burden.measure_name.isin(pa.OUTCOMES)
        & (burden.metric_name == "Rate")
        & (burden.age_name == pa.ASR)
    ]
    sizes = primary.groupby(["location_name", "sex_name", "measure_name"]).year.nunique()
    assert len(sizes) == 12
    assert (sizes == 34).all()


def test_missing_primary_summary_panel_fails_early(burden):
    missing = burden.drop(
        burden[
            (burden.location_name == "China")
            & (burden.sex_name == "Female")
            & (burden.measure_name == "Incidence")
            & (burden.age_name == pa.ALL_AGES)
            & (burden.metric_name == "Number")
            & (burden.year == 1990)
        ].index
    )
    with pytest.raises(ValueError, match="all-age Number"):
        pa.validate_required_burden_summaries(missing)


def test_yld_daly_are_numerically_identical(burden):
    result = pa.verify_yld_daly_identity(burden).iloc[0]
    assert bool(result.audit_passed)
    assert result.audit_status == "verified_identical"
    assert bool(result.numerically_identical)
    assert bool(result.complete_expected_panel)
    assert result.matched_cells == result.expected_cells == 3944
    assert result.duplicate_cells == 0
    for field in ("val", "lower", "upper"):
        assert result[f"max_relative_difference_{field}"] < pa.YLD_DALY_IDENTITY_TOLERANCE


def test_yld_panel_is_optional_but_partial_yld_input_fails(burden):
    without_yld = burden[burden.measure_name != "YLDs"].copy()
    absent = pa.verify_yld_daly_identity(without_yld).iloc[0]
    assert bool(absent.audit_passed)
    assert not bool(absent.yld_panel_available)
    assert absent.audit_status == "not_available"
    assert not bool(absent.numerically_identical)

    yld_index = burden.index[burden.measure_name == "YLDs"][0]
    partial = pa.verify_yld_daly_identity(burden.drop(index=yld_index)).iloc[0]
    assert not bool(partial.audit_passed)
    assert bool(partial.yld_panel_available)
    assert partial.audit_status == "incomplete_or_nonidentical"


def test_document_context_tracks_production_inputs_without_proxy_wording():
    tables = {
        "decomposition": pd.DataFrame({"age_group_count": [20]}),
        "apc_summary": pd.DataFrame({"age_coverage": ["10-14 to 65-69"]}),
        "yld_daly_identity": pd.DataFrame({"yld_panel_available": [False]}),
    }
    context = bd.document_build_context({"submission_ready": True}, tables)
    assert context == {
        "submission_ready": True,
        "age_count": 20,
        "age_span": "0-4 through 95+ years",
        "apc_coverage": "10-14 to 65-69",
        "yld_available": False,
    }


def _write_synthetic_fine_production_inputs(tmp_path, burden, proxy_population):
    population_rows = []
    burden_rows = []
    age_weights = np.linspace(1.0, 2.0, len(pa.FINE_DECOMPOSITION_AGES))
    profiles = {
        "Incidence": np.exp(-((np.arange(20) - 4.0) / 3.2) ** 2) + 0.08,
        "Prevalence": np.exp(-((np.arange(20) - 8.0) / 6.0) ** 2) + 0.18,
        "DALYs": np.exp(-((np.arange(20) - 7.0) / 5.5) ** 2) + 0.16,
    }
    total_population = proxy_population.groupby(
        ["location_name", "sex_name", "year"], as_index=False
    ).population.sum()

    for population_cell in total_population.itertuples(index=False):
        year_shift = 1.0 + 0.001 * (population_cell.year - 1990) * np.linspace(-1, 1, 20)
        shares = age_weights * year_shift
        shares = shares / shares.sum()
        populations = float(population_cell.population) * shares
        for age_name, population in zip(pa.FINE_DECOMPOSITION_AGES, populations):
            population_rows.append({
                "location_name": population_cell.location_name,
                "sex_name": population_cell.sex_name,
                "age_name": age_name,
                "year": population_cell.year,
                "population": population,
                "gbd_release": "GBD 2023",
            })

        for outcome in pa.OUTCOMES:
            all_age = burden[
                burden.location_name.eq(population_cell.location_name)
                & burden.sex_name.eq(population_cell.sex_name)
                & burden.measure_name.eq(outcome)
                & burden.age_name.eq(pa.ALL_AGES)
                & burden.metric_name.eq("Number")
                & burden.year.eq(population_cell.year)
            ].iloc[0]
            raw_counts = populations * profiles[outcome]
            counts = float(all_age.val) * raw_counts / raw_counts.sum()
            rates = 100000.0 * counts / populations
            for age_name, count, rate in zip(pa.FINE_DECOMPOSITION_AGES, counts, rates):
                for metric_name, value in (("Number", count), ("Rate", rate)):
                    burden_rows.append({
                        "location_name": population_cell.location_name,
                        "sex_name": population_cell.sex_name,
                        "age_name": age_name,
                        "measure_name": outcome,
                        "metric_name": metric_name,
                        "cause_name": "Schizophrenia",
                        "year": population_cell.year,
                        "val": value,
                        "lower": value * 0.9,
                        "upper": value * 1.1,
                    })

    summary = burden[
        burden.measure_name.isin(pa.OUTCOMES)
        & (
            (burden.age_name.eq(pa.ALL_AGES) & burden.metric_name.eq("Number"))
            | (burden.age_name.eq(pa.ASR) & burden.metric_name.eq("Rate"))
        )
    ][[
        "location_name", "sex_name", "age_name", "measure_name", "metric_name",
        "cause_name", "year", "val", "lower", "upper",
    ]]
    burden_frame = pd.concat([summary, pd.DataFrame(burden_rows)], ignore_index=True)
    population_frame = pd.DataFrame(population_rows)
    burden_path = tmp_path / "synthetic_fine_burden.csv"
    population_path = tmp_path / "synthetic_population.csv"
    burden_frame.to_csv(burden_path, index=False)
    population_frame.to_csv(population_path, index=False)

    common_dimensions = {
        "locations": ["China", "United States of America"],
        "sexes": ["Female", "Male"],
        "years": "1990-2023",
    }
    metadata_paths = []
    for role, data_path, extra_dimensions in (
        (
            "burden",
            burden_path,
            {
                "ages": [*pa.FINE_DECOMPOSITION_AGES, pa.ALL_AGES, pa.ASR],
                "measures": list(pa.OUTCOMES),
                "metrics": ["Number", "Rate"],
            },
        ),
        (
            "population",
            population_path,
            {
                "ages": list(pa.FINE_DECOMPOSITION_AGES),
                "measures": ["Population"],
                "metrics": ["Number"],
            },
        ),
    ):
        metadata = {
            "export_role": role,
            "gbd_release": "GBD 2023",
            "retrieval_date": "2026-09-02",
            "export_id": f"SYNTHETIC-TEST-{role.upper()}",
            "source_url": "https://vizhub.healthdata.org/gbd-results/",
            "query_dimensions": {**common_dimensions, **extra_dimensions},
            "raw_files": [{"file": data_path.name, "sha256": pa.file_sha256(data_path)}],
        }
        metadata_path = tmp_path / f"{role}_metadata.json"
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        metadata_paths.append(metadata_path)
    return burden_path, population_path, *metadata_paths


def test_full_synthetic_production_pipeline_without_optional_yld(
    tmp_path, monkeypatch, burden, proxy_population
):
    monkeypatch.setattr(pa, "ROOT", tmp_path)
    burden_path, population_path, burden_metadata, population_metadata = (
        _write_synthetic_fine_production_inputs(tmp_path, burden, proxy_population)
    )
    output = tmp_path / "production_output"
    args = SimpleNamespace(
        output_dir=output,
        burden_csv=burden_path,
        population_csv=population_path,
        burden_metadata_json=burden_metadata,
        population_metadata_json=population_metadata,
        population_release="GBD 2023",
        allow_proxy_population=False,
        nci_results_csv=None,
    )
    result = pa.run(args)
    assert result["metadata"]["submission_ready"]
    assert result["metadata"]["fine_age_burden_validated"]
    assert result["metadata"]["population_status"] == "official_GBD_2023"
    assert result["tables"]["yld_daly_identity"].iloc[0].audit_status == "not_available"
    assert result["tables"]["apc_summary"].age_coverage.eq("10-14 to 65-69").all()
    assert result["tables"]["decomposition"].age_group_count.eq(20).all()
    assert result["tables"]["population_source_comparison"].relative_difference_pct.abs().max() < 1e-10
    assert "YLD" not in result["tables"]["provenance"].iloc[0].dimensions

    documents = output / "documents"
    documents.mkdir()
    manuscript = bd.build_manuscript(output, documents, result["metadata"], result["tables"])
    text = "\n".join(paragraph.text for paragraph in bd.Document(manuscript).paragraphs)
    assert "PROVISIONAL ANALYTICAL DRAFT" not in text
    assert "official GBD 2023 population estimates" in text
    assert "did not include YLDs" in text


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


def test_proxy_population_backfills_one_zero_burden_age_from_all_age_rate(burden):
    modified = burden.copy()
    zero_age = pa.PROVISIONAL_DECOMPOSITION_AGES[0]
    mask = (
        modified.age_name.eq(zero_age)
        & modified.measure_name.isin((*pa.OUTCOMES, "YLDs"))
        & modified.metric_name.isin(("Number", "Rate"))
    )
    modified.loc[mask, ["val", "lower", "upper"]] = 0.0

    population = pa.infer_proxy_population(modified)
    pa.validate_population(population)
    backfilled = population[population.age_name.eq(zero_age)]
    assert len(backfilled) == len(pa.LOCATIONS) * len(pa.SEXES) * len(pa.YEARS)
    assert backfilled.population.gt(0).all()


def test_proxy_population_refuses_to_invent_split_for_two_zero_burden_ages(
    burden, proxy_population
):
    modified = burden.copy()
    zero_ages = pa.PROVISIONAL_DECOMPOSITION_AGES[:2]
    mask = (
        modified.age_name.isin(zero_ages)
        & modified.measure_name.isin((*pa.OUTCOMES, "YLDs"))
        & modified.metric_name.isin(("Number", "Rate"))
    )
    modified.loc[mask, ["val", "lower", "upper"]] = 0.0

    with pytest.raises(ValueError, match="missing_ages"):
        pa.infer_proxy_population(modified)

    partial = pa.infer_proxy_population(modified, allow_undefined=True)
    unavailable = partial.population.isna()
    assert unavailable.sum() == (
        len(zero_ages) * len(pa.LOCATIONS) * len(pa.SEXES) * len(pa.YEARS)
    )
    comparison = pa.compare_population_sources(proxy_population, partial)
    assert comparison.reconstruction_available.eq(~unavailable).all()
    assert set(comparison.loc[~comparison.reconstruction_available, "comparison_status"]) == {
        "unavailable_zero_burden_number_and_rate"
    }


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


def test_official_population_requires_complete_fine_age_keys(tmp_path):
    rows = []
    for location in pa.LOCATIONS:
        for sex in pa.SEXES:
            for age in pa.FINE_DECOMPOSITION_AGES:
                for year in pa.YEARS:
                    rows.append({
                        "location_name": location,
                        "sex_name": sex,
                        "age_name": age,
                        "year": year,
                        "population": 100_000 + year,
                        "gbd_release": "GBD 2023",
                    })
    path = tmp_path / "fine_population.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    loaded = pa.load_official_population(path, "GBD 2023")
    assert len(loaded) == 2720
    assert set(loaded.age_name) == set(pa.FINE_DECOMPOSITION_AGES)
    assert set(loaded.population_source) == {"official_GBD_2023"}

    pd.DataFrame(rows[:-1]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="expected 2720"):
        pa.load_official_population(path, "GBD 2023")


def test_apc_primary_periods_are_equal_width_and_recorded():
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
    assert {
        "population_size_change",
        "age_structure_change",
        "age_specific_rate_change",
    } <= set(out)


def test_decomposition_age_bin_sensitivity_is_explicit(burden, proxy_population):
    sensitivity = pa.decomposition_age_bin_sensitivity(burden, proxy_population)
    assert len(sensitivity) == 12
    assert sensitivity.finest_age_group_count.eq(13).all()
    assert sensitivity.collapsed_age_group_count.eq(4).all()
    assert sensitivity.maximum_component_shift_pct_of_total_change.ge(0).all()
    assert sensitivity.material_age_bin_sensitivity.dtype == bool


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


def test_segmented_model_specification_sensitivity_is_descriptive(burden):
    sensitivity = pa.segmented_specification_sensitivity(burden)
    assert len(sensitivity) == 12 * 6
    assert sensitivity.groupby(
        ["location_name", "sex_name", "measure_name"]
    ).specification.nunique().eq(6).all()
    assert set(sensitivity.model_scale) == {"log_rate", "rate"}
    assert not sensitivity.formal_inference_performed.any()
    assert sensitivity.direction_stable_vs_primary.dtype == bool


def test_endpoint_outputs_exclude_percent_and_duplicate_yld(burden):
    table = pa.endpoint_table(burden)
    assert set(table.measure_name) == set(pa.OUTCOMES)
    assert not table.metric_name.str.contains("Percent", case=False).any()


def test_secondary_apc_dimensions(apc_results):
    apc = apc_results
    assert set(apc["summary"].measure_name) == set(pa.OUTCOMES)
    assert len(apc["summary"]) == 12
    assert apc["cells"].period.nunique() == 6
    assert apc["cells"].age_name.nunique() == 11
    assert (apc["summary"].design_rank == apc["summary"].design_columns).all()
    assert {"period_rr", "cohort_rr"} <= set(apc)
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


def test_apc_window_sensitivity_compares_every_outcome(
    burden, proxy_population
):
    primary = pa.run_secondary_apc(burden, proxy_population)["summary"]
    sensitivity = pa.run_secondary_apc(
        burden, proxy_population, include_last_period=False
    )["summary"]
    comparison = pa.compare_apc_windows(primary, sensitivity)
    assert len(comparison) == 12
    assert set(comparison.measure_name) == set(pa.OUTCOMES)
    assert comparison.direction_agreement.dtype == bool
    assert comparison.comparison_note.str.contains("not an inferential test").all()


def test_cross_analysis_table_and_contradiction_audit(
    burden, proxy_population, descriptive_trends, apc_results
):
    segmented, _, _ = descriptive_trends
    endpoints = pa.endpoint_table(burden)
    decomposition = pa.run_decomposition(burden, proxy_population)
    consistency = pa.build_cross_analysis_consistency(
        burden,
        proxy_population,
        endpoints,
        segmented,
        apc_results,
        decomposition,
    )
    contradictions = pa.investigate_cross_method_contradictions(
        burden, proxy_population, consistency
    )
    assert len(consistency) == 12
    assert not consistency.methods_are_independent_replications.any()
    assert consistency.apc_net_drift_1994_2023.notna().all()
    assert set(consistency.apc_net_drift_direction) <= {"increase", "decrease", "stable"}
    assert not contradictions.empty
    assert set(contradictions.measure_name) <= set(pa.OUTCOMES)
    assert not contradictions.implementation_failure_indicated.any()
    assert contradictions.likely_explanatory_factors.str.len().gt(0).all()


def test_pandemic_sensitivity_windows_are_unambiguous(
    burden, proxy_population
):
    trend, _, _ = pa.run_segmented(burden, end_year=2019)
    assert (trend.start_year == 1990).all()
    assert (trend.end_year == 2019).all()
    apc = pa.run_secondary_apc(burden, proxy_population, include_last_period=False)
    assert (apc["summary"].start_year == 1990).all()
    assert (apc["summary"].end_year == 2019).all()
    assert set(apc["cells"].period) == {
        "1990-1994", "1995-1999", "2000-2004", "2005-2009",
        "2010-2014", "2015-2019",
    }


def test_nci_version_is_recorded_without_hard_coding_one_release():
    versions = pd.Series(["NCI Joinpoint 6.1.0", "6.1.0"])
    assert pa._consistent_software_version(versions) == "6.1.0"
    with pytest.raises(ValueError, match="mix software versions"):
        pa._consistent_software_version(pd.Series(["6.0.1", "6.1.0"]))


def test_portable_workspace_paths_use_forward_slashes():
    assert pa.portable_path(pa.DEFAULT_BURDEN) == "prepared_inputs/cause_all.csv"


def test_external_paths_are_redacted_to_portable_names():
    external = pa.ROOT.parent / "__outside_repository_test__" / "official_population.csv"
    assert pa.portable_path(external) == "external/official_population.csv"


def test_sha256_fingerprint_is_content_based(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("abc", encoding="utf-8")
    assert pa.file_sha256(source) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_export_metadata_verifies_preserved_raw_hashes(tmp_path, monkeypatch):
    monkeypatch.setattr(pa, "ROOT", tmp_path)
    raw = tmp_path / "data" / "raw" / "export.zip"
    raw.parent.mkdir(parents=True)
    raw.write_bytes(b"preserved raw export")
    metadata = {
        "export_role": "burden",
        "gbd_release": "GBD 2023",
        "retrieval_date": "2026-06-10",
        "export_id": "IHME-test-export",
        "source_url": "https://vizhub.healthdata.org/gbd-results/",
        "query_dimensions": {
            "locations": list(pa.LOCATIONS),
            "sexes": list(pa.SEXES),
            "years": "1990-2023",
            "ages": [*pa.FINE_DECOMPOSITION_AGES, pa.ALL_AGES, pa.ASR],
            "measures": list(pa.OUTCOMES),
            "metrics": ["Number", "Rate"],
        },
        "raw_files": [
            {
                "file": "data/raw/export.zip",
                "sha256": pa.file_sha256(raw),
            }
        ],
    }
    path = tmp_path / "metadata.json"
    path.write_text(json.dumps(metadata), encoding="utf-8")
    loaded = pa.load_export_metadata(path, "burden")
    assert loaded["export_id"] == "IHME-test-export"
    assert len(loaded["metadata_sha256"]) == 64

    incomplete = dict(metadata)
    incomplete["query_dimensions"] = {"years": "1990-2023"}
    path.write_text(json.dumps(incomplete), encoding="utf-8")
    with pytest.raises(ValueError, match="query_dimensions is missing"):
        pa.load_export_metadata(path, "burden")

    metadata["raw_files"][0]["sha256"] = "0" * 64
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        pa.load_export_metadata(path, "burden")


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


def test_analysis_rejects_nonempty_output_directory(tmp_path):
    (tmp_path / "stale.csv").write_text("obsolete", encoding="utf-8")
    args = SimpleNamespace(
        output_dir=tmp_path,
        burden_csv=pa.DEFAULT_BURDEN,
        population_csv=None,
        population_release="GBD 2023",
        allow_proxy_population=True,
        nci_results_csv=None,
    )
    with pytest.raises(ValueError, match="must be absent or empty"):
        pa.run(args)


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
    assert "apc_sensitivity_summary_1990_2019" in result["tables"]
    assert "apc_window_sensitivity" in result["tables"]
    assert "apc_primary_direction_agreement" in result["tables"]
    assert "decomposition_age_bin_sensitivity" in result["tables"]
    assert "segmented_specification_sensitivity" in result["tables"]
    assert "cross_analysis_consistency" in result["tables"]
    assert "cross_method_contradictions" in result["tables"]
    assert (tmp_path / "qa" / "methodological_notes.md").exists()
    assert len(result["metadata"]["burden_csv_sha256"]) == 64
    validation = json.loads((tmp_path / "qa" / "validation_summary.json").read_text(encoding="utf-8"))
    assert validation["all_age_count_reconstruction_rows"] == 408
    assert validation["all_age_count_reconstruction_within_tolerance"]
    assert (tmp_path / "qa" / "validation_summary.json").exists()

