from __future__ import annotations

import argparse
import math
import re
from itertools import combinations, permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import apc_analysis as apc


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "GBD_data"
DEFAULT_BURDEN = DATA_DIR / "GBD_2023_schizophrenia_fine_age_China_US.csv"
DEFAULT_POPULATION = DATA_DIR / "GBD_2023_population_China_US.csv"
DEFAULT_OUTPUT = ROOT / "output"

LOCATIONS = ("China", "United States of America")
SEXES = ("Female", "Male")
OUTCOMES = ("Incidence", "Prevalence", "DALYs")
FINE_DECOMPOSITION_AGES = (
    "0-4 years", "5-9 years", "10-14 years", "15-19 years", "20-24 years",
    "25-29 years", "30-34 years", "35-39 years", "40-44 years",
    "45-49 years", "50-54 years", "55-59 years", "60-64 years",
    "65-69 years", "70-74 years", "75-79 years", "80-84 years",
    "85-89 years", "90-94 years", "95+ years",
)
APC_AGES = apc.BASE_APC_AGES
ALL_AGES = "All ages"
ASR = "Age-standardized"
YEARS = tuple(range(1990, 2024))
ENDPOINTS = (1990, 2023)
PRACTICAL_STABILITY_THRESHOLD_PCT_PER_YEAR = 0.05
PRACTICAL_STABILITY_SENSITIVITY_THRESHOLDS = (0.02, 0.05, 0.10)
COLORS = {"China": "#0072B2", "United States of America": "#D55E00"}
SEX_LINE = {"Female": "-", "Male": "--"}
COMPONENT_COLORS = {
    "population_size_change": "#009E73",
    "age_structure_change": "#E69F00",
    "age_specific_rate_change": "#CC79A7",
}

DESCRIPTIVE_INFERENCE_NOTE = (
    "Descriptive analysis of annual GBD posterior means. Primary analyses report no "
    "hypothesis tests or confidence intervals because posterior draws and cross-year/"
    "cross-stratum correlations are unavailable."
)


def clean_measure(value: str) -> str:
    return {
        "DALYs (Disability-Adjusted Life Years)": "DALYs",
        "YLDs (Years Lived with Disability)": "YLDs",
    }.get(str(value), str(value))


def load_burden(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    required = {
        "location_name", "sex_name", "age_name", "measure_name", "metric_name",
        "cause_name", "year", "val", "lower", "upper",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Burden file is missing columns: {sorted(missing)}")
    df = df.copy()
    causes = set(df["cause_name"].dropna().astype(str).str.strip())
    if causes != {"Schizophrenia"}:
        raise ValueError(f"Burden file must contain only cause_name='Schizophrenia'; found {sorted(causes)}")
    df["measure_name"] = df["measure_name"].map(clean_measure)
    for col in ("year", "val", "lower", "upper"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[
        df["location_name"].isin(LOCATIONS)
        & df["sex_name"].isin(SEXES)
        & df["year"].isin(YEARS)
        & df["measure_name"].isin(OUTCOMES)
        & df["metric_name"].isin(("Number", "Rate"))
    ].copy()
    return df


def _normalise_population_columns(df: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "location": "location_name", "sex": "sex_name", "age": "age_name",
        "population": "population", "pop": "population", "value": "population",
        "val": "population",
    }
    rename = {}
    for col in df.columns:
        key = col.strip().lower()
        if key in aliases and aliases[key] not in df.columns:
            rename[col] = aliases[key]
    return df.rename(columns=rename)


def select_decomposition_ages(available_ages: set[str]) -> tuple[str, ...]:
    missing = sorted(set(FINE_DECOMPOSITION_AGES) - available_ages)
    if missing:
        raise ValueError(f"Missing required decomposition age groups: {missing}")
    return FINE_DECOMPOSITION_AGES


def validate_required_burden_summaries(df: pd.DataFrame) -> None:
    """Require the all-age counts and ASRs used by the frozen primary analysis."""
    keys = ["location_name", "sex_name", "measure_name", "age_name", "metric_name", "year"]
    required = (
        (ALL_AGES, "Number", "all-age Number"),
        (ASR, "Rate", "age-standardized Rate"),
    )
    expected = len(LOCATIONS) * len(SEXES) * len(OUTCOMES) * len(YEARS)
    problems = []
    for age_name, metric_name, label in required:
        panel = df[
            df.measure_name.isin(OUTCOMES)
            & df.age_name.eq(age_name)
            & df.metric_name.eq(metric_name)
        ]
        duplicate_count = int(panel.duplicated(keys).sum())
        if len(panel) != expected or duplicate_count:
            problems.append(
                f"{label}: found {len(panel)} rows and {duplicate_count} duplicates; "
                f"expected {expected} unique rows"
            )
    if problems:
        raise ValueError(
            "Burden input is incomplete for the primary descriptive analysis: "
            + "; ".join(problems)
        )


def load_official_population(path: Path, release: str) -> pd.DataFrame:
    df = _normalise_population_columns(pd.read_csv(path, low_memory=False))
    required = {"location_name", "sex_name", "age_name", "year", "population", "gbd_release"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Population file is missing columns: {sorted(missing)}")
    if "2023" not in str(release).lower():
        raise ValueError("--population-release must explicitly identify the GBD 2023 release.")
    file_releases = df["gbd_release"].dropna().astype(str).str.strip()
    file_releases = file_releases[file_releases.ne("")]
    if file_releases.empty or not file_releases.str.contains("2023", regex=False).all():
        found = sorted(file_releases.unique().tolist())
        raise ValueError(f"Population file gbd_release values must consistently identify GBD 2023; found {found}")
    out = df[
        df["location_name"].isin(LOCATIONS)
        & df["sex_name"].isin(SEXES)
        & df["age_name"].isin(FINE_DECOMPOSITION_AGES)
        & pd.to_numeric(df["year"], errors="coerce").isin(YEARS)
    ].copy()
    out["year"] = pd.to_numeric(out["year"], errors="raise").astype(int)
    out["population"] = pd.to_numeric(out["population"], errors="raise")
    out["population_source"] = "official_GBD_2023"
    out = out[["location_name", "sex_name", "age_name", "year", "population", "population_source"]]
    validate_population(out, required_ages=FINE_DECOMPOSITION_AGES)
    return out


def validate_population(
    pop: pd.DataFrame, required_ages: tuple[str, ...] | None = None
) -> None:
    keys = ["location_name", "sex_name", "age_name", "year"]
    if pop.duplicated(keys).any():
        raise ValueError("Population file contains duplicate dimensional keys.")
    ages = required_ages or select_decomposition_ages(
        set(pop.age_name.dropna().astype(str))
    )
    if set(pop.age_name) != set(ages):
        unexpected = sorted(set(pop.age_name) - set(ages))
        missing = sorted(set(ages) - set(pop.age_name))
        raise ValueError(
            f"Population ages do not match the required partition; missing={missing}, "
            f"unexpected={unexpected}."
        )
    expected = len(LOCATIONS) * len(SEXES) * len(ages) * len(YEARS)
    if len(pop) != expected:
        raise ValueError(f"Population file has {len(pop)} rows; expected {expected} complete rows.")
    if not np.isfinite(pop["population"]).all() or (pop["population"] <= 0).any():
        raise ValueError("Population values must be finite and positive.")


def percent_change(start: float, end: float) -> float:
    return 100.0 * (end / start - 1.0) if start > 0 else np.nan


def endpoint_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for loc in LOCATIONS:
        for sex in SEXES:
            for outcome in OUTCOMES:
                p = df[(df.location_name == loc) & (df.sex_name == sex) & (df.measure_name == outcome)]
                for metric, age in (("Number", ALL_AGES), ("Rate", ASR)):
                    s = p[(p.metric_name == metric) & (p.age_name == age)].set_index("year")
                    a, b = s.loc[1990], s.loc[2023]
                    rows.append({
                        "location_name": loc, "sex_name": sex, "measure_name": outcome,
                        "metric_name": "All-age count" if metric == "Number" else "Age-standardized rate per 100,000",
                        "value_1990": a.val, "lower_1990": a.lower, "upper_1990": a.upper,
                        "value_2023": b.val, "lower_2023": b.lower, "upper_2023": b.upper,
                        "absolute_change": b.val - a.val,
                        "percent_change_point_estimate": percent_change(a.val, b.val),
                        "interval_dominance_change": "increase" if b.lower > a.upper else ("decrease" if b.upper < a.lower else "overlap/inconclusive"),
                        "uncertainty_note": "Native endpoint UIs; change is a point estimate, not a GBD draw-derived UI.",
                    })
    return pd.DataFrame(rows)


def contrast_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    asr = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    country_rows, sex_rows = [], []
    for year in ENDPOINTS:
        for sex in SEXES:
            for outcome in OUTCOMES:
                p = asr[(asr.year == year) & (asr.sex_name == sex) & (asr.measure_name == outcome)].set_index("location_name")
                a, b = p.loc["China"], p.loc["United States of America"]
                country_rows.append({
                    "year": year, "sex_name": sex, "measure_name": outcome,
                    "china_asr": a.val, "us_asr": b.val, "china_us_rate_ratio": a.val / b.val,
                    "absolute_rate_difference": a.val - b.val,
                    "interval_dominance": "China lower" if a.upper < b.lower else ("China higher" if a.lower > b.upper else "overlap/inconclusive"),
                    "derived_uncertainty_note": "Point-estimate ratio/difference; no posterior-draw UI.",
                })
        for loc in LOCATIONS:
            for outcome in OUTCOMES:
                p = asr[(asr.year == year) & (asr.location_name == loc) & (asr.measure_name == outcome)].set_index("sex_name")
                m, f = p.loc["Male"], p.loc["Female"]
                sex_rows.append({
                    "year": year, "location_name": loc, "measure_name": outcome,
                    "male_asr": m.val, "female_asr": f.val, "male_female_rate_ratio": m.val / f.val,
                    "absolute_rate_difference": m.val - f.val,
                    "interval_dominance": "Male higher" if m.lower > f.upper else ("Male lower" if m.upper < f.lower else "overlap/inconclusive"),
                    "derived_uncertainty_note": "Point-estimate ratio/difference; no posterior-draw UI.",
                })
    return pd.DataFrame(country_rows), pd.DataFrame(sex_rows)


def segmented_design(years: np.ndarray, knots: tuple[int, ...]) -> np.ndarray:
    x = years.astype(float) - float(years.min())
    return np.column_stack([np.ones(len(years)), x, *[np.maximum(0.0, years - k) for k in knots]])


def candidate_knots(years: np.ndarray, count: int, min_segment: int = 4) -> list[tuple[int, ...]]:
    if count == 0:
        return [()]
    candidates = [int(x) for x in years if x - years.min() >= min_segment and years.max() - x >= min_segment]
    return [k for k in combinations(candidates, count) if all(b - a >= min_segment for a, b in zip((int(years.min()), *k), (*k, int(years.max()))))]


def fit_linear(y: np.ndarray, X: np.ndarray, weights: np.ndarray | None = None) -> dict:
    if weights is None:
        Xw, yw = X, y
    else:
        root = np.sqrt(np.asarray(weights, float))
        Xw, yw = X * root[:, None], y * root
    beta = np.linalg.pinv(Xw.T @ Xw) @ (Xw.T @ yw)
    residual = y - X @ beta
    rss = float(np.sum(residual**2 if weights is None else weights * residual**2))
    dof = max(len(y) - X.shape[1], 1)
    sigma2 = rss / dof
    cov = sigma2 * np.linalg.pinv(Xw.T @ Xw)
    return {"beta": beta, "residual": residual, "rss": rss, "cov": cov, "dof": dof, "fitted": X @ beta}


def best_segmented_fit(
    years: np.ndarray,
    y: np.ndarray,
    knot_count: int,
    weights: np.ndarray | None = None,
    min_segment: int = 4,
) -> dict:
    best = None
    for knots in candidate_knots(years, knot_count, min_segment=min_segment):
        X = segmented_design(years, knots)
        fit = fit_linear(y, X, weights)
        if best is None or fit["rss"] < best["rss"]:
            best = fit | {"knots": knots, "X": X}
    if best is None:
        raise ValueError("No valid segmented model candidate.")
    return best


def select_segmented_bic(
    years: np.ndarray,
    y: np.ndarray,
    max_joinpoints: int = 2,
    min_segment: int = 4,
) -> dict:
    """Select a descriptive segmented curve using a conservative BIC heuristic.

    A breakpoint contributes both a slope-change coefficient and an estimated
    location to the parameter penalty. BIC is used only to summarize shape; it
    is not treated as formal evidence for a change point.
    """
    candidates = []
    for knot_count in range(max_joinpoints + 1):
        fit = best_segmented_fit(
            years, y, knot_count, min_segment=min_segment
        )
        n = len(y)
        parameter_count = fit["X"].shape[1] + knot_count
        bic = n * math.log(max(fit["rss"] / n, np.finfo(float).tiny)) + parameter_count * math.log(n)
        candidates.append(fit | {"bic": float(bic), "bic_parameter_count": parameter_count})
    selected = min(candidates, key=lambda item: item["bic"])
    selected["candidate_bic"] = {len(item["knots"]): item["bic"] for item in candidates}
    return selected


def fit_ar1_gls(
    y: np.ndarray, X: np.ndarray, max_iterations: int = 100, tolerance: float = 1e-10
) -> dict:
    """Fit a stationary AR(1) feasible-GLS model by Prais-Winsten iteration."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    beta = np.linalg.pinv(X) @ y
    rho = 0.0
    for _ in range(max_iterations):
        residual = y - X @ beta
        denominator = float(residual[:-1] @ residual[:-1])
        updated_rho = (
            float(residual[1:] @ residual[:-1] / denominator)
            if denominator > np.finfo(float).tiny
            else 0.0
        )
        updated_rho = float(np.clip(updated_rho, -0.95, 0.95))
        scale0 = math.sqrt(max(1.0 - updated_rho**2, np.finfo(float).eps))
        transformed_y = np.r_[scale0 * y[0], y[1:] - updated_rho * y[:-1]]
        transformed_X = np.vstack(
            [scale0 * X[0], X[1:] - updated_rho * X[:-1]]
        )
        updated_beta = np.linalg.pinv(transformed_X) @ transformed_y
        converged = bool(
            abs(updated_rho - rho) < tolerance
            and np.max(np.abs(updated_beta - beta)) < tolerance
        )
        rho, beta = updated_rho, updated_beta
        if converged:
            break
    residual = y - X @ beta
    innovations = np.r_[
        math.sqrt(max(1.0 - rho**2, np.finfo(float).eps)) * residual[0],
        residual[1:] - rho * residual[:-1],
    ]
    innovation_rss = float(innovations @ innovations)
    return {
        "beta": beta,
        "fitted": X @ beta,
        "residual": residual,
        "rho": rho,
        "innovation_rss": innovation_rss,
        "iterations_converged": converged,
    }


def select_segmented_ar1_bic(
    years: np.ndarray,
    y: np.ndarray,
    max_joinpoints: int = 2,
    min_segment: int = 4,
) -> dict:
    """Select a segmented AR(1) GLS sensitivity using Gaussian likelihood BIC."""
    candidates = []
    n = len(y)
    for knot_count in range(max_joinpoints + 1):
        best = None
        for knots in candidate_knots(years, knot_count, min_segment=min_segment):
            X = segmented_design(years, knots)
            fit = fit_ar1_gls(y, X)
            parameter_count = X.shape[1] + knot_count + 1  # includes AR(1) rho
            bic = (
                n * math.log(max(fit["innovation_rss"] / n, np.finfo(float).tiny))
                - math.log(max(1.0 - fit["rho"] ** 2, np.finfo(float).eps))
                + parameter_count * math.log(n)
            )
            item = fit | {
                "knots": knots,
                "X": X,
                "bic": float(bic),
                "bic_parameter_count": parameter_count,
            }
            if best is None or item["bic"] < best["bic"]:
                best = item
        if best is None:
            raise ValueError("No valid segmented AR(1) model candidate.")
        candidates.append(best)
    selected = min(candidates, key=lambda item: item["bic"])
    selected["candidate_bic"] = {
        len(item["knots"]): item["bic"] for item in candidates
    }
    return selected


def segmented_ar1_sensitivity(
    df: pd.DataFrame, primary_summary: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare primary segmented curves with AR(1) feasible-GLS sensitivity fits."""
    trend = df[
        df.measure_name.isin(OUTCOMES)
        & df.metric_name.eq("Rate")
        & df.age_name.eq(ASR)
    ]
    primary = primary_summary.set_index(["location_name", "sex_name", "measure_name"])
    summaries = []
    segments = []
    for keys, panel in trend.groupby(
        ["location_name", "sex_name", "measure_name"], sort=True
    ):
        panel = panel.sort_values("year")
        years = panel.year.to_numpy(int)
        selected = select_segmented_ar1_bic(years, np.log(panel.val.to_numpy(float)))
        span = float(years.max() - years.min())
        contrast = np.zeros(len(selected["beta"]))
        contrast[1] = 1.0
        for index, knot in enumerate(selected["knots"]):
            contrast[2 + index] = (years.max() - knot) / span
        aapc = 100.0 * (math.exp(float(contrast @ selected["beta"])) - 1.0)
        primary_row = primary.loc[keys]
        summaries.append(
            {
                "location_name": keys[0],
                "sex_name": keys[1],
                "measure_name": keys[2],
                "primary_joinpoint_count": int(primary_row.joinpoint_count),
                "primary_joinpoint_years": primary_row.joinpoint_years,
                "primary_aapc": float(primary_row.aapc),
                "ar1_joinpoint_count": len(selected["knots"]),
                "ar1_joinpoint_years": ",".join(map(str, selected["knots"])),
                "ar1_aapc": aapc,
                "ar1_minus_primary_aapc": aapc - float(primary_row.aapc),
                "ar1_rho": float(selected["rho"]),
                "ar1_bic": float(selected["bic"]),
                "iterations_converged": bool(selected["iterations_converged"]),
                "primary_practical_label": trend_direction(float(primary_row.aapc)),
                "ar1_practical_label": trend_direction(aapc),
                "practical_label_agreement": trend_direction(float(primary_row.aapc)) == trend_direction(aapc),
                "interpretation": (
                    "AR(1) feasible-GLS robustness analysis of posterior-mean series; "
                    "not GBD posterior uncertainty."
                ),
            }
        )
        breaks = (int(years.min()), *selected["knots"], int(years.max()))
        for segment_index, (start, end) in enumerate(zip(breaks[:-1], breaks[1:]), 1):
            slope_contrast = np.zeros(len(selected["beta"]))
            slope_contrast[1] = 1.0
            for index, knot in enumerate(selected["knots"]):
                slope_contrast[2 + index] = 1.0 if start >= knot else 0.0
            segments.append(
                {
                    "location_name": keys[0],
                    "sex_name": keys[1],
                    "measure_name": keys[2],
                    "segment_index": segment_index,
                    "start_year": start,
                    "end_year": end,
                    "ar1_segment_apc": 100.0 * (
                        math.exp(float(slope_contrast @ selected["beta"])) - 1.0
                    ),
                    "ar1_rho": float(selected["rho"]),
                    "formal_inference_performed": False,
                }
            )
    return pd.DataFrame(summaries), pd.DataFrame(segments)


def residual_diagnostics(residual: np.ndarray) -> dict[str, float | bool]:
    residual = np.asarray(residual, float)
    denominator = float(np.sum(residual**2))
    if len(residual) < 2 or denominator <= np.finfo(float).tiny:
        return {
            "lag1_residual_autocorrelation": np.nan,
            "durbin_watson": np.nan,
            "material_residual_autocorrelation": False,
        }
    left, right = residual[:-1], residual[1:]
    lag1 = float(np.corrcoef(left, right)[0, 1]) if np.std(left) > 0 and np.std(right) > 0 else np.nan
    return {
        "lag1_residual_autocorrelation": lag1,
        "durbin_watson": float(np.sum(np.diff(residual) ** 2) / denominator),
        "material_residual_autocorrelation": bool(np.isfinite(lag1) and abs(lag1) >= 0.3),
    }


def _annualized_endpoint_change(years: np.ndarray, values: np.ndarray) -> float:
    span = int(years[-1] - years[0])
    if span <= 0 or values[0] <= 0 or values[-1] <= 0:
        return np.nan
    return 100.0 * (math.exp(math.log(values[-1] / values[0]) / span) - 1.0)


def segmented_summary(panel: pd.DataFrame) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Fit a descriptive segmented curve to annual posterior means."""
    panel = panel.sort_values("year").copy()
    years = panel.year.to_numpy(int)
    values = panel.val.to_numpy(float)
    y = np.log(values)
    selected = select_segmented_bic(years, y)

    contrast = np.zeros(len(selected["beta"]))
    contrast[1] = 1.0
    for j, knot in enumerate(selected["knots"]):
        contrast[2 + j] = (years.max() - knot) / (years.max() - years.min())
    slope = float(contrast @ selected["beta"])
    aapc = 100.0 * (math.exp(slope) - 1.0)
    diagnostics = residual_diagnostics(selected["residual"])
    summary = {
        "location_name": panel.location_name.iloc[0], "sex_name": panel.sex_name.iloc[0],
        "measure_name": panel.measure_name.iloc[0], "start_year": int(years.min()), "end_year": int(years.max()),
        "joinpoint_count": len(selected["knots"]), "joinpoint_years": ",".join(map(str, selected["knots"])),
        "selection_method": "minimum BIC; breakpoint location counted in penalty; descriptive only",
        "selected_bic": selected["bic"],
        "candidate_bic_0_joinpoints": selected["candidate_bic"][0],
        "candidate_bic_1_joinpoint": selected["candidate_bic"][1],
        "candidate_bic_2_joinpoints": selected["candidate_bic"][2],
        "aapc": aapc,
        "observed_annualized_endpoint_change_pct": _annualized_endpoint_change(years, values),
        "formal_inference_performed": False,
        "model_label": "descriptive BIC-selected segmented log-linear curve; not Joinpoint",
        "inference_note": DESCRIPTIVE_INFERENCE_NOTE,
        **diagnostics,
    }
    segments = []
    breaks = (int(years.min()), *selected["knots"], int(years.max()))
    for index, (start, end) in enumerate(zip(breaks[:-1], breaks[1:])):
        c = np.zeros(len(selected["beta"]))
        c[1] = 1.0
        for j, knot in enumerate(selected["knots"]):
            c[2 + j] = 1.0 if start >= knot else 0.0
        seg_slope = float(c @ selected["beta"])
        segments.append({
            "location_name": summary["location_name"], "sex_name": summary["sex_name"], "measure_name": summary["measure_name"],
            "segment_index": index + 1, "start_year": start, "end_year": end,
            "apc": 100.0 * (math.exp(seg_slope) - 1.0),
            "formal_inference_performed": False, "inference_note": DESCRIPTIVE_INFERENCE_NOTE,
        })
    fitted = panel[["location_name", "sex_name", "measure_name", "year", "val", "lower", "upper"]].copy()
    fitted["fitted"] = np.exp(selected["fitted"])
    fitted["joinpoint_years"] = summary["joinpoint_years"]
    return summary, pd.DataFrame(segments), fitted


def segmented_specification_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    """Assess qualitative stability across reasonable descriptive curve choices."""
    specifications = (
        ("primary", 2, 4, 2023, "log_rate"),
        ("max_1_breakpoint", 1, 4, 2023, "log_rate"),
        ("max_3_breakpoints", 3, 4, 2023, "log_rate"),
        ("minimum_5_year_segments", 2, 5, 2023, "log_rate"),
        ("prepandemic_1990_2019", 2, 4, 2019, "log_rate"),
        ("linear_rate_scale", 2, 4, 2023, "rate"),
    )
    trend = df[
        df.measure_name.isin(OUTCOMES)
        & df.metric_name.eq("Rate")
        & df.age_name.eq(ASR)
    ]
    rows = []
    for (location, sex, outcome), full_panel in trend.groupby(
        ["location_name", "sex_name", "measure_name"], sort=True
    ):
        primary_direction = None
        group_rows = []
        for specification, max_joinpoints, min_segment, end_year, scale in specifications:
            panel = full_panel[full_panel.year.le(end_year)].sort_values("year")
            years = panel.year.to_numpy(int)
            values = panel.val.to_numpy(float)
            response = np.log(values) if scale == "log_rate" else values
            selected = select_segmented_bic(
                years,
                response,
                max_joinpoints=max_joinpoints,
                min_segment=min_segment,
            )
            fitted_values = (
                np.exp(selected["fitted"])
                if scale == "log_rate"
                else selected["fitted"]
            )
            annualized = _annualized_endpoint_change(years, fitted_values)
            direction = trend_direction(annualized)
            if specification == "primary":
                primary_direction = direction
            group_rows.append(
                {
                    "location_name": location,
                    "sex_name": sex,
                    "measure_name": outcome,
                    "specification": specification,
                    "model_scale": scale,
                    "start_year": int(years.min()),
                    "end_year": int(years.max()),
                    "maximum_breakpoints": max_joinpoints,
                    "minimum_segment_years": min_segment,
                    "selected_breakpoint_count": len(selected["knots"]),
                    "selected_breakpoint_years": ",".join(
                        map(str, selected["knots"])
                    ),
                    "fitted_annualized_endpoint_change_pct": annualized,
                    "trajectory_direction": direction,
                    "formal_inference_performed": False,
                    "inference_note": DESCRIPTIVE_INFERENCE_NOTE,
                }
            )
        for row in group_rows:
            row["primary_direction"] = primary_direction
            row["direction_stable_vs_primary"] = (
                row["trajectory_direction"] == primary_direction
            )
            rows.append(row)
    return pd.DataFrame(rows)


def run_segmented(df: pd.DataFrame, end_year: int = 2023):
    """Run deterministic descriptive trend summaries."""
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR) & (df.year <= end_year)]
    summaries, segments, fitted = [], [], []
    for _, panel in trend.groupby(["location_name", "sex_name", "measure_name"], sort=True):
        a, b, c = segmented_summary(panel)
        summaries.append(a)
        segments.append(b)
        fitted.append(c)
    return pd.DataFrame(summaries), pd.concat(segments, ignore_index=True), pd.concat(fitted, ignore_index=True)


def weighted_trend_sensitivity(df: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    rows = []
    for _, row in primary.iterrows():
        p = trend[(trend.location_name == row.location_name) & (trend.sex_name == row.sex_name) & (trend.measure_name == row.measure_name)].sort_values("year")
        years = p.year.to_numpy(int)
        y = np.log(p.val.to_numpy(float))
        se = (np.log(p.upper.to_numpy(float)) - np.log(p.lower.to_numpy(float))) / (2 * 1.96)
        weights = 1.0 / np.maximum(se**2, 1e-12)
        knots = tuple(int(x) for x in str(row.joinpoint_years).split(",") if str(x).strip())
        fit = fit_linear(y, segmented_design(years, knots), weights)
        aapc = 100.0 * (math.exp((fit["fitted"][-1] - fit["fitted"][0]) / (years[-1] - years[0])) - 1.0)
        rows.append({"location_name": row.location_name, "sex_name": row.sex_name, "measure_name": row.measure_name,
                     "primary_aapc": row.aapc, "ui_weighted_fixed_knot_aapc": aapc,
                     "difference": aapc - row.aapc, "interpretation": "Sensitivity only; GBD UIs are not independent sampling SEs."})
    return pd.DataFrame(rows)


def trend_direction(
    value: float, tolerance: float = PRACTICAL_STABILITY_THRESHOLD_PCT_PER_YEAR
) -> str:
    """Return a descriptive trajectory label using the practical-stability band."""
    if not np.isfinite(value):
        return "not available"
    if value > tolerance:
        return "increase"
    if value < -tolerance:
        return "decrease"
    return "practically stable"


def practical_stability_sensitivity(segmented: pd.DataFrame) -> pd.DataFrame:
    """Show how compact trajectory labels change across descriptive thresholds."""
    rows = []
    for row in segmented.itertuples(index=False):
        item = {
            "location_name": row.location_name,
            "sex_name": row.sex_name,
            "measure_name": row.measure_name,
            "aapc_pct_per_year": float(row.aapc),
        }
        for threshold in PRACTICAL_STABILITY_SENSITIVITY_THRESHOLDS:
            item[f"label_at_{threshold:.2f}_pct_per_year"] = trend_direction(
                float(row.aapc), tolerance=threshold
            )
        rows.append(item)
    out = pd.DataFrame(rows)
    out["primary_threshold_pct_per_year"] = PRACTICAL_STABILITY_THRESHOLD_PCT_PER_YEAR
    out["interpretation"] = (
        "Practical-description sensitivity only; labels are not statistical equivalence tests."
    )
    return out


def descriptive_trajectory_contrast(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Compare annual point-estimate trajectories without inferential claims."""
    a = a.sort_values("year")
    b = b.sort_values("year")
    years = a.year.to_numpy(int)
    if not np.array_equal(years, b.year.to_numpy(int)):
        raise ValueError("Trajectory contrasts require identical years.")
    values_a, values_b = a.val.to_numpy(float), b.val.to_numpy(float)
    changes_a = 100.0 * np.diff(np.log(values_a))
    changes_b = 100.0 * np.diff(np.log(values_b))
    correlation = (
        float(np.corrcoef(changes_a, changes_b)[0, 1])
        if np.std(changes_a) > 0 and np.std(changes_b) > 0
        else np.nan
    )
    endpoint_a = _annualized_endpoint_change(years, values_a)
    endpoint_b = _annualized_endpoint_change(years, values_b)
    return {
        "group_a_annualized_endpoint_change_pct": endpoint_a,
        "group_b_annualized_endpoint_change_pct": endpoint_b,
        "annualized_endpoint_change_difference_b_minus_a_pct_points": endpoint_b - endpoint_a,
        "rms_annual_log_change_difference": float(np.sqrt(np.mean((changes_b - changes_a) ** 2))),
        "annual_log_change_correlation": correlation,
        "same_annualized_endpoint_change_direction": trend_direction(endpoint_a) == trend_direction(endpoint_b),
    }


def build_trajectory_contrasts(df: pd.DataFrame) -> pd.DataFrame:
    """Create descriptive country and sex trajectory contrasts."""
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    rows = []
    for family in ("country", "sex"):
        if family == "country":
            specs = [(sex, outcome) for sex in SEXES for outcome in OUTCOMES]
            for sex, outcome in specs:
                a = trend[(trend.location_name == LOCATIONS[0]) & (trend.sex_name == sex) & (trend.measure_name == outcome)]
                b = trend[(trend.location_name == LOCATIONS[1]) & (trend.sex_name == sex) & (trend.measure_name == outcome)]
                contrast = descriptive_trajectory_contrast(a, b)
                rows.append({"comparison_family": family, "stratum": sex, "measure_name": outcome,
                             "group_a": LOCATIONS[0], "group_b": LOCATIONS[1], **contrast})
        else:
            specs = [(loc, outcome) for loc in LOCATIONS for outcome in OUTCOMES]
            for loc, outcome in specs:
                a = trend[(trend.location_name == loc) & (trend.sex_name == SEXES[0]) & (trend.measure_name == outcome)]
                b = trend[(trend.location_name == loc) & (trend.sex_name == SEXES[1]) & (trend.measure_name == outcome)]
                contrast = descriptive_trajectory_contrast(a, b)
                rows.append({"comparison_family": family, "stratum": loc, "measure_name": outcome,
                             "group_a": SEXES[0], "group_b": SEXES[1], **contrast})
    out = pd.DataFrame(rows)
    out["formal_inference_performed"] = False
    out["method_note"] = DESCRIPTIVE_INFERENCE_NOTE
    return out


def decompose_change(pop0: np.ndarray, pop1: np.ndarray, rate0: np.ndarray, rate1: np.ndarray) -> dict:
    n0, n1 = pop0.sum(), pop1.sum()
    s0, s1 = pop0 / n0, pop1 / n1
    base = {"N": n0, "S": s0, "R": rate0}
    end = {"N": n1, "S": s1, "R": rate1}
    names = {
        "N": "population_size_change",
        "S": "age_structure_change",
        "R": "age_specific_rate_change",
    }
    out = dict.fromkeys(names.values(), 0.0)
    def total(x): return float(x["N"] * np.sum(x["S"] * x["R"]))
    for order in permutations(("N", "S", "R")):
        cur = {k: np.array(v, copy=True) if isinstance(v, np.ndarray) else v for k, v in base.items()}
        for factor in order:
            before = total(cur)
            cur[factor] = (
                np.array(end[factor], copy=True)
                if isinstance(end[factor], np.ndarray)
                else end[factor]
            )
            out[names[factor]] += (total(cur) - before) / 6.0
    out["total_change"] = total(end) - total(base)
    out["component_sum"] = sum(out[name] for name in names.values())
    out["closure_error"] = out["component_sum"] - out["total_change"]
    return out


def run_decomposition(
    df: pd.DataFrame,
    pop: pd.DataFrame,
    windows=((1990, 2023), (2000, 2023), (2010, 2023)),
    ages: tuple[str, ...] | None = None,
    outcomes: tuple[str, ...] = OUTCOMES,
) -> pd.DataFrame:
    ages = ages or select_decomposition_ages(set(pop.age_name.dropna().astype(str)))
    rates = df[(df.measure_name.isin(outcomes)) & (df.metric_name == "Rate") & df.age_name.isin(ages)][
        ["location_name", "sex_name", "measure_name", "age_name", "year", "val"]
    ].rename(columns={"val": "rate"})
    merged = rates.merge(pop, on=["location_name", "sex_name", "age_name", "year"], validate="many_to_one")
    rows = []
    for loc in LOCATIONS:
        for sex in SEXES:
            for outcome in outcomes:
                p = merged[(merged.location_name == loc) & (merged.sex_name == sex) & (merged.measure_name == outcome)]
                for start, end in windows:
                    a = p[p.year == start].set_index("age_name").reindex(ages)
                    b = p[p.year == end].set_index("age_name").reindex(ages)
                    if a[["population", "rate"]].isna().any().any() or b[["population", "rate"]].isna().any().any():
                        raise ValueError(
                            f"Decomposition is missing age-specific cells for {loc}, {sex}, "
                            f"{outcome}, {start}-{end}."
                        )
                    result = decompose_change(a.population.to_numpy(), b.population.to_numpy(), a.rate.to_numpy() / 100000.0, b.rate.to_numpy() / 100000.0)
                    base_count = float(np.sum(a.population * a.rate / 100000.0))
                    end_count = float(np.sum(b.population * b.rate / 100000.0))
                    row = {"location_name": loc, "sex_name": sex, "measure_name": outcome,
                           "start_year": start, "end_year": end, "all_age_count_start_reconstructed": base_count,
                           "all_age_count_end_reconstructed": end_count, "population_source": pop.population_source.iloc[0],
                           "age_partition": (
                               "fine_5_year" if ages == FINE_DECOMPOSITION_AGES
                               else "supported_10_79" if ages == tuple(FINE_DECOMPOSITION_AGES[2:16])
                               else "custom_or_provisional"
                           ),
                           "age_group_count": len(ages), **result}
                    for component in COMPONENT_COLORS:
                        row[f"{component}_pct_of_total"] = 100.0 * row[component] / row["total_change"] if row["total_change"] else np.nan
                    rows.append(row)
    return pd.DataFrame(rows)


def incidence_supported_age_decomposition(
    df: pd.DataFrame, pop: pd.DataFrame, start_year: int = 1990, end_year: int = 2023
) -> pd.DataFrame:
    """Sensitivity restricted to incidence ages with positive source-export rates."""
    ages = tuple(FINE_DECOMPOSITION_AGES[2:16])  # 10-14 through 75-79
    burden_ages = set(df.loc[df.measure_name.eq("Incidence"), "age_name"].astype(str))
    population_ages = set(pop.age_name.astype(str))
    if not set(ages).issubset(burden_ages & population_ages):
        return pd.DataFrame(columns=[
            "location_name", "sex_name", "measure_name", "start_year", "end_year",
            "supported_age_quantity_start_reconstructed",
            "supported_age_quantity_end_reconstructed", "population_source",
            "age_partition", "quantity_scope", "age_group_count",
            "population_size_change", "age_structure_change",
            "age_specific_rate_change", "total_change", "component_sum",
            "closure_error", "population_size_change_pct_of_total",
            "age_structure_change_pct_of_total", "age_specific_rate_change_pct_of_total",
        ])
    result = run_decomposition(
        df,
        pop,
        windows=((start_year, end_year),),
        ages=ages,
        outcomes=("Incidence",),
    ).rename(
        columns={
            "all_age_count_start_reconstructed": "supported_age_quantity_start_reconstructed",
            "all_age_count_end_reconstructed": "supported_age_quantity_end_reconstructed",
        }
    )
    result.insert(
        result.columns.get_loc("age_partition") + 1,
        "quantity_scope",
        "incidence ages 10-79; sensitivity only; not all-age burden",
    )
    return result


def _collapsed_age_group(age_name: str) -> str:
    match = re.match(r"(\d+)", str(age_name))
    if not match:
        raise ValueError(f"Cannot collapse unrecognized age label {age_name!r}.")
    lower = int(match.group(1))
    if lower < 20:
        return "0-19 years"
    if lower < 40:
        return "20-39 years"
    if lower < 60:
        return "40-59 years"
    return "60+ years"


def decomposition_age_bin_sensitivity(
    df: pd.DataFrame, pop: pd.DataFrame, start_year: int = 1990, end_year: int = 2023
) -> pd.DataFrame:
    """Compare the finest available decomposition with four collapsed age groups."""
    ages = select_decomposition_ages(set(pop.age_name.dropna().astype(str)))
    baseline = run_decomposition(df, pop, windows=((start_year, end_year),))
    rates = df[
        df.measure_name.isin(OUTCOMES)
        & df.metric_name.eq("Rate")
        & df.age_name.isin(ages)
        & df.year.isin((start_year, end_year))
    ][["location_name", "sex_name", "measure_name", "age_name", "year", "val"]].rename(
        columns={"val": "rate"}
    )
    merged = rates.merge(
        pop[["location_name", "sex_name", "age_name", "year", "population"]],
        on=["location_name", "sex_name", "age_name", "year"],
        validate="many_to_one",
    )
    merged["collapsed_age_group"] = merged.age_name.map(_collapsed_age_group)
    merged["count"] = merged.population * merged.rate / 100_000.0
    collapsed = merged.groupby(
        ["location_name", "sex_name", "measure_name", "year", "collapsed_age_group"],
        as_index=False,
        observed=True,
    ).agg(population=("population", "sum"), count=("count", "sum"))
    collapsed["rate"] = collapsed["count"] / collapsed.population

    component_names = tuple(COMPONENT_COLORS)
    rows = []
    for (location, sex, outcome), panel in collapsed.groupby(
        ["location_name", "sex_name", "measure_name"], sort=True
    ):
        age_order = ("0-19 years", "20-39 years", "40-59 years", "60+ years")
        first = panel[panel.year.eq(start_year)].set_index("collapsed_age_group").reindex(age_order)
        last = panel[panel.year.eq(end_year)].set_index("collapsed_age_group").reindex(age_order)
        if first[["population", "rate"]].isna().any().any() or last[["population", "rate"]].isna().any().any():
            raise ValueError("Collapsed decomposition sensitivity has missing age-year cells.")
        coarse = decompose_change(
            first.population.to_numpy(),
            last.population.to_numpy(),
            first.rate.to_numpy(),
            last.rate.to_numpy(),
        )
        fine = baseline[
            baseline.location_name.eq(location)
            & baseline.sex_name.eq(sex)
            & baseline.measure_name.eq(outcome)
        ].iloc[0]
        fine_rank = tuple(
            sorted(component_names, key=lambda name: abs(float(fine[name])), reverse=True)
        )
        coarse_rank = tuple(
            sorted(component_names, key=lambda name: abs(float(coarse[name])), reverse=True)
        )
        total_scale = max(abs(float(fine.total_change)), np.finfo(float).eps)
        shifts = {
            name: abs(float(fine[name]) - float(coarse[name])) / total_scale * 100.0
            for name in component_names
        }
        sign_stability = {
            name: bool(np.sign(float(fine[name])) == np.sign(float(coarse[name])))
            for name in component_names
        }
        row = {
            "location_name": location,
            "sex_name": sex,
            "measure_name": outcome,
            "start_year": start_year,
            "end_year": end_year,
            "finest_age_partition": (
                "fine_5_year" if ages == FINE_DECOMPOSITION_AGES else "provisional_broad_end_groups"
            ),
            "finest_age_group_count": len(ages),
            "collapsed_age_group_count": len(age_order),
            "finest_component_rank": " > ".join(fine_rank),
            "collapsed_component_rank": " > ".join(coarse_rank),
            "component_rank_stable": fine_rank == coarse_rank,
            "maximum_component_shift_pct_of_total_change": max(shifts.values()),
        }
        for name in component_names:
            row[f"finest_{name}"] = float(fine[name])
            row[f"collapsed_{name}"] = float(coarse[name])
            row[f"{name}_sign_stable"] = sign_stability[name]
            row[f"{name}_shift_pct_of_total_change"] = shifts[name]
        row["material_age_bin_sensitivity"] = bool(
            not all(sign_stability.values())
            or fine_rank != coarse_rank
            or max(shifts.values()) >= 10.0
        )
        row["sensitivity_rule"] = (
            "Flagged if a component changes sign, absolute-magnitude rank changes, "
            "or any component shifts by at least 10% of total count change."
        )
        rows.append(row)
    return pd.DataFrame(rows)


def chained_decomposition(df: pd.DataFrame, pop: pd.DataFrame, step: int) -> pd.DataFrame:
    years = list(range(1990, 2024, step))
    if years[-1] != 2023: years.append(2023)
    pieces = []
    for start, end in zip(years[:-1], years[1:]):
        pieces.append(run_decomposition(df, pop, windows=((start, end),)))
    out = pd.concat(pieces, ignore_index=True).sort_values(["location_name", "sex_name", "measure_name", "end_year"])
    for col in (*COMPONENT_COLORS, "total_change"):
        out[f"cumulative_{col}"] = out.groupby(["location_name", "sex_name", "measure_name"])[col].cumsum()
    out["chain_step_years"] = step
    return out


def decomposition_path_sensitivity_summary(
    endpoint: pd.DataFrame, annual: pd.DataFrame, fiveyear: pd.DataFrame
) -> pd.DataFrame:
    """Compare endpoint Shapley allocations with chained allocation paths."""
    keys = ["location_name", "sex_name", "measure_name"]
    components = tuple(COMPONENT_COLORS)
    baseline = endpoint[endpoint.start_year.eq(1990) & endpoint.end_year.eq(2023)][
        keys + ["total_change", *components]
    ].copy()
    out = baseline
    for label, chained in (("annual", annual), ("fiveyear", fiveyear)):
        last = chained.sort_values("end_year").groupby(keys, as_index=False).tail(1)
        columns = {
            f"cumulative_{component}": f"{label}_chained_{component}"
            for component in components
        }
        selected = last[keys + list(columns)].rename(columns=columns)
        out = out.merge(selected, on=keys, validate="one_to_one")
        for component in components:
            out[f"{label}_{component}_shift_pct_of_total_change"] = (
                100.0
                * (out[f"{label}_chained_{component}"] - out[component]).abs()
                / out.total_change.abs().clip(lower=np.finfo(float).eps)
            )
    shift_columns = [column for column in out if column.endswith("shift_pct_of_total_change")]
    out["maximum_path_shift_pct_of_total_change"] = out[shift_columns].max(axis=1)
    out["interpretation"] = (
        "Difference in accounting allocation under endpoint, annual-chain, and five-year-chain "
        "replacement paths; not statistical uncertainty."
    )
    return out


def apc_period(year: int) -> str | None:
    """Compatibility wrapper for the primary six-period APC window."""
    return apc.PRIMARY_WINDOW.period_for_year(year)


def run_secondary_apc(
    df: pd.DataFrame,
    pop: pd.DataFrame,
    include_last_period: bool = True,
    weighting: str = "population",
) -> dict[str, pd.DataFrame]:
    """Run the primary 1994--2023 or sensitivity 1990--2019 APC model."""
    window = apc.PRIMARY_WINDOW if include_last_period else apc.SENSITIVITY_WINDOW
    return apc.run_apc(
        df, pop, window, LOCATIONS, SEXES, measures=OUTCOMES, weighting=weighting
    )


def compare_apc_weighting(
    weighted: dict[str, pd.DataFrame], unweighted: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    """Summarize the sensitivity of APC point estimates to population weighting."""
    keys = ["location_name", "sex_name", "measure_name"]
    left = weighted["summary"][keys + ["net_drift"]].rename(
        columns={"net_drift": "population_weighted_global_period_slope"}
    )
    right = unweighted["summary"][keys + ["net_drift"]].rename(
        columns={"net_drift": "equal_weight_global_period_slope"}
    )
    out = left.merge(right, on=keys, validate="one_to_one")
    out["equal_minus_population_weighted_global_period_slope"] = (
        out.equal_weight_global_period_slope - out.population_weighted_global_period_slope
    )
    curve_specs = (
        ("local_drift", "age_name", "local_drift", "maximum_absolute_age_specific_slope_difference"),
        ("age_curve", "age_name", "longitudinal_age_rr", "maximum_absolute_log_age_rate_ratio_difference"),
        ("period_rr", "period", "period_rr", "maximum_absolute_log_period_rate_ratio_difference"),
        ("cohort_rr", "cohort_midpoint", "cohort_rr", "maximum_absolute_log_cohort_rate_ratio_difference"),
    )
    for table, coordinate, value, output_name in curve_specs:
        joined = weighted[table][keys + [coordinate, value]].merge(
            unweighted[table][keys + [coordinate, value]],
            on=keys + [coordinate],
            suffixes=("_population", "_equal"),
            validate="one_to_one",
        )
        if value == "local_drift":
            joined["difference"] = (
                joined[f"{value}_equal"] - joined[f"{value}_population"]
            ).abs()
        else:
            joined["difference"] = np.abs(
                np.log(joined[f"{value}_equal"] / joined[f"{value}_population"])
            )
        maximum = joined.groupby(keys, as_index=False).difference.max().rename(
            columns={"difference": output_name}
        )
        out = out.merge(maximum, on=keys, validate="one_to_one")
    out["interpretation"] = (
        "Descriptive weighting sensitivity; differences are point-estimate changes and "
        "are not GBD posterior uncertainty."
    )
    return out


def format_apc_tables(
    results: dict[str, pd.DataFrame], prefix: str
) -> dict[str, pd.DataFrame]:
    """Relabel custom APC outputs so they are not confused with NCI estimable functions."""
    method = (
        "Custom descriptive weighted least-squares age-period-cohort summaries; "
        "not asserted equivalent to conventional NCI APC estimable functions."
    )
    role = np.where(
        results["summary"].measure_name.eq("Incidence"),
        "principal secondary APC outcome",
        "exploratory rate-surface extension",
    )
    summary = results["summary"].rename(
        columns={
            "net_drift": "global_period_slope_pct_per_year",
            "weighted_log_rate_rss": "log_rate_objective_value",
        }
    ).copy()
    summary["analysis_role"] = role
    summary["method_label"] = method
    local = results["local_drift"].rename(
        columns={"local_drift": "age_specific_slope_pct_per_year"}
    ).copy()
    local["method_label"] = method
    age = results["age_curve"].rename(
        columns={"longitudinal_age_rr": "descriptive_age_rate_ratio"}
    ).copy()
    age["method_label"] = method
    period = results["period_rr"].rename(
        columns={"period_rr": "period_curvature_rate_ratio"}
    ).copy()
    period["method_label"] = method
    cohort = results["cohort_rr"].rename(
        columns={"cohort_rr": "cohort_curvature_rate_ratio"}
    ).copy()
    cohort["method_label"] = method
    cells = results["cells"].rename(
        columns={"rate": "observed_rate", "fitted_rate": "custom_fitted_rate"}
    ).copy()
    cells["method_label"] = method
    return {
        f"{prefix}_summary": summary,
        f"{prefix}_age_specific_slopes": local,
        f"{prefix}_age_curve": age,
        f"{prefix}_period_curvature": period,
        f"{prefix}_cohort_curvature": cohort,
        f"{prefix}_cells": cells,
    }


def compare_primary_apc_directions(
    df: pd.DataFrame,
    primary: pd.DataFrame,
    apc_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Compare incidence trends with the custom APC global period slope."""
    incidence = df[
        (df.measure_name == "Incidence")
        & (df.metric_name == "Rate")
        & (df.age_name == ASR)
    ]
    endpoint_rows = []
    for (loc, sex), panel in incidence.groupby(["location_name", "sex_name"], sort=True):
        panel = panel.sort_values("year")
        endpoint_rows.append({
            "location_name": loc,
            "sex_name": sex,
            "primary_observed_annualized_endpoint_change_pct": _annualized_endpoint_change(
                panel.year.to_numpy(int), panel.val.to_numpy(float)
            ),
        })
    endpoint = pd.DataFrame(endpoint_rows)
    primary_incidence = primary[primary.measure_name == "Incidence"][
        ["location_name", "sex_name", "aapc"]
    ].rename(columns={"aapc": "primary_segmented_aapc"})
    apc_incidence = apc_summary[apc_summary.measure_name == "Incidence"][
        ["location_name", "sex_name", "net_drift"]
    ].rename(columns={"net_drift": "apc_global_period_slope"})
    out = primary_incidence.merge(endpoint, on=["location_name", "sex_name"], validate="one_to_one")
    out = out.merge(apc_incidence, on=["location_name", "sex_name"], validate="one_to_one")
    out.insert(2, "measure_name", "Incidence")
    out["primary_segmented_direction"] = out.primary_segmented_aapc.map(trend_direction)
    out["primary_observed_direction"] = out.primary_observed_annualized_endpoint_change_pct.map(trend_direction)
    out["apc_global_period_slope_label"] = out.apc_global_period_slope.map(trend_direction)
    out["apc_vs_segmented_direction_agreement"] = (
        out.apc_global_period_slope_label == out.primary_segmented_direction
    )
    out["apc_vs_observed_direction_agreement"] = (
        out.apc_global_period_slope_label == out.primary_observed_direction
    )
    out["apc_minus_segmented_aapc_pct_points"] = out.apc_global_period_slope - out.primary_segmented_aapc
    out["apc_minus_observed_annualized_endpoint_change_pct_points"] = (
        out.apc_global_period_slope - out.primary_observed_annualized_endpoint_change_pct
    )
    out["comparison_note"] = (
        "Labels apply the 0.05%/year practical-stability band to descriptive point estimates; "
        "agreement is not a significance or equivalence test."
    )
    return out


def compare_apc_windows(
    primary_summary: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Compare custom global period slopes across the two complete windows."""
    keys = ["location_name", "sex_name", "measure_name"]
    primary = primary_summary[keys + ["net_drift"]].rename(
        columns={"net_drift": "primary_global_period_slope_1994_2023"}
    )
    sensitivity = sensitivity_summary[keys + ["net_drift"]].rename(
        columns={"net_drift": "sensitivity_global_period_slope_1990_2019"}
    )
    out = primary.merge(sensitivity, on=keys, validate="one_to_one")
    out["primary_direction"] = out.primary_global_period_slope_1994_2023.map(trend_direction)
    out["sensitivity_direction"] = out.sensitivity_global_period_slope_1990_2019.map(
        trend_direction
    )
    out["direction_agreement"] = out.primary_direction.eq(out.sensitivity_direction)
    out["global_period_slope_difference_pct_points"] = (
        out.primary_global_period_slope_1994_2023 - out.sensitivity_global_period_slope_1990_2019
    )
    out["comparison_note"] = (
        "Labels compare descriptive point estimates using the practical-stability band across "
        "two six-period windows; this is not an inferential or equivalence test."
    )
    return out


def _extreme_pattern(
    frame: pd.DataFrame, value_column: str, label_column: str, digits: int = 3
) -> str:
    if frame.empty:
        return "not analyzed"
    low = frame.loc[frame[value_column].idxmin()]
    high = frame.loc[frame[value_column].idxmax()]
    return (
        f"minimum {low[value_column]:.{digits}f} at {low[label_column]}; "
        f"maximum {high[value_column]:.{digits}f} at {high[label_column]}"
    )


def build_cross_analysis_consistency(
    burden: pd.DataFrame,
    population: pd.DataFrame,
    endpoints: pd.DataFrame,
    segmented: pd.DataFrame,
    apc_results: dict[str, pd.DataFrame],
    decomposition: pd.DataFrame,
) -> pd.DataFrame:
    """Place complementary estimands side by side without treating them as replications."""
    rows = []
    for location in LOCATIONS:
        for sex in SEXES:
            for outcome in OUTCOMES:
                endpoint = endpoints[
                    endpoints.location_name.eq(location)
                    & endpoints.sex_name.eq(sex)
                    & endpoints.measure_name.eq(outcome)
                    & endpoints.metric_name.eq("Age-standardized rate per 100,000")
                ].iloc[0]
                count_endpoint = endpoints[
                    endpoints.location_name.eq(location)
                    & endpoints.sex_name.eq(sex)
                    & endpoints.measure_name.eq(outcome)
                    & endpoints.metric_name.eq("All-age count")
                ].iloc[0]
                trend = segmented[
                    segmented.location_name.eq(location)
                    & segmented.sex_name.eq(sex)
                    & segmented.measure_name.eq(outcome)
                ].iloc[0]
                decomp = decomposition[
                    decomposition.location_name.eq(location)
                    & decomposition.sex_name.eq(sex)
                    & decomposition.measure_name.eq(outcome)
                    & decomposition.start_year.eq(1990)
                    & decomposition.end_year.eq(2023)
                ].iloc[0]
                row = {
                    "location_name": location,
                    "sex_name": sex,
                    "measure_name": outcome,
                    "asr_endpoint_percent_change_1990_2023": float(
                        endpoint.percent_change_point_estimate
                    ),
                    "all_age_count_percent_change_1990_2023": float(
                        count_endpoint.percent_change_point_estimate
                    ),
                    "segmented_aapc_1990_2023": float(trend.aapc),
                    "segmented_overall_direction": trend_direction(float(trend.aapc)),
                    "segmented_breakpoint_years": trend.joinpoint_years,
                    "population_size_change": float(decomp.population_size_change),
                    "age_structure_change": float(decomp.age_structure_change),
                    "age_specific_rate_change": float(decomp.age_specific_rate_change),
                    "decomposition_total_change": float(decomp.total_change),
                    "decomposition_age_partition": decomp.age_partition,
                    "apc_global_period_slope_1994_2023": np.nan,
                    "apc_global_period_slope_label": "not analyzed",
                    "important_age_specific_slopes": "not analyzed",
                    "age_specific_slopes_include_opposing_directions": False,
                    "period_curvature_pattern": "not analyzed",
                    "cohort_curvature_pattern": "not analyzed",
                    "methods_are_independent_replications": False,
                    "uncertainty_note": DESCRIPTIVE_INFERENCE_NOTE,
                }
                panel_filter = (
                    apc_results["summary"].location_name.eq(location)
                    & apc_results["summary"].sex_name.eq(sex)
                    & apc_results["summary"].measure_name.eq(outcome)
                )
                apc_row = apc_results["summary"][panel_filter].iloc[0]
                local = apc_results["local_drift"]
                local = local[
                    local.location_name.eq(location)
                    & local.sex_name.eq(sex)
                    & local.measure_name.eq(outcome)
                ]
                period = apc_results["period_rr"]
                period = period[
                    period.location_name.eq(location)
                    & period.sex_name.eq(sex)
                    & period.measure_name.eq(outcome)
                ]
                cohort = apc_results["cohort_rr"]
                cohort = cohort[
                    cohort.location_name.eq(location)
                    & cohort.sex_name.eq(sex)
                    & cohort.measure_name.eq(outcome)
                ]
                row.update(
                    {
                        "apc_global_period_slope_1994_2023": float(apc_row.net_drift),
                        "apc_global_period_slope_label": trend_direction(float(apc_row.net_drift)),
                        "important_age_specific_slopes": _extreme_pattern(
                            local, "local_drift", "age_name"
                        ),
                        "age_specific_slopes_include_opposing_directions": bool(
                            local.local_drift.min() < 0 < local.local_drift.max()
                        ),
                        "period_curvature_pattern": _extreme_pattern(period, "period_rr", "period"),
                        "cohort_curvature_pattern": _extreme_pattern(
                            cohort, "cohort_rr", "cohort_midpoint"
                        ),
                    }
                )
                rows.append(row)
    return pd.DataFrame(rows)


def investigate_cross_method_contradictions(
    burden: pd.DataFrame,
    population: pd.DataFrame,
    consistency: pd.DataFrame,
) -> pd.DataFrame:
    """Diagnose apparent cross-method contradictions using aligned contrasts."""
    rows = []
    ages = apc.select_apc_ages(set(burden.age_name.dropna().astype(str)))
    for _, item in consistency.iterrows():
        location, sex, outcome = item.location_name, item.sex_name, item.measure_name
        asr = burden[
            burden.location_name.eq(location)
            & burden.sex_name.eq(sex)
            & burden.measure_name.eq(outcome)
            & burden.metric_name.eq("Rate")
            & burden.age_name.eq(ASR)
            & burden.year.between(1994, 2023)
        ].sort_values("year")
        asr_aligned = _annualized_endpoint_change(
            asr.year.to_numpy(int), asr.val.to_numpy(float)
        )
        keys = ["location_name", "sex_name", "age_name", "year"]
        age_rates = burden[
            burden.location_name.eq(location)
            & burden.sex_name.eq(sex)
            & burden.measure_name.eq(outcome)
            & burden.metric_name.eq("Rate")
            & burden.age_name.isin(ages)
            & burden.year.isin((1994, 2023))
        ][keys + ["val"]].rename(columns={"val": "rate"})
        age_population = population[
            population.location_name.eq(location)
            & population.sex_name.eq(sex)
            & population.age_name.isin(ages)
            & population.year.isin((1994, 2023))
        ][keys + ["population"]]
        crude = age_rates.merge(age_population, on=keys, validate="one_to_one")
        crude["events"] = crude.population * crude.rate / 100_000.0
        crude = crude.groupby("year", as_index=False).agg(
            events=("events", "sum"), population=("population", "sum")
        )
        crude["rate"] = crude.events / crude.population * 100_000.0
        crude_aligned = _annualized_endpoint_change(
            crude.year.to_numpy(int), crude.rate.to_numpy(float)
        )
        directions = {
            "segmented_1990_2023": item.segmented_overall_direction,
            "asr_endpoint_1994_2023": trend_direction(asr_aligned),
            "selected_age_crude_1994_2023": trend_direction(crude_aligned),
            "apc_global_period_slope_1994_2023": item.apc_global_period_slope_label,
        }
        disagreement = len(set(directions.values())) > 1
        local_opposition = bool(item.age_specific_slopes_include_opposing_directions)
        if not disagreement and not local_opposition:
            continue
        likely_factors = []
        if directions["segmented_1990_2023"] != directions["asr_endpoint_1994_2023"]:
            likely_factors.append("calendar window and piecewise-versus-endpoint estimand")
        if directions["asr_endpoint_1994_2023"] != directions["selected_age_crude_1994_2023"]:
            likely_factors.append("age coverage, age standardization, and population weighting")
        if directions["selected_age_crude_1994_2023"] != directions["apc_global_period_slope_1994_2023"]:
            likely_factors.append("custom APC global period slope versus crude endpoint change and model constraints")
        if local_opposition:
            likely_factors.append("opposing age-specific slopes hidden by aggregate summaries")
        rows.append(
            {
                "location_name": location,
                "sex_name": sex,
                "measure_name": outcome,
                "segmented_direction_1990_2023": directions["segmented_1990_2023"],
                "asr_endpoint_direction_1994_2023": directions[
                    "asr_endpoint_1994_2023"
                ],
                "selected_age_crude_direction_1994_2023": directions[
                    "selected_age_crude_1994_2023"
                ],
                "apc_global_period_slope_label_1994_2023": directions[
                    "apc_global_period_slope_1994_2023"
                ],
                "asr_annualized_endpoint_change_1994_2023": asr_aligned,
                "selected_age_crude_annualized_change_1994_2023": crude_aligned,
                "apc_global_period_slope_1994_2023": item.apc_global_period_slope_1994_2023,
                "opposing_age_specific_slopes": local_opposition,
                "likely_explanatory_factors": "; ".join(likely_factors),
                "implementation_failure_indicated": False,
                "interpretation": (
                    "The compared methods target different age coverage, weighting, windows, "
                    "and model functions. The discrepancy is retained as a substantive "
                    "cross-method finding unless synthetic recovery or input QA fails."
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_asr(df: pd.DataFrame, path: Path) -> None:
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    fig, axes = plt.subplots(1,3,figsize=(14,4.5),sharex=True)
    for ax,outcome in zip(axes,OUTCOMES):
        for loc in LOCATIONS:
            for sex in SEXES:
                p=trend[(trend.location_name==loc)&(trend.sex_name==sex)&(trend.measure_name==outcome)].sort_values("year")
                ax.plot(p.year,p.val,color=COLORS[loc],ls=SEX_LINE[sex],lw=2,label=f"{loc.replace('United States of America','United States')} {sex}")
                ax.fill_between(p.year,p.lower,p.upper,color=COLORS[loc],alpha=.08)
        ax.set_title(outcome); ax.set_xlabel("Year"); ax.set_ylabel("ASR per 100,000"); ax.grid(alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels(); fig.legend(handles,labels,loc="lower center",ncol=4,frameon=False)
    fig.suptitle("Sex-specific schizophrenia age-standardized rates, 1990-2023",fontweight="bold")
    fig.tight_layout(rect=(0,.12,1,.93)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_segmented(df: pd.DataFrame, fitted: pd.DataFrame, path: Path) -> None:
    fig,axes=plt.subplots(3,2,figsize=(12,11),sharex=True)
    for r,outcome in enumerate(OUTCOMES):
        for c,sex in enumerate(SEXES):
            ax=axes[r,c]
            for loc in LOCATIONS:
                p=fitted[(fitted.location_name==loc)&(fitted.sex_name==sex)&(fitted.measure_name==outcome)].sort_values("year")
                ax.scatter(p.year,p.val,s=8,color=COLORS[loc],alpha=.45)
                ax.plot(p.year,p.fitted,color=COLORS[loc],lw=2,label=loc.replace("United States of America","United States"))
            ax.set_title(f"{outcome}: {sex}"); ax.set_ylabel("ASR per 100,000"); ax.grid(alpha=.2)
    for ax in axes[-1]: ax.set_xlabel("Year")
    h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=2,frameon=False)
    fig.suptitle("Descriptive BIC-selected segmented rate trends",fontweight="bold")
    fig.tight_layout(rect=(0,.04,1,.96)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_age_patterns(df: pd.DataFrame, path: Path) -> None:
    ages = select_decomposition_ages(set(df.age_name.dropna().astype(str)))
    p=df[(df.measure_name.isin(OUTCOMES))&(df.metric_name=="Rate")&df.age_name.isin(ages)&(df.year==2023)].copy()
    p["age_index"]=p.age_name.map({x:i for i,x in enumerate(ages)})
    fig,axes=plt.subplots(3,2,figsize=(12,11),sharex=True)
    for r,outcome in enumerate(OUTCOMES):
        for c,loc in enumerate(LOCATIONS):
            ax=axes[r,c]
            for sex in SEXES:
                s=p[(p.measure_name==outcome)&(p.location_name==loc)&(p.sex_name==sex)].sort_values("age_index")
                ax.plot(s.age_index,s.val,marker="o",ls=SEX_LINE[sex],lw=2,label=sex)
                ax.fill_between(s.age_index,s.lower,s.upper,alpha=.12)
            ax.set_title(f"{outcome}: {loc.replace('United States of America','United States')}"); ax.set_ylabel("Rate per 100,000"); ax.grid(alpha=.2)
    for ax in axes[-1]: ax.set_xticks(range(len(ages))); ax.set_xticklabels([x.replace(" years","") for x in ages],rotation=45,ha="right")
    h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=2,frameon=False)
    fig.suptitle("Age-specific schizophrenia burden in 2023",fontweight="bold"); fig.tight_layout(rect=(0,.06,1,.96)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_decomposition(decomp: pd.DataFrame, path: Path) -> None:
    d=decomp[(decomp.start_year==1990)&(decomp.end_year==2023)]
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    x = np.arange(2)
    width = .22
    for r,outcome in enumerate(OUTCOMES):
        for c,sex in enumerate(SEXES):
            ax = axes[r,c]
            p = d[(d.measure_name==outcome)&(d.sex_name==sex)].set_index("location_name").reindex(LOCATIONS)
            for j,(component,color) in enumerate(COMPONENT_COLORS.items()):
                ax.bar(x+(j-1)*width,p[component],width,label=component.replace("_"," ").title(),color=color)
            ax.axhline(0,color="black",lw=.7)
            ax.set_xticks(x)
            ax.set_xticklabels(["China","United States"])
            ax.set_title(f"{outcome}: {sex}")
            ax.set_ylabel("Change in all-age quantity")
            ax.grid(axis="y",alpha=.2)
    h,l = axes[0,0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",ncol=3,frameon=False)
    fig.suptitle("Shapley decomposition of all-age burden change, 1990-2023",fontweight="bold")
    fig.tight_layout(rect=(0,.05,1,.96))
    fig.savefig(path,dpi=300,bbox_inches="tight")
    plt.close(fig)


def plot_counts(df: pd.DataFrame, path: Path) -> None:
    p=df[(df.measure_name.isin(OUTCOMES))&(df.metric_name=="Number")&(df.age_name==ALL_AGES)]
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),sharex=True)
    for ax,outcome in zip(axes,OUTCOMES):
        for loc in LOCATIONS:
            for sex in SEXES:
                s=p[(p.measure_name==outcome)&(p.location_name==loc)&(p.sex_name==sex)].sort_values("year")
                ax.plot(s.year,s.val,color=COLORS[loc],ls=SEX_LINE[sex],lw=2,label=f"{loc.replace('United States of America','United States')} {sex}")
        ax.set_title(outcome); ax.set_ylabel("All-age count"); ax.set_xlabel("Year"); ax.grid(alpha=.2)
    h,l=axes[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=4,frameon=False); fig.tight_layout(rect=(0,.12,1,1)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_apc(apc: dict[str,pd.DataFrame], path: Path) -> None:
    specifications = [
        ("local_drift", "age_midpoint", "local_drift", "Age-specific log-rate slope"),
        ("age_curve", "age_midpoint", "longitudinal_age_rr", "Descriptive age rate ratio"),
        ("period_rr", "period_midpoint", "period_rr", "Period-curvature rate ratio"),
        ("cohort_rr", "cohort_midpoint", "cohort_rr", "Cohort-curvature rate ratio"),
    ]
    fig, axes = plt.subplots(len(OUTCOMES), len(specifications), figsize=(16, 11))
    for row_index, outcome in enumerate(OUTCOMES):
        for column_index, (key, xcol, ycol, title) in enumerate(specifications):
            ax = axes[row_index, column_index]
            data = apc[key]
            for loc in LOCATIONS:
                for sex in SEXES:
                    s = data[
                        data.location_name.eq(loc)
                        & data.sex_name.eq(sex)
                        & data.measure_name.eq(outcome)
                    ].sort_values(xcol)
                    ax.plot(
                        s[xcol], s[ycol], color=COLORS[loc], ls=SEX_LINE[sex], lw=1.8,
                        label=f"{loc.replace('United States of America','United States')} {sex}",
                    )
            ax.axhline(0 if key == "local_drift" else 1, color="black", lw=.6)
            if row_index == 0:
                ax.set_title(title)
            if column_index == 0:
                ax.set_ylabel(f"{outcome}\n%/year")
            else:
                ax.set_ylabel(f"{outcome}\nrate ratio")
            ax.grid(alpha=.2)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Custom descriptive age-period-cohort summaries", fontweight="bold")
    fig.tight_layout(rect=(0, .05, 1, .96))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_tables(tables: dict[str, pd.DataFrame], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(out / f"{name}.csv", index=False)


def run(args) -> dict[str, pd.DataFrame]:
    out = Path(args.output_dir)
    tables_dir = out / "tables"
    main_fig = out / "figures" / "main"
    supp_fig = out / "figures" / "supplement"
    for directory in (tables_dir, main_fig, supp_fig):
        directory.mkdir(parents=True, exist_ok=True)

    burden = load_burden(Path(args.burden_csv))
    validate_required_burden_summaries(burden)
    population = load_official_population(Path(args.population_csv), "GBD 2023")

    endpoints = endpoint_table(burden)
    country, sex = contrast_tables(burden)
    segmented, segments, fitted = run_segmented(burden)
    ar1_sensitivity, ar1_segments = segmented_ar1_sensitivity(burden, segmented)
    stability_sensitivity = practical_stability_sensitivity(segmented)
    specification_sensitivity = segmented_specification_sensitivity(burden)
    weighted_sensitivity = weighted_trend_sensitivity(burden, segmented)
    through_2019, _, _ = run_segmented(burden, end_year=2019)
    trajectory_contrasts = build_trajectory_contrasts(burden)

    decomposition = run_decomposition(burden, population)
    annual_decomposition = chained_decomposition(burden, population, 1)
    fiveyear_decomposition = chained_decomposition(burden, population, 5)
    path_sensitivity = decomposition_path_sensitivity_summary(
        decomposition, annual_decomposition, fiveyear_decomposition
    )
    supported_age_decomposition = incidence_supported_age_decomposition(burden, population)
    age_bin_sensitivity = decomposition_age_bin_sensitivity(burden, population)

    apc_primary = run_secondary_apc(burden, population, True)
    apc_equal_weight = run_secondary_apc(burden, population, True, weighting="equal")
    apc_earlier_window = run_secondary_apc(burden, population, False)
    apc_weighting = compare_apc_weighting(apc_primary, apc_equal_weight)
    apc_agreement = compare_primary_apc_directions(burden, segmented, apc_primary["summary"])
    apc_windows = compare_apc_windows(apc_primary["summary"], apc_earlier_window["summary"])
    cross_analysis = build_cross_analysis_consistency(
        burden, population, endpoints, segmented, apc_primary, decomposition
    )
    contradictions = investigate_cross_method_contradictions(burden, population, cross_analysis)

    tables = {
        "endpoint_summary": endpoints,
        "country_contrasts": country,
        "sex_contrasts": sex,
        "segmented_summary": segmented,
        "segmented_segments": segments,
        "segmented_fitted": fitted,
        "segmented_ar1_sensitivity": ar1_sensitivity,
        "segmented_ar1_segments": ar1_segments,
        "trend_practical_stability_sensitivity": stability_sensitivity,
        "segmented_specification_sensitivity": specification_sensitivity,
        "ui_weighted_sensitivity": weighted_sensitivity,
        "trend_excluding_2020_2023": through_2019,
        "trajectory_contrasts": trajectory_contrasts,
        "decomposition": decomposition,
        "decomposition_age_bin_sensitivity": age_bin_sensitivity,
        "incidence_supported_age_decomposition": supported_age_decomposition,
        "annual_chained_decomposition": annual_decomposition,
        "fiveyear_chained_decomposition": fiveyear_decomposition,
        "decomposition_path_sensitivity": path_sensitivity,
        **format_apc_tables(apc_primary, "apc_descriptive"),
        **format_apc_tables(apc_equal_weight, "apc_unweighted"),
        **format_apc_tables(apc_earlier_window, "apc_sensitivity_1990_2019"),
        "apc_weighting_sensitivity": apc_weighting,
        "apc_window_sensitivity": apc_windows,
        "apc_primary_direction_agreement": apc_agreement,
        "cross_analysis_consistency": cross_analysis,
        "cross_method_contradictions": contradictions,
    }
    write_tables(tables, tables_dir)
    plot_asr(burden, main_fig / "figure_1_asr_trends.png")
    plot_segmented(burden, fitted, main_fig / "figure_2_segmented_trends.png")
    plot_age_patterns(burden, main_fig / "figure_3_age_patterns.png")
    plot_decomposition(decomposition, main_fig / "figure_4_decomposition.png")
    plot_counts(burden, supp_fig / "figure_s1_counts.png")
    plot_apc(apc_primary, supp_fig / "figure_s2_custom_apc_summaries.png")
    return tables


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze China-US schizophrenia burden using GBD 2023 data."
    )
    parser.add_argument("--burden-csv", type=Path, default=DEFAULT_BURDEN)
    parser.add_argument("--population-csv", type=Path, default=DEFAULT_POPULATION)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    result = run(arguments)
    print(f"Wrote {len(result)} result tables to {arguments.output_dir}")
