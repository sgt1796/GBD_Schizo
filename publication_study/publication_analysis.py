from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import date
from itertools import combinations, permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.stats import f as f_dist
from scipy.stats import t as t_dist


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BURDEN = ROOT / "prepared_inputs" / "cause_all.csv"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output"

LOCATIONS = ("China", "United States of America")
SEXES = ("Female", "Male")
OUTCOMES = ("Incidence", "Prevalence", "DALYs")
ADULT_AGES = (
    "15-19 years", "20-24 years", "25-29 years", "30-34 years",
    "35-39 years", "40-44 years", "45-49 years", "50-54 years",
    "55-59 years", "60-64 years", "65-69 years", "70+ years",
)
APC_AGES = ADULT_AGES[:-1]
AGE_MID = {age: 17.0 + 5.0 * i for i, age in enumerate(APC_AGES)} | {"70+ years": 72.5}
ALL_AGES = "All ages"
ASR = "Age-standardized"
YEARS = tuple(range(1990, 2024))
ENDPOINTS = (1990, 2023)
COLORS = {"China": "#0072B2", "United States of America": "#D55E00"}
SEX_LINE = {"Female": "-", "Male": "--"}
COMPONENT_COLORS = {
    "population_growth": "#009E73",
    "population_aging": "#E69F00",
    "rate_change": "#CC79A7",
}


def clean_measure(value: str) -> str:
    return {
        "DALYs (Disability-Adjusted Life Years)": "DALYs",
        "YLDs (Years Lived with Disability)": "YLDs",
    }.get(str(value), str(value))


def load_burden(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    required = {
        "location_name", "sex_name", "age_name", "measure_name", "metric_name",
        "year", "val", "lower", "upper",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Burden file is missing columns: {sorted(missing)}")
    df = df.copy()
    df["measure_name"] = df["measure_name"].map(clean_measure)
    for col in ("year", "val", "lower", "upper"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[
        df["location_name"].isin(LOCATIONS)
        & df["sex_name"].isin(SEXES)
        & df["year"].isin(YEARS)
        & df["measure_name"].isin((*OUTCOMES, "YLDs"))
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


def load_official_population(path: Path, release: str) -> pd.DataFrame:
    df = _normalise_population_columns(pd.read_csv(path, low_memory=False))
    required = {"location_name", "sex_name", "age_name", "year", "population"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Population file is missing columns: {sorted(missing)}")
    release_text = " ".join(
        str(x) for x in ([release] + (df["gbd_release"].dropna().astype(str).unique().tolist() if "gbd_release" in df else []))
    ).lower()
    if "2023" not in release_text:
        raise ValueError("Population provenance must explicitly identify the GBD 2023 release.")
    out = df[
        df["location_name"].isin(LOCATIONS)
        & df["sex_name"].isin(SEXES)
        & df["age_name"].isin(ADULT_AGES)
        & pd.to_numeric(df["year"], errors="coerce").isin(YEARS)
    ].copy()
    out["year"] = pd.to_numeric(out["year"], errors="raise").astype(int)
    out["population"] = pd.to_numeric(out["population"], errors="raise")
    out["population_source"] = "official_GBD_2023"
    return out[["location_name", "sex_name", "age_name", "year", "population", "population_source"]]


def infer_proxy_population(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["location_name", "sex_name", "age_name", "year", "measure_name"]
    num = df[(df.metric_name == "Number") & df.age_name.isin(ADULT_AGES)][keys + ["val"]].rename(columns={"val": "number"})
    rate = df[(df.metric_name == "Rate") & df.age_name.isin(ADULT_AGES)][keys + ["val"]].rename(columns={"val": "rate"})
    merged = num.merge(rate, on=keys, validate="one_to_one")
    merged["population"] = merged["number"] / merged["rate"] * 100000.0
    out = merged.groupby(keys[:-1], as_index=False)["population"].median()
    out["population_source"] = "derived_proxy_NOT_OFFICIAL"
    return out


def validate_population(pop: pd.DataFrame) -> None:
    keys = ["location_name", "sex_name", "age_name", "year"]
    if pop.duplicated(keys).any():
        raise ValueError("Population file contains duplicate dimensional keys.")
    expected = len(LOCATIONS) * len(SEXES) * len(ADULT_AGES) * len(YEARS)
    if len(pop) != expected:
        raise ValueError(f"Population file has {len(pop)} rows; expected {expected} complete rows.")
    if not np.isfinite(pop["population"]).all() or (pop["population"] <= 0).any():
        raise ValueError("Population values must be finite and positive.")


def audit_burden(df: pd.DataFrame, pop: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    key_cols = ["location_name", "sex_name", "age_name", "measure_name", "metric_name", "year"]
    rows = []
    for outcome in OUTCOMES:
        for loc in LOCATIONS:
            for sex in SEXES:
                panel = df[(df.measure_name == outcome) & (df.location_name == loc) & (df.sex_name == sex)]
                counts = panel[(panel.metric_name == "Number") & (panel.age_name == ALL_AGES)]
                rates = panel[(panel.metric_name == "Rate") & (panel.age_name == ASR)]
                adult = panel[panel.age_name.isin(ADULT_AGES)]
                rows.append({
                    "location_name": loc, "sex_name": sex, "measure_name": outcome,
                    "all_age_count_years": counts.year.nunique(), "asr_years": rates.year.nunique(),
                    "adult_age_groups": adult.age_name.nunique(),
                    "adult_age_year_cells": adult[["age_name", "year"]].drop_duplicates().shape[0],
                    "missing_values": int(panel[["val", "lower", "upper"]].isna().sum().sum()),
                    "invalid_ui_rows": int(((panel.lower > panel.val) | (panel.val > panel.upper)).sum()),
                    "nonpositive_rows": int((panel.val <= 0).sum()),
                })
    audit = pd.DataFrame(rows)
    duplicate_audit = pd.DataFrame([{
        "rows": len(df), "duplicate_dimensional_keys": int(df.duplicated(key_cols).sum()),
        "population_rows": len(pop),
        "population_source": pop.population_source.iloc[0],
    }])

    # Reconstruction check uses official population when available; otherwise it is explicitly provisional.
    age_rates = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & df.age_name.isin(ADULT_AGES)][
        ["location_name", "sex_name", "age_name", "measure_name", "year", "val"]
    ].rename(columns={"val": "rate"})
    age_counts = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Number") & df.age_name.isin(ADULT_AGES)][
        ["location_name", "sex_name", "age_name", "measure_name", "year", "val"]
    ].rename(columns={"val": "reported_count"})
    recon = age_rates.merge(pop, on=["location_name", "sex_name", "age_name", "year"], validate="many_to_one")
    recon = recon.merge(age_counts, on=["location_name", "sex_name", "age_name", "measure_name", "year"], validate="one_to_one")
    recon["reconstructed_count"] = recon.population * recon.rate / 100000.0
    recon["absolute_error"] = recon.reconstructed_count - recon.reported_count
    recon["relative_error_pct"] = 100.0 * recon.absolute_error / recon.reported_count
    return audit, duplicate_audit, recon


def verify_yld_daly_identity(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["location_name", "sex_name", "age_name", "metric_name", "year"]
    a = df[df.measure_name == "DALYs"][keys + ["val", "lower", "upper"]]
    b = df[df.measure_name == "YLDs"][keys + ["val", "lower", "upper"]]
    m = a.merge(b, on=keys, suffixes=("_daly", "_yld"), validate="one_to_one")
    result = {"matched_cells": len(m)}
    for col in ("val", "lower", "upper"):
        delta = np.abs(m[f"{col}_daly"] - m[f"{col}_yld"])
        denom = np.maximum(np.abs(m[f"{col}_daly"]), 1e-12)
        result[f"max_abs_difference_{col}"] = float(delta.max())
        result[f"max_relative_difference_{col}"] = float((delta / denom).max())
    result["numerically_identical"] = bool(result["max_relative_difference_val"] < 1e-7)
    return pd.DataFrame([result])


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


def best_segmented_fit(years: np.ndarray, y: np.ndarray, knot_count: int, weights: np.ndarray | None = None) -> dict:
    best = None
    for knots in candidate_knots(years, knot_count):
        X = segmented_design(years, knots)
        fit = fit_linear(y, X, weights)
        if best is None or fit["rss"] < best["rss"]:
            best = fit | {"knots": knots, "X": X}
    if best is None:
        raise ValueError("No valid segmented model candidate.")
    return best


def _candidate_rss_matrix(years: np.ndarray, y_matrix: np.ndarray, knot_count: int) -> np.ndarray:
    total = np.sum(y_matrix**2, axis=0)
    best = np.full(y_matrix.shape[1], np.inf)
    for knots in candidate_knots(years, knot_count):
        q, _ = np.linalg.qr(segmented_design(years, knots), mode="reduced")
        rss = np.maximum(total - np.sum((q.T @ y_matrix) ** 2, axis=0), 0.0)
        best = np.minimum(best, rss)
    return best


def permutation_add_joinpoint(years: np.ndarray, y: np.ndarray, null_fit: dict, alt_count: int, n_perm: int, rng: np.random.Generator) -> float:
    alt = best_segmented_fit(years, y, alt_count)
    observed = (null_fit["rss"] - alt["rss"]) / max(alt["rss"], 1e-16)
    indexes = np.vstack([rng.permutation(len(y)) for _ in range(n_perm)]).T
    y_sim = null_fit["fitted"][:, None] + null_fit["residual"][indexes]
    q0, _ = np.linalg.qr(null_fit["X"], mode="reduced")
    null_rss = np.sum(y_sim**2, axis=0) - np.sum((q0.T @ y_sim) ** 2, axis=0)
    alt_rss = _candidate_rss_matrix(years, y_sim, alt_count)
    simulated = (null_rss - alt_rss) / np.maximum(alt_rss, 1e-16)
    return float((1 + np.sum(simulated >= observed)) / (n_perm + 1))


def segmented_summary(panel: pd.DataFrame, n_perm: int, rng: np.random.Generator) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    panel = panel.sort_values("year").copy()
    years = panel.year.to_numpy(int)
    y = np.log(panel.val.to_numpy(float))
    fit0 = best_segmented_fit(years, y, 0)
    p01 = permutation_add_joinpoint(years, y, fit0, 1, n_perm, rng)
    selected = fit0
    p12 = np.nan
    if p01 < 0.05:
        fit1 = best_segmented_fit(years, y, 1)
        p12 = permutation_add_joinpoint(years, y, fit1, 2, n_perm, rng)
        selected = best_segmented_fit(years, y, 2) if p12 < 0.05 else fit1

    contrast = np.zeros(len(selected["beta"]))
    contrast[1] = 1.0
    for j, knot in enumerate(selected["knots"]):
        contrast[2 + j] = (years.max() - knot) / (years.max() - years.min())
    slope = float(contrast @ selected["beta"])
    slope_se = math.sqrt(max(float(contrast @ selected["cov"] @ contrast), 0.0))
    crit = t_dist.ppf(0.975, selected["dof"])
    aapc = 100.0 * (math.exp(slope) - 1.0)
    aapc_lo = 100.0 * (math.exp(slope - crit * slope_se) - 1.0)
    aapc_hi = 100.0 * (math.exp(slope + crit * slope_se) - 1.0)
    summary = {
        "location_name": panel.location_name.iloc[0], "sex_name": panel.sex_name.iloc[0],
        "measure_name": panel.measure_name.iloc[0], "joinpoint_count": len(selected["knots"]),
        "joinpoint_years": ",".join(map(str, selected["knots"])), "permutation_p_0_vs_1": p01,
        "permutation_p_1_vs_2": p12, "aapc": aapc, "aapc_lower_model_ci": aapc_lo,
        "aapc_upper_model_ci": aapc_hi, "model_label": "independent permutation-selected segmented regression; not NCI Joinpoint",
        "permutations": n_perm,
    }
    segments = []
    breaks = (int(years.min()), *selected["knots"], int(years.max()))
    for index, (start, end) in enumerate(zip(breaks[:-1], breaks[1:])):
        c = np.zeros(len(selected["beta"])); c[1] = 1.0
        for j, knot in enumerate(selected["knots"]):
            c[2 + j] = 1.0 if start >= knot else 0.0
        seg_slope = float(c @ selected["beta"])
        seg_se = math.sqrt(max(float(c @ selected["cov"] @ c), 0.0))
        segments.append({
            "location_name": summary["location_name"], "sex_name": summary["sex_name"], "measure_name": summary["measure_name"],
            "segment_index": index + 1, "start_year": start, "end_year": end,
            "apc": 100.0 * (math.exp(seg_slope) - 1.0),
            "apc_lower_model_ci": 100.0 * (math.exp(seg_slope - crit * seg_se) - 1.0),
            "apc_upper_model_ci": 100.0 * (math.exp(seg_slope + crit * seg_se) - 1.0),
        })
    fitted = panel[["location_name", "sex_name", "measure_name", "year", "val", "lower", "upper"]].copy()
    fitted["fitted"] = np.exp(selected["fitted"])
    fitted["joinpoint_years"] = summary["joinpoint_years"]
    return summary, pd.DataFrame(segments), fitted


def run_segmented(df: pd.DataFrame, n_perm: int, seed: int, end_year: int = 2023):
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR) & (df.year <= end_year)]
    rng = np.random.default_rng(seed)
    summaries, segments, fitted = [], [], []
    for _, panel in trend.groupby(["location_name", "sex_name", "measure_name"], sort=True):
        a, b, c = segmented_summary(panel, n_perm, rng)
        summaries.append(a); segments.append(b); fitted.append(c)
    return pd.DataFrame(summaries), pd.concat(segments, ignore_index=True), pd.concat(fitted, ignore_index=True)


def weighted_trend_sensitivity(df: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    rows = []
    for _, row in primary.iterrows():
        p = trend[(trend.location_name == row.location_name) & (trend.sex_name == row.sex_name) & (trend.measure_name == row.measure_name)].sort_values("year")
        years = p.year.to_numpy(int); y = np.log(p.val.to_numpy(float))
        se = (np.log(p.upper.to_numpy(float)) - np.log(p.lower.to_numpy(float))) / (2 * 1.96)
        weights = 1.0 / np.maximum(se**2, 1e-12)
        knots = tuple(int(x) for x in str(row.joinpoint_years).split(",") if str(x).strip())
        fit = fit_linear(y, segmented_design(years, knots), weights)
        aapc = 100.0 * (math.exp((fit["fitted"][-1] - fit["fitted"][0]) / (years[-1] - years[0])) - 1.0)
        rows.append({"location_name": row.location_name, "sex_name": row.sex_name, "measure_name": row.measure_name,
                     "primary_aapc": row.aapc, "ui_weighted_fixed_knot_aapc": aapc,
                     "difference": aapc - row.aapc, "interpretation": "Sensitivity only; GBD UIs are not independent sampling SEs."})
    return pd.DataFrame(rows)


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    p = np.asarray(pvalues, float)
    order = np.argsort(p); ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted); out[order] = np.minimum(adjusted, 1.0)
    return out


def parallelism_test(a: pd.DataFrame, b: pd.DataFrame, knots: tuple[int, ...]) -> dict:
    a = a.sort_values("year"); b = b.sort_values("year")
    years = a.year.to_numpy(int)
    if not np.array_equal(years, b.year.to_numpy(int)):
        raise ValueError("Parallelism panels require identical years.")
    B = segmented_design(years, knots)
    y = np.r_[np.log(a.val.to_numpy(float)), np.log(b.val.to_numpy(float))]
    group = np.r_[np.zeros(len(years)), np.ones(len(years))]
    base = np.vstack([B, B])
    restricted = np.column_stack([base, group])
    interactions = np.vstack([np.zeros((len(years), B.shape[1] - 1)), B[:, 1:]])
    full = np.column_stack([restricted, interactions])
    rfit, ffit = fit_linear(y, restricted), fit_linear(y, full)
    df1 = full.shape[1] - restricted.shape[1]
    df2 = len(y) - full.shape[1]
    stat = ((rfit["rss"] - ffit["rss"]) / df1) / max(ffit["rss"] / df2, 1e-16)
    return {"f_statistic": stat, "df1": df1, "df2": df2, "p_value": float(f_dist.sf(stat, df1, df2))}


def pairwise_tests(df: pd.DataFrame, segmented: pd.DataFrame) -> pd.DataFrame:
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    lookup = segmented.set_index(["location_name", "sex_name", "measure_name"])
    rows = []
    for family in ("country", "sex"):
        if family == "country":
            specs = [(sex, outcome) for sex in SEXES for outcome in OUTCOMES]
            for sex, outcome in specs:
                a = trend[(trend.location_name == LOCATIONS[0]) & (trend.sex_name == sex) & (trend.measure_name == outcome)]
                b = trend[(trend.location_name == LOCATIONS[1]) & (trend.sex_name == sex) & (trend.measure_name == outcome)]
                knot_text = [lookup.loc[(loc, sex, outcome), "joinpoint_years"] for loc in LOCATIONS]
                knots = tuple(sorted({int(x) for s in knot_text for x in str(s).split(",") if x.strip()}))
                test = parallelism_test(a, b, knots)
                rows.append({"comparison_family": family, "stratum": sex, "measure_name": outcome,
                             "group_a": LOCATIONS[0], "group_b": LOCATIONS[1], "common_knot_union": ",".join(map(str, knots)), **test})
        else:
            specs = [(loc, outcome) for loc in LOCATIONS for outcome in OUTCOMES]
            for loc, outcome in specs:
                a = trend[(trend.location_name == loc) & (trend.sex_name == SEXES[0]) & (trend.measure_name == outcome)]
                b = trend[(trend.location_name == loc) & (trend.sex_name == SEXES[1]) & (trend.measure_name == outcome)]
                knot_text = [lookup.loc[(loc, sex, outcome), "joinpoint_years"] for sex in SEXES]
                knots = tuple(sorted({int(x) for s in knot_text for x in str(s).split(",") if x.strip()}))
                test = parallelism_test(a, b, knots)
                rows.append({"comparison_family": family, "stratum": loc, "measure_name": outcome,
                             "group_a": SEXES[0], "group_b": SEXES[1], "common_knot_union": ",".join(map(str, knots)), **test})
    out = pd.DataFrame(rows)
    out["q_value_bh_within_family"] = out.groupby("comparison_family")["p_value"].transform(lambda s: bh_adjust(s))
    out["parallel_trends_at_q_0_05"] = out.q_value_bh_within_family >= 0.05
    out["method_note"] = "OLS spline-basis parallelism F test; not the NCI Joinpoint pairwise comparison test."
    return out


def decompose_change(pop0: np.ndarray, pop1: np.ndarray, rate0: np.ndarray, rate1: np.ndarray) -> dict:
    n0, n1 = pop0.sum(), pop1.sum(); s0, s1 = pop0 / n0, pop1 / n1
    base = {"N": n0, "S": s0, "R": rate0}; end = {"N": n1, "S": s1, "R": rate1}
    names = {"N": "population_growth", "S": "population_aging", "R": "rate_change"}
    out = dict.fromkeys(names.values(), 0.0)
    def total(x): return float(x["N"] * np.sum(x["S"] * x["R"]))
    for order in permutations(("N", "S", "R")):
        cur = {k: np.array(v, copy=True) if isinstance(v, np.ndarray) else v for k, v in base.items()}
        for factor in order:
            before = total(cur); cur[factor] = np.array(end[factor], copy=True) if isinstance(end[factor], np.ndarray) else end[factor]
            out[names[factor]] += (total(cur) - before) / 6.0
    out["total_change"] = total(end) - total(base)
    out["component_sum"] = out["population_growth"] + out["population_aging"] + out["rate_change"]
    out["closure_error"] = out["component_sum"] - out["total_change"]
    return out


def run_decomposition(df: pd.DataFrame, pop: pd.DataFrame, windows=((1990, 2023), (2000, 2023), (2010, 2023))) -> pd.DataFrame:
    rates = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & df.age_name.isin(ADULT_AGES)][
        ["location_name", "sex_name", "measure_name", "age_name", "year", "val"]
    ].rename(columns={"val": "rate"})
    merged = rates.merge(pop, on=["location_name", "sex_name", "age_name", "year"], validate="many_to_one")
    rows = []
    for loc in LOCATIONS:
        for sex in SEXES:
            for outcome in OUTCOMES:
                p = merged[(merged.location_name == loc) & (merged.sex_name == sex) & (merged.measure_name == outcome)]
                for start, end in windows:
                    a = p[p.year == start].set_index("age_name").reindex(ADULT_AGES)
                    b = p[p.year == end].set_index("age_name").reindex(ADULT_AGES)
                    result = decompose_change(a.population.to_numpy(), b.population.to_numpy(), a.rate.to_numpy() / 100000.0, b.rate.to_numpy() / 100000.0)
                    base_count = float(np.sum(a.population * a.rate / 100000.0)); end_count = float(np.sum(b.population * b.rate / 100000.0))
                    row = {"location_name": loc, "sex_name": sex, "measure_name": outcome,
                           "start_year": start, "end_year": end, "adult_count_start_reconstructed": base_count,
                           "adult_count_end_reconstructed": end_count, "population_source": pop.population_source.iloc[0], **result}
                    for component in COMPONENT_COLORS:
                        row[f"{component}_pct_of_total"] = 100.0 * row[component] / row["total_change"] if row["total_change"] else np.nan
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


def apc_period(year: int) -> str | None:
    for start in range(1994, 2020, 5):
        if start <= year <= start + 4:
            return f"{start}-{start + 4}"
    return None


def curvature_basis(n: int) -> np.ndarray:
    # Orthogonal complement to intercept and linear trend: identifiable nonlinear curvature.
    z = np.column_stack([np.ones(n), np.arange(n, dtype=float)])
    return null_space(z.T)


def run_secondary_apc(df: pd.DataFrame, pop: pd.DataFrame, include_last_period: bool = True) -> dict[str, pd.DataFrame]:
    keys = ["location_name", "sex_name", "age_name", "year"]
    rates = df[(df.measure_name == "Incidence") & (df.metric_name == "Rate") & df.age_name.isin(APC_AGES)][keys + ["val"]].rename(columns={"val": "rate"})
    annual = rates.merge(pop[keys + ["population"]], on=keys, validate="one_to_one")
    annual = annual[annual.year.between(1994, 2023)].copy(); annual["period"] = annual.year.map(apc_period)
    if not include_last_period: annual = annual[annual.period != "2019-2023"]
    annual["age_midpoint"] = annual.age_name.map(AGE_MID)
    periods = sorted(annual.period.dropna().unique())
    p_mid = {p: np.mean([int(p[:4]), int(p[-4:])]) for p in periods}
    annual["cohort_midpoint"] = annual.apply(lambda r: p_mid[r.period] - r.age_midpoint, axis=1)
    # Recompute pooled five-year rates using person-year weights.
    temp = annual.assign(events=lambda x: x.population * x.rate / 100000.0).groupby(
        ["location_name", "sex_name", "age_name", "age_midpoint", "period", "cohort_midpoint"], as_index=False
    ).agg(events=("events", "sum"), population=("population", "sum"))
    temp["rate"] = temp.events / temp.population * 100000.0
    cells = temp
    summaries, local_rows, age_rows, period_rows, cohort_rows = [], [], [], [], []
    for (loc, sex), cell in cells.groupby(["location_name", "sex_name"]):
        ann = annual[(annual.location_name == loc) & (annual.sex_name == sex)]
        net = fit_linear(np.log(ann.rate.to_numpy()), np.column_stack([np.ones(len(ann)), ann.year - ann.year.mean(), pd.get_dummies(ann.age_name).to_numpy()[:, 1:]]), ann.population.to_numpy())
        net_slope = float(net["beta"][1]); net_se = math.sqrt(max(float(net["cov"][1, 1]), 0.0)); crit = 1.96
        summaries.append({"location_name": loc, "sex_name": sex, "measure_name": "Incidence", "start_year": int(ann.year.min()), "end_year": int(ann.year.max()),
                          "net_drift": 100*(math.exp(net_slope)-1), "net_drift_lower_model_ci": 100*(math.exp(net_slope-crit*net_se)-1),
                          "net_drift_upper_model_ci": 100*(math.exp(net_slope+crit*net_se)-1), "interpretation": "Estimable net drift; model CI does not propagate GBD posterior uncertainty."})
        for age in APC_AGES:
            aa = ann[ann.age_name == age].sort_values("year")
            fit = fit_linear(np.log(aa.rate.to_numpy()), np.column_stack([np.ones(len(aa)), aa.year - aa.year.mean()]), aa.population.to_numpy())
            slope = float(fit["beta"][1]); se = math.sqrt(max(float(fit["cov"][1,1]),0.0))
            local_rows.append({"location_name":loc,"sex_name":sex,"age_name":age,"age_midpoint":AGE_MID[age],
                               "local_drift":100*(math.exp(slope)-1),"local_drift_lower_model_ci":100*(math.exp(slope-1.96*se)-1),
                               "local_drift_upper_model_ci":100*(math.exp(slope+1.96*se)-1)})
        ages = list(APC_AGES); period_levels = periods
        cohorts = sorted(cell.cohort_midpoint.unique())
        a_index = cell.age_name.map({x:i for i,x in enumerate(ages)}).to_numpy(int)
        p_index = cell.period.map({x:i for i,x in enumerate(period_levels)}).to_numpy(int)
        c_index = cell.cohort_midpoint.map({x:i for i,x in enumerate(cohorts)}).to_numpy(int)
        A, P, C = curvature_basis(len(ages)), curvature_basis(len(period_levels)), curvature_basis(len(cohorts))
        X = np.column_stack([np.ones(len(cell)), a_index - a_index.mean(), p_index - p_index.mean(), A[a_index], P[p_index], C[c_index]])
        fit = fit_linear(np.log(cell.rate.to_numpy()), X, cell.population.to_numpy())
        beta = fit["beta"]; ia0=3; ip0=ia0+A.shape[1]; ic0=ip0+P.shape[1]
        age_dev=A@beta[ia0:ip0]; per_dev=P@beta[ip0:ic0]; coh_dev=C@beta[ic0:]
        ref_age=ages.index("40-44 years"); ref_period=len(period_levels)//2; ref_cohort=len(cohorts)//2
        for i,age in enumerate(ages): age_rows.append({"location_name":loc,"sex_name":sex,"age_name":age,"age_midpoint":AGE_MID[age],"age_rr_curvature":math.exp(age_dev[i]-age_dev[ref_age])})
        for i,p in enumerate(period_levels): period_rows.append({"location_name":loc,"sex_name":sex,"period":p,"period_midpoint":p_mid[p],"period_rr_curvature":math.exp(per_dev[i]-per_dev[ref_period])})
        for i,c in enumerate(cohorts): cohort_rows.append({"location_name":loc,"sex_name":sex,"cohort_midpoint":c,"cohort_rr_curvature":math.exp(coh_dev[i]-coh_dev[ref_cohort])})
    return {"summary":pd.DataFrame(summaries),"local_drift":pd.DataFrame(local_rows),"age_curve":pd.DataFrame(age_rows),
            "period_curvature":pd.DataFrame(period_rows),"cohort_curvature":pd.DataFrame(cohort_rows),"cells":cells}


def export_nci_inputs(df: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    trend = df[(df.measure_name.isin(OUTCOMES)) & (df.metric_name == "Rate") & (df.age_name == ASR)]
    manifest = []
    for (loc, sex, outcome), p in trend.groupby(["location_name", "sex_name", "measure_name"]):
        slug = f"{loc}_{sex}_{outcome}".lower().replace(" ", "_")
        x = p.sort_values("year")[["year", "val", "lower", "upper"]].rename(columns={"val":"rate"})
        x["log_se_from_ui_sensitivity_only"] = (np.log(x.upper)-np.log(x.lower))/(2*1.96)
        x.to_csv(out_dir / f"{slug}.csv", index=False)
        manifest.append({"file":f"{slug}.csv","location":loc,"sex":sex,"outcome":outcome})
    settings = {"software":"NCI Joinpoint 6.0.1","log_transformation":True,"min_joinpoints":0,"max_joinpoints":2,
                "minimum_observations_per_segment":4,"permutations":4499,"alpha":0.05,
                "primary_error_model":"homoscedastic","sensitivity_error_model":"log SE approximated from native GBD UI",
                "status":"requires user registration and execution in official NCI software"}
    (out_dir/"analysis_settings.json").write_text(json.dumps(settings,indent=2),encoding="utf-8")
    pd.DataFrame(manifest).to_csv(out_dir/"input_manifest.csv",index=False)


def load_normalized_nci_results(path: Path) -> dict[str, pd.DataFrame]:
    raw = pd.read_csv(path)
    required = {"analysis_type", "software_version", "location_name", "sex_name", "measure_name"}
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"Normalized NCI result file is missing columns: {sorted(missing)}")
    if not raw.software_version.astype(str).str.contains("6.0.1", regex=False).all():
        raise ValueError("Every normalized NCI result row must identify software version 6.0.1.")
    trend = raw[raw.analysis_type == "trend"].copy()
    segments = raw[raw.analysis_type == "segment"].copy()
    comparisons = raw[raw.analysis_type == "comparison"].copy()
    fitted = raw[raw.analysis_type == "fitted"].copy()
    keys = ["location_name", "sex_name", "measure_name"]
    if trend[keys].drop_duplicates().shape[0] != 12:
        raise ValueError("NCI normalized results require 12 unique trend rows.")
    if segments[keys].drop_duplicates().shape[0] != 12:
        raise ValueError("NCI normalized results require segment rows for all 12 series.")
    if len(comparisons) != 12:
        raise ValueError("NCI normalized results require 12 pairwise comparison rows.")
    if fitted[keys + ["year"]].drop_duplicates().shape[0] != 12 * 34:
        raise ValueError("NCI normalized results require 408 fitted annual values.")
    trend_out = trend.rename(columns={"aapc_lower_ci":"aapc_lower_model_ci","aapc_upper_ci":"aapc_upper_model_ci"})
    trend_out["model_label"] = "official NCI Joinpoint 6.0.1"
    segment_out = segments.rename(columns={"apc_lower_ci":"apc_lower_model_ci","apc_upper_ci":"apc_upper_model_ci"})
    pair_out = comparisons.rename(columns={"parallelism_p_value":"p_value"})
    pair_out["q_value_bh_within_family"] = pair_out.groupby("comparison_family")["p_value"].transform(lambda s: bh_adjust(s))
    pair_out["parallel_trends_at_q_0_05"] = pair_out.q_value_bh_within_family >= .05
    pair_out["method_note"] = "Official NCI Joinpoint 6.0.1 pairwise comparison."
    return {"summary":trend_out,"segments":segment_out,"pairwise":pair_out,"fitted":fitted}


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
    fig.suptitle("Permutation-selected segmented rate trends (independent implementation)",fontweight="bold")
    fig.tight_layout(rect=(0,.04,1,.96)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_age_patterns(df: pd.DataFrame, path: Path) -> None:
    p=df[(df.measure_name.isin(OUTCOMES))&(df.metric_name=="Rate")&df.age_name.isin(ADULT_AGES)&(df.year==2023)].copy()
    p["age_index"]=p.age_name.map({x:i for i,x in enumerate(ADULT_AGES)})
    fig,axes=plt.subplots(3,2,figsize=(12,11),sharex=True)
    for r,outcome in enumerate(OUTCOMES):
        for c,loc in enumerate(LOCATIONS):
            ax=axes[r,c]
            for sex in SEXES:
                s=p[(p.measure_name==outcome)&(p.location_name==loc)&(p.sex_name==sex)].sort_values("age_index")
                ax.plot(s.age_index,s.val,marker="o",ls=SEX_LINE[sex],lw=2,label=sex)
                ax.fill_between(s.age_index,s.lower,s.upper,alpha=.12)
            ax.set_title(f"{outcome}: {loc.replace('United States of America','United States')}"); ax.set_ylabel("Rate per 100,000"); ax.grid(alpha=.2)
    for ax in axes[-1]: ax.set_xticks(range(len(ADULT_AGES))); ax.set_xticklabels([x.replace(" years","") for x in ADULT_AGES],rotation=45,ha="right")
    h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=2,frameon=False)
    fig.suptitle("Age-specific schizophrenia burden in 2023",fontweight="bold"); fig.tight_layout(rect=(0,.06,1,.96)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_decomposition(decomp: pd.DataFrame, path: Path) -> None:
    d=decomp[(decomp.start_year==1990)&(decomp.end_year==2023)]
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    x=np.arange(2); width=.22
    for r,outcome in enumerate(OUTCOMES):
        for c,sex in enumerate(SEXES):
            ax=axes[r,c]; p=d[(d.measure_name==outcome)&(d.sex_name==sex)].set_index("location_name").reindex(LOCATIONS)
            for j,(component,color) in enumerate(COMPONENT_COLORS.items()): ax.bar(x+(j-1)*width,p[component],width,label=component.replace("_"," ").title(),color=color)
            ax.axhline(0,color="black",lw=.7); ax.set_xticks(x); ax.set_xticklabels(["China","United States"]); ax.set_title(f"{outcome}: {sex}"); ax.set_ylabel("Change in adult count"); ax.grid(axis="y",alpha=.2)
    h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=3,frameon=False)
    fig.suptitle("Shapley decomposition of adult burden change, 1990-2023",fontweight="bold"); fig.tight_layout(rect=(0,.05,1,.96)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


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
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,(key,xcol,ycol,title) in zip(axes.ravel(),[
        ("local_drift","age_midpoint","local_drift","Local drift by age"),
        ("age_curve","age_midpoint","age_rr_curvature","Age curvature RR"),
        ("period_curvature","period_midpoint","period_rr_curvature","Period curvature RR"),
        ("cohort_curvature","cohort_midpoint","cohort_rr_curvature","Cohort curvature RR")]):
        data=apc[key]
        for loc in LOCATIONS:
            for sex in SEXES:
                s=data[(data.location_name==loc)&(data.sex_name==sex)].sort_values(xcol)
                ax.plot(s[xcol],s[ycol],color=COLORS[loc],ls=SEX_LINE[sex],lw=2,label=f"{loc.replace('United States of America','United States')} {sex}")
        ax.axhline(0 if key=="local_drift" else 1,color="black",lw=.6); ax.set_title(title); ax.grid(alpha=.2)
    h,l=axes[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=4,frameon=False)
    fig.suptitle("Secondary incidence age-period-cohort estimable summaries",fontweight="bold"); fig.tight_layout(rect=(0,.07,1,.95)); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def write_tables(tables: dict[str,pd.DataFrame], out: Path) -> None:
    out.mkdir(parents=True,exist_ok=True)
    for name,table in tables.items(): table.to_csv(out/f"{name}.csv",index=False)
    with pd.ExcelWriter(out/"publication_tables.xlsx",engine="openpyxl") as writer:
        for name,table in tables.items(): table.to_excel(writer,sheet_name=name[:31],index=False)


def write_data_dictionary(tables: dict[str, pd.DataFrame], path: Path) -> None:
    descriptions = {
        "val":"GBD posterior mean point estimate", "lower":"Lower native 95% GBD uncertainty bound",
        "upper":"Upper native 95% GBD uncertainty bound", "aapc":"Average annual percentage change",
        "p_value":"Unadjusted model comparison P value", "q_value_bh_within_family":"Benjamini-Hochberg adjusted P value",
        "population_source":"Provenance status of population denominator", "closure_error":"Component sum minus total reconstructed change",
    }
    rows=[]
    for table_name,frame in tables.items():
        for col in frame.columns:
            rows.append({"table":table_name,"column":col,"dtype":str(frame[col].dtype),
                         "description":descriptions.get(col,col.replace("_"," ").capitalize()),
                         "uncertainty_class":"native GBD UI" if col in {"lower","upper"} else ("model-based" if "model_ci" in col else "none/not applicable")})
    pd.DataFrame(rows).to_csv(path,index=False)


def provenance_table(burden_path: Path, population_path: Path | None, population_source: str) -> pd.DataFrame:
    return pd.DataFrame([
        {"source_role":"burden estimates","gbd_release":"GBD 2023","file":str(burden_path),"retrieval_date":"not recorded in source archive","dimensions":"China; United States; 1990-2023; Female; Male; incidence; prevalence; YLD; DALY; counts and rates","status":"included"},
        {"source_role":"population denominators","gbd_release":"GBD 2023" if population_source=="official_GBD_2023" else "UNVERIFIED - PROVISIONAL ONLY","file":str(population_path or "derived from burden count/rate pairs"),"retrieval_date":date.today().isoformat(),"dimensions":"China; United States; 1990-2023; Female; Male; ages 15-19 through 70+","status":population_source},
        {"source_role":"GBD percent metric","gbd_release":"GBD 2023","file":str(burden_path),"retrieval_date":"not applicable","dimensions":"mixed measure-specific denominators","status":"excluded"},
        {"source_role":"probability of death","gbd_release":"GBD 2023","file":"prepared_inputs/GBD_1990_2023_ProbabilityOfDeath_ChinaUS_Schizophrenia.csv","retrieval_date":"not applicable","dimensions":"not a schizophrenia-specific causal outcome","status":"excluded"},
        {"source_role":"risk factor extract","gbd_release":"GBD 2023","file":"schizo/IHME-GBD_2023_DATA-5ef7a575-1.zip","retrieval_date":"not applicable","dimensions":"single sexual-violence risk branch","status":"excluded"},
    ])


def run(args) -> dict:
    out=Path(args.output_dir); tables_dir=out/"tables"; main_fig=out/"figures"/"main"; supp_fig=out/"figures"/"supplement"
    for p in (tables_dir,main_fig,supp_fig,out/"qa",out/"nci_joinpoint_inputs"): p.mkdir(parents=True,exist_ok=True)
    burden=load_burden(Path(args.burden_csv))
    if args.population_csv:
        pop=load_official_population(Path(args.population_csv),args.population_release)
    elif args.allow_proxy_population:
        pop=infer_proxy_population(burden)
    else:
        raise SystemExit("A matching official GBD 2023 --population-csv is required. Use --allow-proxy-population only for a visibly provisional build.")
    validate_population(pop)
    audit,duplicate,reconstruction=audit_burden(burden,pop)
    yld_identity=verify_yld_daly_identity(burden)
    endpoints=endpoint_table(burden); country,sex=contrast_tables(burden)
    seg,segments,fitted=run_segmented(burden,args.permutations,args.seed)
    weighted=weighted_trend_sensitivity(burden,seg)
    prepandemic,_,_=run_segmented(burden,min(args.permutations,999),args.seed+1,end_year=2019)
    pairwise=pairwise_tests(burden,seg)
    nci_valid=False
    if args.nci_results_csv:
        official=load_normalized_nci_results(Path(args.nci_results_csv))
        seg,segments,pairwise=official["summary"],official["segments"],official["pairwise"]
        fitted=official["fitted"].merge(
            burden[(burden.metric_name=="Rate")&(burden.age_name==ASR)&burden.measure_name.isin(OUTCOMES)][["location_name","sex_name","measure_name","year","val","lower","upper"]],
            on=["location_name","sex_name","measure_name","year"],how="left",validate="one_to_one")
        nci_valid=True
    decomp=run_decomposition(burden,pop); annual=chained_decomposition(burden,pop,1); fiveyear=chained_decomposition(burden,pop,5)
    apc=run_secondary_apc(burden,pop,True); apc_pre=run_secondary_apc(burden,pop,False)
    provenance=provenance_table(Path(args.burden_csv),Path(args.population_csv) if args.population_csv else None,pop.population_source.iloc[0])
    tables={"data_audit":audit,"duplicate_audit":duplicate,"population_reconstruction":reconstruction,"yld_daly_identity":yld_identity,
            "endpoint_summary":endpoints,"country_contrasts":country,"sex_contrasts":sex,"segmented_summary":seg,"segmented_segments":segments,
            "segmented_fitted":fitted,"pairwise_parallelism":pairwise,"ui_weighted_sensitivity":weighted,"prepandemic_trend_sensitivity":prepandemic,
            "decomposition":decomp,"annual_chained_decomposition":annual,"fiveyear_chained_decomposition":fiveyear,"apc_summary":apc["summary"],
            "apc_local_drift":apc["local_drift"],"apc_age_curve":apc["age_curve"],"apc_period_curvature":apc["period_curvature"],
            "apc_cohort_curvature":apc["cohort_curvature"],"apc_cells":apc["cells"],"apc_excluding_2020_2023":apc_pre["summary"],"provenance":provenance}
    write_tables(tables,tables_dir); write_data_dictionary(tables,out/"data_dictionary.csv"); export_nci_inputs(burden,out/"nci_joinpoint_inputs")
    plot_asr(burden,main_fig/"figure_1_asr_trends.png"); plot_segmented(burden,fitted,main_fig/"figure_2_segmented_trends.png")
    plot_age_patterns(burden,main_fig/"figure_3_age_patterns.png"); plot_decomposition(decomp,main_fig/"figure_4_decomposition.png")
    plot_counts(burden,supp_fig/"figure_s1_counts.png"); plot_apc(apc,supp_fig/"figure_s2_apc_incidence.png")
    metadata={"build_date":date.today().isoformat(),"population_status":pop.population_source.iloc[0],"submission_ready":bool(pop.population_source.iloc[0]=="official_GBD_2023" and nci_valid),
              "trend_status":"official NCI Joinpoint 6.0.1 imported" if nci_valid else "independent permutation implementation; official NCI results not yet imported","permutations":args.permutations,"seed":args.seed,
              "burden_csv":str(Path(args.burden_csv).resolve()),"population_csv":str(Path(args.population_csv).resolve()) if args.population_csv else None,
              "limitations":["No posterior draws","Cross-year GBD correlation unavailable","70+ is a coarse terminal age group","Ecological modeled estimates","No causal health-system inference"]}
    (out/"build_metadata.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    validation={
        "all_primary_panels_complete_34_years":bool((audit.all_age_count_years.eq(34)&audit.asr_years.eq(34)).all()),
        "duplicate_dimensional_keys":int(duplicate.duplicate_dimensional_keys.iloc[0]),
        "invalid_ui_rows":int(audit.invalid_ui_rows.sum()),
        "nonpositive_rows":int(audit.nonpositive_rows.sum()),
        "yld_daly_numerically_identical":bool(yld_identity.numerically_identical.iloc[0]),
        "population_reconstruction_p99_absolute_relative_error_pct":float(reconstruction.relative_error_pct.abs().quantile(.99)),
        "maximum_absolute_decomposition_closure_error":float(decomp.closure_error.abs().max()),
        "primary_outputs_exclude_percent_metric":True,
        "population_is_official_gbd_2023":bool(pop.population_source.iloc[0]=="official_GBD_2023"),
        "official_nci_results_imported":nci_valid,
        "submission_ready":metadata["submission_ready"],
    }
    (out/"qa"/"validation_summary.json").write_text(json.dumps(validation,indent=2),encoding="utf-8")
    return {"output":out,"tables":tables,"metadata":metadata}


def parse_args():
    p=argparse.ArgumentParser(description="Build the publication-strength China-US schizophrenia GBD 2023 study package.")
    p.add_argument("--burden-csv",default=DEFAULT_BURDEN,type=Path); p.add_argument("--population-csv",type=Path)
    p.add_argument("--population-release",default="GBD 2023"); p.add_argument("--output-dir",default=DEFAULT_OUTPUT,type=Path)
    p.add_argument("--allow-proxy-population",action="store_true",help="Provisional build only; outputs are barred from submission.")
    p.add_argument("--nci-results-csv",type=Path,help="Normalized official NCI output; presence marks trend dependency as supplied.")
    p.add_argument("--permutations",type=int,default=4499); p.add_argument("--seed",type=int,default=20230707)
    return p.parse_args()


if __name__=="__main__":
    result=run(parse_args()); print(json.dumps(result["metadata"],indent=2)); print(f"Outputs: {result['output']}")
