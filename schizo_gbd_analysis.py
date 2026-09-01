from __future__ import annotations

import argparse
import json
import math
from itertools import combinations, permutations
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, t


ROOT = Path(__file__).resolve().parent
PREP = ROOT / "prepared_inputs"
DER = ROOT / "derived"
FIG = ROOT / "figures"
EXP = DER / "plot_explanations.json"

CAUSE, CAUSE_NAME = 559, "Schizophrenia"
LOC = {6: "China", 102: "United States of America"}
LOCS = ("China", "United States of America")
SLOC = {"China": "China", "United States of America": "United States"}
SEX = {1: "Male", 2: "Female"}
SEXES = ("Female", "Male")
MEAS = {
    2: "DALYs",
    3: "YLDs",
    5: "Prevalence",
    6: "Incidence",
    27: "Probability of death",
}
MEASA = {
    "DALYs (Disability-Adjusted Life Years)": "DALYs",
    "YLDs (Years Lived with Disability)": "YLDs",
}
MET = {1: "Number", 2: "Percent", 3: "Rate", 8: "Probability of death"}

ALL, ASR = "All ages", "Age-standardized"
TMEAS = ("Incidence", "Prevalence", "YLDs", "DALYs")
MAIN = ("Incidence", "Prevalence", "DALYs")
TREND_PATTERN_MEASURES = TMEAS
APC_MEASURES = ("Incidence", "Prevalence", "DALYs")
APC_PLOT_MEASURES = ("Incidence", "DALYs")
DECOMP_INTERVAL_YEARS = (1990, 1995, 2000, 2005, 2010, 2015, 2020, 2023)

AGES = {
    8: "15-19 years",
    9: "20-24 years",
    10: "25-29 years",
    11: "30-34 years",
    12: "35-39 years",
    13: "40-44 years",
    14: "45-49 years",
    15: "50-54 years",
    16: "55-59 years",
    17: "60-64 years",
    18: "65-69 years",
    22: ALL,
    26: "70+ years",
    27: ASR,
    39: "0-14 years",
}
ADULT = (
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
    "70+ years",
)
APCAGES = ADULT[:-1]
APCREF = "40-44 years"
MID = {
    "15-19 years": 17.0,
    "20-24 years": 22.0,
    "25-29 years": 27.0,
    "30-34 years": 32.0,
    "35-39 years": 37.0,
    "40-44 years": 42.0,
    "45-49 years": 47.0,
    "50-54 years": 52.0,
    "55-59 years": 57.0,
    "60-64 years": 62.0,
    "65-69 years": 67.0,
    "70+ years": 72.5,
}

COL = {"China": "#1b9e77", "United States of America": "#d95f02"}
SEX_STYLE = {"Female": "-", "Male": "--"}
SEX_MARKER = {"Female": "o", "Male": "s"}
YCOL = {"1990": "#4c78a8", "2023": "#f58518"}

PLOT = {
    "counts": "01_trend_counts.png",
    "percent": "01b_trend_percent.png",
    "asr": "02_trend_asr.png",
    "segmented": "03_segmented_trends.png",
    "age_rate": "04_age_specific_rate_1990_vs_2023.png",
    "age_rate_difference": "05_age_specific_rate_difference.png",
    "age_number": "05b_age_specific_numbers_1990_vs_2023.png",
    "age_number_difference": "05c_age_specific_number_difference.png",
    "age_percent": "05d_age_specific_percent_1990_vs_2023.png",
    "age_percent_difference": "05e_age_specific_percent_difference.png",
    "age_number_heatmap": "05f_age_specific_number_heatmaps_1990_2023.png",
    "age_percent_heatmap": "05g_age_specific_percent_heatmaps_1990_2023.png",
    "age_rate_heatmap": "05h_age_specific_rate_heatmaps_1990_2023.png",
    "decomposition": "06_decomposition.png",
    "apc": "07_apc_rate_heatmaps.png",
    "annual_chain": "08_annual_chained_decomposition.png",
    "interval_chain": "09_interval_chained_decomposition.png",
    "apc_marginal": "11_apc_marginal_effects.png",
    "apc_cohort": "12_apc_cohort_rr.png",
}
TABLE = {
    "burden_audit": "burden_data_audit_1990_2023.csv",
    "burden_summary": "burden_summary_1990_2023.csv",
    "comparative_metrics": "comparative_metrics_1990_2023.csv",
    "segmented_summary": "segmented_summary_1990_2023.csv",
    "segmented_segments": "segmented_segment_details_1990_2023.csv",
    "segmented_fitted": "segmented_fitted_values_1990_2023.csv",
    "decomposition": "decomposition_1990_2023.csv",
    "decomposition_interval": "decomposition_interval_1990_2023.csv",
    "annual_chain_decomposition": "annual_chain_decomposition.csv",
    "annual_chain_cumulative": "annual_chain_cumulative.csv",
    "interval_chain_decomposition": "interval_chain_decomposition.csv",
    "interval_chain_cumulative": "interval_chain_cumulative.csv",
    "annual_chain_vs_endpoint": "annual_chain_vs_endpoint_1990_2023.csv",
    "interval_chain_vs_endpoint": "interval_chain_vs_endpoint_1990_2023.csv",
    "annual_chain_final_summary": "annual_chain_final_summary.csv",
    "age_specific_comparison": "age_specific_comparison_1990_2023.csv",
    "age_difference": "age_specific_china_us_difference_1990_2023.csv",
    "trend_counts": "trend_counts_source_1990_2023.csv",
    "trend_asr": "trend_asr_source_1990_2023.csv",
    "trend_percent": "trend_percent_source_1990_2023.csv",
    "apc_summary": "apc_summary_1990_2023.csv",
    "apc_age_curve": "apc_age_curve_1990_2023.csv",
    "apc_local_drift": "apc_local_drift_1990_2023.csv",
    "apc_period_rr": "apc_period_rr_1990_2023.csv",
    "apc_cohort_rr": "apc_cohort_rr_1990_2023.csv",
    "apc_cells": "apc_input_cells_1990_2023.csv",
}

STALE_OUTPUTS = (
    DER / "probability_of_death_audit_1990_2023.csv",
    DER / "probability_of_death_summary_1990_2023.csv",
    FIG / "08_probability_of_death.png",
    DER / "descriptive_apc_summary_1990_2023.csv",
    DER / "descriptive_apc_age_curve_1990_2023.csv",
    DER / "descriptive_apc_period_curve_1990_2023.csv",
    DER / "descriptive_apc_cohort_curve_1990_2023.csv",
    DER / "descriptive_apc_input_cells_1990_2023.csv",
    DER / "event_impact_summary_1990_2023.csv",
    DER / "event_impact_years_1990_2023.csv",
    FIG / "07_descriptive_apc_rate_heatmaps.png",
    FIG / "10_event_impact_panel.png",
)

EVENTS = (
    {
        "event_name": "September 11 attacks",
        "event_year": 2001,
        "locations": ("United States of America",),
        "description": "United States terrorist attacks and aftermath",
    },
    {
        "event_name": "SARS outbreak",
        "event_year": 2003,
        "locations": ("China",),
        "description": "Severe acute respiratory syndrome outbreak concentrated in China/Hong Kong region",
    },
    {
        "event_name": "Global financial crisis",
        "event_year": 2008,
        "locations": LOCS,
        "description": "Macroeconomic shock with China and United States exposure",
    },
    {
        "event_name": "COVID-19 pandemic",
        "event_year": 2020,
        "locations": LOCS,
        "description": "Pandemic-period disruption visible in 2020-2023 GBD years",
    },
)

APC_YEAR_START, APC_YEAR_END, APC_PERIOD_WIDTH = 1990, 2023, 5
ENDPOINT_YEARS = (1990, 2023)
YEARS = set(range(1990, 2024))
VALUE_COLS = ["val", "lower", "upper"]
ID_COLS = ["location_name", "sex_name", "measure_name"]


def subset(df: pd.DataFrame, **filters) -> pd.DataFrame:
    keep = pd.Series(True, index=df.index)
    for col, values in filters.items():
        if isinstance(values, (list, tuple, set, range)):
            keep &= df[col].isin(values)
        else:
            keep &= df[col].eq(values)
    return df.loc[keep]


def configure_style() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "semibold",
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.frameon": False,
            "font.size": 9,
        }
    )


def ensure_output_dirs() -> None:
    DER.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)


def clean(x) -> str:
    return "" if pd.isna(x) else " ".join(str(x).split())


def pct(start, end):
    if pd.isna(start) or pd.isna(end) or start == 0:
        return np.nan
    return 100.0 * ((end / start) - 1.0)


def sf(x):
    return np.nan if pd.isna(x) else float(x)


def age_order_index(age: str) -> int:
    order = {value: index for index, value in enumerate(("0-14 years", *ADULT, ASR, ALL))}
    return order.get(age, 999)


def save_fig(fig, key: str) -> Path:
    path = FIG / PLOT[key]
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def normalize_gbd_frame(df: pd.DataFrame) -> pd.DataFrame:
    ids = ["cause_id", "location_id", "sex_id", "measure_id", "metric_id", "age_id", "year"]
    df = df.rename(
        columns={
            "cause": "cause_id",
            "location": "location_id",
            "sex": "sex_id",
            "measure": "measure_id",
            "metric": "metric_id",
            "age": "age_id",
        }
    ).copy()
    for col in ids:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in VALUE_COLS:
        df[col] = pd.to_numeric(df[col] if col in df else np.nan, errors="coerce")
    name_maps = [
        ("cause_name", "cause_id", {CAUSE: CAUSE_NAME}),
        ("location_name", "location_id", LOC),
        ("sex_name", "sex_id", SEX),
        ("measure_name", "measure_id", MEAS),
        ("metric_name", "metric_id", MET),
        ("age_name", "age_id", AGES),
    ]
    for name, source, mapping in name_maps:
        if name not in df:
            df[name] = df[source].map(mapping)
    for col in ["cause_name", "location_name", "sex_name", "measure_name", "metric_name", "age_name"]:
        df[col] = df[col].map(clean)
    df["measure_name"] = df["measure_name"].replace(MEASA)
    df["age_name"] = df["age_name"].replace(
        {"Age standardized": ASR, "Age-standardised": ASR, "Age-standardized": ASR, "All ages": ALL}
    )
    df = df.dropna(subset=ids + ["val"]).copy()
    df[ids] = df[ids].astype(int)
    return df


def load_inputs(burden_csv=None):
    burden_path = Path(burden_csv) if burden_csv else PREP / "cause_all.csv"
    burden_df = normalize_gbd_frame(pd.read_csv(burden_path, low_memory=False))
    burden_df = subset(
        burden_df,
        cause_id=CAUSE,
        location_name=LOCS,
        sex_name=SEXES,
        measure_name=TMEAS,
        metric_name=("Number", "Percent", "Rate"),
        age_name=(ALL, ASR, *ADULT),
    )
    burden_df = burden_df[burden_df["year"].between(1990, 2023)].copy()
    return burden_df, burden_path


def renamed_values(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    return df.rename(columns={col: f"{prefix}_{'value' if col == 'val' else col}" for col in VALUE_COLS})


def select_preferred_percent_rows(df: pd.DataFrame) -> pd.DataFrame:
    percent = df[df["metric_name"] == "Percent"].copy()
    percent["basis_priority"] = np.select(
        [percent["age_name"].eq(ASR), percent["age_name"].eq(ALL)],
        [0, 1],
        default=2,
    )
    return (
        percent.sort_values(["location_name", "sex_name", "measure_name", "year", "basis_priority"])
        .groupby(["location_name", "sex_name", "measure_name", "year"], as_index=False)
        .head(1)
        .drop(columns=["basis_priority"])
    )


def compute_eapc(series: pd.DataFrame, value_col: str = "val") -> tuple[float, float, float]:
    s = series.sort_values("year")
    s = s[np.isfinite(s[value_col]) & (s[value_col] > 0)]
    if len(s) < 3:
        return np.nan, np.nan, np.nan
    reg = linregress(s["year"], np.log(s[value_col]))
    tcrit = t.ppf(0.975, len(s) - 2) if len(s) > 2 else np.nan
    return (
        100.0 * (math.exp(reg.slope) - 1.0),
        100.0 * (math.exp(reg.slope - tcrit * reg.stderr) - 1.0),
        100.0 * (math.exp(reg.slope + tcrit * reg.stderr) - 1.0),
    )


def value_at(df: pd.DataFrame, year: int, col: str = "val") -> float:
    s = df[df["year"] == year]
    return sf(s.iloc[0][col]) if not s.empty else np.nan


def text_at(df: pd.DataFrame, year: int, col: str) -> str:
    s = df[df["year"] == year]
    return clean(s.iloc[0][col]) if not s.empty else ""


def build_burden_audit_table(df: pd.DataFrame, path: Path, sexes=SEXES) -> pd.DataFrame:
    preferred_percent = select_preferred_percent_rows(df)
    adult = set(ADULT)
    rows = []
    for location in LOCS:
        for sex in sexes:
            for measure in TMEAS:
                panel = subset(df, location_name=location, sex_name=sex, measure_name=measure)
                count_years = set(subset(panel, metric_name="Number", age_name=ALL)["year"])
                rate_years = set(subset(panel, metric_name="Rate", age_name=ASR)["year"])
                percent_years = set(
                    subset(preferred_percent, location_name=location, sex_name=sex, measure_name=measure)["year"]
                )
                adult_rates = subset(panel, metric_name="Rate", year=ENDPOINT_YEARS, age_name=adult)
                adult_counts = subset(panel, metric_name="Number", year=ENDPOINT_YEARS, age_name=adult)
                rows.append(
                    {
                        "source_file": str(path),
                        "location_name": location,
                        "sex_name": sex,
                        "measure_name": measure,
                        "count_years_present": len(count_years),
                        "rate_years_present": len(rate_years),
                        "percent_years_present": len(percent_years),
                        "missing_count_years": ",".join(map(str, sorted(YEARS - count_years))),
                        "missing_rate_years": ",".join(map(str, sorted(YEARS - rate_years))),
                        "missing_percent_years": ",".join(map(str, sorted(YEARS - percent_years))),
                        "adult_rate_groups_1990": adult_rates[adult_rates["year"] == 1990]["age_name"].nunique(),
                        "adult_rate_groups_2023": adult_rates[adult_rates["year"] == 2023]["age_name"].nunique(),
                        "adult_count_groups_1990": adult_counts[adult_counts["year"] == 1990]["age_name"].nunique(),
                        "adult_count_groups_2023": adult_counts[adult_counts["year"] == 2023]["age_name"].nunique(),
                    }
                )
    return pd.DataFrame(rows)


def build_burden_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    preferred_percent = select_preferred_percent_rows(df)
    rows = []
    for location in LOCS:
        for sex in SEXES:
            for measure in TMEAS:
                panel = subset(df, location_name=location, sex_name=sex, measure_name=measure)
                count = subset(panel, metric_name="Number", age_name=ALL).sort_values("year")
                asr = subset(panel, metric_name="Rate", age_name=ASR).sort_values("year")
                percent = subset(preferred_percent, location_name=location, sex_name=sex, measure_name=measure).sort_values("year")
                eapc, eapc_lower, eapc_upper = compute_eapc(asr)
                rows.append(
                    {
                        "location_name": location,
                        "sex_name": sex,
                        "measure_name": measure,
                        "count_1990": value_at(count, 1990),
                        "count_1990_lower": value_at(count, 1990, "lower"),
                        "count_1990_upper": value_at(count, 1990, "upper"),
                        "count_2023": value_at(count, 2023),
                        "count_2023_lower": value_at(count, 2023, "lower"),
                        "count_2023_upper": value_at(count, 2023, "upper"),
                        "count_change_pct": pct(value_at(count, 1990), value_at(count, 2023)),
                        "asr_1990": value_at(asr, 1990),
                        "asr_1990_lower": value_at(asr, 1990, "lower"),
                        "asr_1990_upper": value_at(asr, 1990, "upper"),
                        "asr_2023": value_at(asr, 2023),
                        "asr_2023_lower": value_at(asr, 2023, "lower"),
                        "asr_2023_upper": value_at(asr, 2023, "upper"),
                        "asr_change_pct": pct(value_at(asr, 1990), value_at(asr, 2023)),
                        "percent_basis": text_at(percent, 1990, "age_name"),
                        "percent_1990": value_at(percent, 1990),
                        "percent_1990_lower": value_at(percent, 1990, "lower"),
                        "percent_1990_upper": value_at(percent, 1990, "upper"),
                        "percent_2023": value_at(percent, 2023),
                        "percent_2023_lower": value_at(percent, 2023, "lower"),
                        "percent_2023_upper": value_at(percent, 2023, "upper"),
                        "percent_change_pct": pct(value_at(percent, 1990), value_at(percent, 2023)),
                        "asr_eapc": eapc,
                        "asr_eapc_lower": eapc_lower,
                        "asr_eapc_upper": eapc_upper,
                    }
                )
    return pd.DataFrame(rows).sort_values(["sex_name", "measure_name", "location_name"]).reset_index(drop=True)


def event_interpretation(max_excess_pct: float, mean_excess_pct: float, slope_change: float) -> str:
    if pd.isna(max_excess_pct):
        return "insufficient data"
    if max_excess_pct >= 2.0 and (pd.isna(slope_change) or slope_change >= 0):
        return "possible upward deviation"
    if mean_excess_pct >= 2.0:
        return "above pre-event baseline"
    if max_excess_pct <= -2.0:
        return "below pre-event baseline"
    return "no clear upward deviation"


def build_event_impact_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    annual_rows = []
    summary_rows = []
    rate_df = subset(df, metric_name="Rate", age_name=ASR, measure_name=TMEAS).copy()
    for event in EVENTS:
        event_year = int(event["event_year"])
        pre_years = list(range(event_year - 5, event_year))
        window_years = list(range(event_year, min(APC_YEAR_END, event_year + 3) + 1))
        for location in event["locations"]:
            for sex in SEXES:
                for measure in TMEAS:
                    s = subset(rate_df, location_name=location, sex_name=sex, measure_name=measure).sort_values("year")
                    train = s[s["year"].isin(pre_years) & np.isfinite(s["val"]) & (s["val"] > 0)].copy()
                    observed = s[s["year"].isin(window_years)].copy()
                    if len(train) < 3 or observed.empty:
                        continue

                    slope, intercept = np.polyfit(train["year"].to_numpy(float), np.log(train["val"].to_numpy(float)), 1)
                    baseline_eapc = 100.0 * (math.exp(float(slope)) - 1.0)
                    post_eapc, _, _ = compute_eapc(observed)
                    slope_change = post_eapc - baseline_eapc if not pd.isna(post_eapc) else np.nan

                    observed = observed.copy()
                    observed["expected_value"] = np.exp(intercept + slope * observed["year"].to_numpy(float))
                    observed["observed_expected_ratio"] = np.where(
                        observed["expected_value"] > 0,
                        observed["val"] / observed["expected_value"],
                        np.nan,
                    )
                    observed["excess_pct_vs_pre_event_baseline"] = 100.0 * (observed["observed_expected_ratio"] - 1.0)
                    observed["event_name"] = event["event_name"]
                    observed["event_year"] = event_year
                    observed["event_description"] = event["description"]
                    observed["pre_event_years"] = f"{pre_years[0]}-{pre_years[-1]}"
                    observed["event_window"] = f"{window_years[0]}-{window_years[-1]}"
                    observed["baseline_eapc"] = baseline_eapc
                    observed["post_window_eapc"] = post_eapc
                    observed["post_minus_pre_eapc"] = slope_change
                    annual_rows.append(
                        observed[
                            [
                                "event_name",
                                "event_year",
                                "event_description",
                                "pre_event_years",
                                "event_window",
                                "location_name",
                                "sex_name",
                                "measure_name",
                                "year",
                                "val",
                                "expected_value",
                                "observed_expected_ratio",
                                "excess_pct_vs_pre_event_baseline",
                                "baseline_eapc",
                                "post_window_eapc",
                                "post_minus_pre_eapc",
                            ]
                        ]
                    )

                    event_year_row = observed[observed["year"] == event_year]
                    event_year_excess = (
                        sf(event_year_row.iloc[0]["excess_pct_vs_pre_event_baseline"])
                        if not event_year_row.empty
                        else np.nan
                    )
                    max_idx = observed["excess_pct_vs_pre_event_baseline"].idxmax()
                    max_row = observed.loc[max_idx]
                    mean_excess = float(observed["excess_pct_vs_pre_event_baseline"].mean())
                    max_excess = sf(max_row["excess_pct_vs_pre_event_baseline"])
                    summary_rows.append(
                        {
                            "event_name": event["event_name"],
                            "event_year": event_year,
                            "event_description": event["description"],
                            "pre_event_years": f"{pre_years[0]}-{pre_years[-1]}",
                            "event_window": f"{window_years[0]}-{window_years[-1]}",
                            "location_name": location,
                            "sex_name": sex,
                            "measure_name": measure,
                            "event_year_excess_pct": event_year_excess,
                            "max_excess_pct": max_excess,
                            "max_excess_year": int(max_row["year"]),
                            "mean_excess_pct": mean_excess,
                            "above_expected_years": int((observed["observed_expected_ratio"] > 1.0).sum()),
                            "window_year_count": int(len(observed)),
                            "baseline_eapc": baseline_eapc,
                            "post_window_eapc": post_eapc,
                            "post_minus_pre_eapc": slope_change,
                            "interpretation": event_interpretation(max_excess, mean_excess, slope_change),
                        }
                    )

    annual = pd.concat(annual_rows, ignore_index=True) if annual_rows else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary = summary.sort_values(["event_year", "event_name", "location_name", "sex_name", "measure_name"]).reset_index(drop=True)
    return summary, annual


def build_comparative_metrics_table(summary: pd.DataFrame) -> pd.DataFrame:
    metrics = {
        "All-age number": "count",
        "Age-standardized rate": "asr",
        "Measure-specific GBD percent": "percent",
    }
    rows = []
    for sex in SEXES:
        for measure in TMEAS:
            panel = subset(summary, sex_name=sex, measure_name=measure)
            china = subset(panel, location_name="China")
            us = subset(panel, location_name="United States of America")
            if china.empty or us.empty:
                continue
            c, u = china.iloc[0], us.iloc[0]
            for metric_label, prefix in metrics.items():
                c0, c1 = c[f"{prefix}_1990"], c[f"{prefix}_2023"]
                u0, u1 = u[f"{prefix}_1990"], u[f"{prefix}_2023"]
                rows.append(
                    {
                        "comparison_type": "China vs United States",
                        "location_name": "China vs United States",
                        "sex_name": sex,
                        "measure_name": measure,
                        "metric_name": metric_label,
                        "numerator": "China",
                        "denominator": "United States of America",
                        "ratio_1990": c0 / u0 if u0 else np.nan,
                        "ratio_2023": c1 / u1 if u1 else np.nan,
                        "difference_1990": c0 - u0,
                        "difference_2023": c1 - u1,
                        "ratio_change_pct": pct(c0 / u0 if u0 else np.nan, c1 / u1 if u1 else np.nan),
                    }
                )
    for location in LOCS:
        for measure in TMEAS:
            panel = subset(summary, location_name=location, measure_name=measure)
            female = subset(panel, sex_name="Female")
            male = subset(panel, sex_name="Male")
            if female.empty or male.empty:
                continue
            f, m = female.iloc[0], male.iloc[0]
            for metric_label, prefix in metrics.items():
                f0, f1 = f[f"{prefix}_1990"], f[f"{prefix}_2023"]
                m0, m1 = m[f"{prefix}_1990"], m[f"{prefix}_2023"]
                rows.append(
                    {
                        "comparison_type": "Male vs Female",
                        "location_name": location,
                        "sex_name": "Male vs Female",
                        "measure_name": measure,
                        "metric_name": metric_label,
                        "numerator": "Male",
                        "denominator": "Female",
                        "ratio_1990": m0 / f0 if f0 else np.nan,
                        "ratio_2023": m1 / f1 if f1 else np.nan,
                        "difference_1990": m0 - f0,
                        "difference_2023": m1 - f1,
                        "ratio_change_pct": pct(m0 / f0 if f0 else np.nan, m1 / f1 if f1 else np.nan),
                    }
                )
    return pd.DataFrame(rows)


def valid_knot_combo(years: np.ndarray, knots: tuple[int, ...], min_segment_years: int = 5) -> bool:
    breaks = [int(years.min()), *knots, int(years.max())]
    return all((b - a) >= min_segment_years for a, b in zip(breaks[:-1], breaks[1:]))


def knot_sets(years: np.ndarray, max_knots: int = 2, min_segment_years: int = 5):
    candidates = [int(y) for y in years[min_segment_years:-min_segment_years]]
    yield ()
    for n_knots in range(1, max_knots + 1):
        for knots in combinations(candidates, n_knots):
            if valid_knot_combo(years, knots, min_segment_years):
                yield knots


def segmented_design(years: np.ndarray, knots: tuple[int, ...]) -> np.ndarray:
    year0 = years.min()
    cols = [np.ones_like(years, dtype=float), years.astype(float) - year0]
    for knot in knots:
        cols.append(np.maximum(0.0, years.astype(float) - float(knot)))
    return np.column_stack(cols)


def fit_segmented_model(s: pd.DataFrame) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    s = s.sort_values("year").copy()
    s = s[np.isfinite(s["val"]) & (s["val"] > 0)].copy()
    years = s["year"].to_numpy()
    y = np.log(s["val"].to_numpy())
    best = None
    for knots in knot_sets(years):
        x = segmented_design(years, knots)
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        fitted_log = x @ beta
        rss = float(np.sum((y - fitted_log) ** 2))
        n = len(y)
        p = x.shape[1]
        bic = n * math.log(max(rss / n, 1e-16)) + p * math.log(n)
        if best is None or bic < best["bic"]:
            best = {"knots": knots, "beta": beta, "fitted_log": fitted_log, "rss": rss, "bic": bic}
    if best is None:
        raise ValueError("No segmented model could be fitted.")

    knots = tuple(best["knots"])
    fitted = s[["location_name", "sex_name", "measure_name", "year", "val", "lower", "upper"]].copy()
    fitted["fitted"] = np.exp(best["fitted_log"])
    fitted["residual_log"] = y - best["fitted_log"]
    fitted["model_knots"] = ",".join(map(str, knots))

    full_eapc, full_lower, full_upper = compute_eapc(s)
    aapc = 100.0 * (math.exp((best["fitted_log"][-1] - best["fitted_log"][0]) / (years[-1] - years[0])) - 1.0)
    summary = {
        "location_name": s.iloc[0]["location_name"],
        "sex_name": s.iloc[0]["sex_name"],
        "measure_name": s.iloc[0]["measure_name"],
        "metric_name": "Rate",
        "age_name": ASR,
        "start_year": int(years[0]),
        "end_year": int(years[-1]),
        "joinpoint_count": len(knots),
        "joinpoint_years": ",".join(map(str, knots)),
        "aapc": aapc,
        "single_slope_eapc": full_eapc,
        "single_slope_eapc_lower": full_lower,
        "single_slope_eapc_upper": full_upper,
        "bic": best["bic"],
        "rss": best["rss"],
    }

    beta = best["beta"]
    breaks = [int(years[0]), *knots, int(years[-1])]
    segments = []
    for idx, (start, end) in enumerate(zip(breaks[:-1], breaks[1:])):
        active = sum(1 for knot in knots if knot <= start)
        slope = float(beta[1] + np.sum(beta[2 : 2 + active]))
        segments.append(
            {
                "location_name": s.iloc[0]["location_name"],
                "sex_name": s.iloc[0]["sex_name"],
                "measure_name": s.iloc[0]["measure_name"],
                "segment_index": idx + 1,
                "segment_start": start,
                "segment_end": end,
                "segment_apc": 100.0 * (math.exp(slope) - 1.0),
            }
        )
    return summary, pd.DataFrame(segments), fitted


def run_segmented_trend_analysis(df: pd.DataFrame):
    rows, details, fitted = [], [], []
    trend = subset(df, metric_name="Rate", age_name=ASR, measure_name=TREND_PATTERN_MEASURES)
    for _, s in trend.groupby(["location_name", "sex_name", "measure_name"], sort=True):
        summary, segment_detail, fitted_values = fit_segmented_model(s)
        rows.append(summary)
        details.append(segment_detail)
        fitted.append(fitted_values)
    return (
        pd.DataFrame(rows).sort_values(["sex_name", "measure_name", "location_name"]).reset_index(drop=True),
        pd.concat(details, ignore_index=True).sort_values(["sex_name", "measure_name", "location_name", "segment_index"]),
        pd.concat(fitted, ignore_index=True).sort_values(["sex_name", "measure_name", "location_name", "year"]),
    )


def infer_population_table(df: pd.DataFrame) -> pd.DataFrame:
    numbers = subset(df, metric_name="Number", age_name=ADULT).rename(columns={"val": "number_val"})
    rates = subset(df, metric_name="Rate", age_name=ADULT).rename(columns={"val": "rate_val"})
    merged = numbers.merge(
        rates[["location_name", "sex_name", "measure_name", "age_name", "year", "rate_val"]],
        on=["location_name", "sex_name", "measure_name", "age_name", "year"],
        how="inner",
        validate="one_to_one",
    )
    merged["population_estimate"] = np.where(
        merged["rate_val"] > 0, merged["number_val"] / merged["rate_val"] * 100000.0, np.nan
    )
    return (
        merged.groupby(["location_name", "sex_name", "age_name", "year"], as_index=False)["population_estimate"]
        .median()
        .dropna()
    )


def replacement_total(total_population, age_shares, rates) -> float:
    return float(total_population * np.sum(age_shares * rates))


def decompose_change(pop_start, pop_end, rate_start, rate_end) -> dict[str, float]:
    nb, ne = float(np.sum(pop_start)), float(np.sum(pop_end))
    sb, se = pop_start / nb, pop_end / ne
    base = {"N": nb, "P": sb, "R": rate_start}
    end = {"N": ne, "P": se, "R": rate_end}
    out = {"population_growth": 0.0, "population_aging": 0.0, "rate_change": 0.0}
    name = {"N": "population_growth", "P": "population_aging", "R": "rate_change"}
    for order in permutations(["N", "P", "R"]):
        cur = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in base.items()}
        for factor in order:
            before = replacement_total(cur["N"], cur["P"], cur["R"])
            cur[factor] = end[factor].copy() if isinstance(end[factor], np.ndarray) else end[factor]
            after = replacement_total(cur["N"], cur["P"], cur["R"])
            out[name[factor]] += (after - before) / 6.0
    out["total_change"] = replacement_total(ne, se, rate_end) - replacement_total(nb, sb, rate_start)
    return out


def age_values(df: pd.DataFrame, source_col: str, output_col: str, **filters) -> pd.DataFrame:
    return subset(df, age_name=ADULT, **filters)[["age_name", source_col]].rename(columns={source_col: output_col})


def decomposition_inputs(df: pd.DataFrame, pop: pd.DataFrame, location: str, sex: str, measure: str, start_year: int, end_year: int) -> pd.DataFrame:
    burden = {"location_name": location, "sex_name": sex, "measure_name": measure}
    population = {"location_name": location, "sex_name": sex}
    parts = [
        age_values(df, "val", "count_start", **burden, metric_name="Number", year=start_year),
        age_values(df, "val", "count_end", **burden, metric_name="Number", year=end_year),
        age_values(df, "val", "rate_start", **burden, metric_name="Rate", year=start_year),
        age_values(df, "val", "rate_end", **burden, metric_name="Rate", year=end_year),
        age_values(pop, "population_estimate", "population_start", **population, year=start_year),
        age_values(pop, "population_estimate", "population_end", **population, year=end_year),
    ]
    merged = pd.DataFrame({"age_name": list(ADULT)})
    for part in parts:
        merged = merged.merge(part, on="age_name", how="inner")
    return merged


def add_contribution_percentages(out: pd.DataFrame) -> pd.DataFrame:
    for col in ["population_growth", "population_aging", "rate_change"]:
        out[f"{col}_pct_of_total"] = np.where(np.abs(out["total_change"]) > 0, 100.0 * out[col] / out["total_change"], np.nan)
    return out


def build_decomposition_table(df: pd.DataFrame, sexes=SEXES) -> pd.DataFrame:
    pop = infer_population_table(df)
    rows = []
    for location in LOCS:
        for sex in sexes:
            for measure in MAIN:
                merged = decomposition_inputs(df, pop, location, sex, measure, *ENDPOINT_YEARS)
                if len(merged) != len(ADULT):
                    continue
                contrib = decompose_change(
                    merged["population_start"].to_numpy(),
                    merged["population_end"].to_numpy(),
                    merged["rate_start"].to_numpy() / 100000.0,
                    merged["rate_end"].to_numpy() / 100000.0,
                )
                rows.append(
                    {
                        **dict(zip(ID_COLS, (location, sex, measure))),
                        "count_1990_adult": merged["count_start"].sum(),
                        "count_2023_adult": merged["count_end"].sum(),
                        **contrib,
                    }
                )
    return add_contribution_percentages(pd.DataFrame(rows)).sort_values(["sex_name", "measure_name", "location_name"]).reset_index(drop=True)


def build_interval_decomposition_table(df: pd.DataFrame, sexes=SEXES, interval_years=DECOMP_INTERVAL_YEARS) -> pd.DataFrame:
    pop = infer_population_table(df)
    rows = []
    for location in LOCS:
        for sex in sexes:
            for measure in MAIN:
                for start_year, end_year in zip(interval_years[:-1], interval_years[1:]):
                    merged = decomposition_inputs(df, pop, location, sex, measure, start_year, end_year)
                    if len(merged) != len(ADULT):
                        continue
                    contrib = decompose_change(
                        merged["population_start"].to_numpy(),
                        merged["population_end"].to_numpy(),
                        merged["rate_start"].to_numpy() / 100000.0,
                        merged["rate_end"].to_numpy() / 100000.0,
                    )
                    rows.append(
                        {
                            **dict(zip(ID_COLS, (location, sex, measure))),
                            "start_year": start_year,
                            "end_year": end_year,
                            "interval_label": f"{start_year}-{end_year}",
                            "count_start_adult": merged["count_start"].sum(),
                            "count_end_adult": merged["count_end"].sum(),
                            **contrib,
                        }
                    )
    return add_contribution_percentages(pd.DataFrame(rows)).sort_values(
        ["sex_name", "measure_name", "location_name", "start_year"]
    ).reset_index(drop=True)


def build_decomposition_interval_row(
    df: pd.DataFrame,
    pop: pd.DataFrame,
    location: str,
    sex: str,
    measure: str,
    start_year: int,
    end_year: int,
) -> dict[str, object] | None:
    merged = decomposition_inputs(df, pop, location, sex, measure, start_year, end_year)
    if len(merged) != len(ADULT):
        return None
    contrib = decompose_change(
        merged["population_start"].to_numpy(),
        merged["population_end"].to_numpy(),
        merged["rate_start"].to_numpy() / 100000.0,
        merged["rate_end"].to_numpy() / 100000.0,
    )
    return {
        **dict(zip(ID_COLS, (location, sex, measure))),
        "start_year": start_year,
        "end_year": end_year,
        "interval_label": f"{start_year}-{end_year}",
        "count_start": merged["count_start"].sum(),
        "count_end": merged["count_end"].sum(),
        **contrib,
    }


def run_decomposition_interval_set(
    df: pd.DataFrame,
    pop: pd.DataFrame,
    row_driver: pd.DataFrame,
    interval_years: tuple[int, ...] | range,
) -> pd.DataFrame:
    rows = []
    years = list(interval_years)
    for _, row in row_driver.iterrows():
        for start_year, end_year in zip(years[:-1], years[1:]):
            interval_row = build_decomposition_interval_row(
                df,
                pop,
                row["location_name"],
                row["sex_name"],
                row["measure_name"],
                int(start_year),
                int(end_year),
            )
            if interval_row is not None:
                rows.append(interval_row)
    out = pd.DataFrame(rows)
    return add_contribution_percentages(out).sort_values(
        ["sex_name", "measure_name", "location_name", "start_year"]
    ).reset_index(drop=True)


def make_cumulative_decomposition(interval_results: pd.DataFrame) -> pd.DataFrame:
    out = interval_results.sort_values(["location_name", "sex_name", "measure_name", "end_year"]).copy()
    for col in ["population_growth", "population_aging", "rate_change", "total_change"]:
        out[f"cumulative_{col}"] = out.groupby(["location_name", "sex_name", "measure_name"])[col].cumsum()
    return out


def compare_chain_to_endpoint(cumulative_results: pd.DataFrame, endpoint_results: pd.DataFrame) -> pd.DataFrame:
    final = (
        cumulative_results.sort_values("end_year")
        .groupby(["location_name", "sex_name", "measure_name"], as_index=False)
        .tail(1)
        .rename(
            columns={
                "cumulative_population_growth": "chain_population_growth",
                "cumulative_population_aging": "chain_population_aging",
                "cumulative_rate_change": "chain_rate_change",
                "cumulative_total_change": "chain_total_change",
            }
        )
    )
    final = final[
        [
            "location_name",
            "sex_name",
            "measure_name",
            "chain_population_growth",
            "chain_population_aging",
            "chain_rate_change",
            "chain_total_change",
        ]
    ]
    out = endpoint_results.merge(final, on=["location_name", "sex_name", "measure_name"], how="inner")
    out["population_growth_difference"] = out["chain_population_growth"] - out["population_growth"]
    out["population_aging_difference"] = out["chain_population_aging"] - out["population_aging"]
    out["rate_change_difference"] = out["chain_rate_change"] - out["rate_change"]
    out["total_change_difference"] = out["chain_total_change"] - out["total_change"]
    return out.sort_values(["sex_name", "measure_name", "location_name"]).reset_index(drop=True)


def build_chained_decomposition_tables(df: pd.DataFrame, endpoint_results: pd.DataFrame, sexes=SEXES):
    pop = infer_population_table(df)
    row_driver = endpoint_results[["location_name", "sex_name", "measure_name"]].drop_duplicates()
    row_driver = row_driver[row_driver["sex_name"].isin(sexes) & row_driver["measure_name"].isin(MAIN)]
    annual = run_decomposition_interval_set(df, pop, row_driver, range(1990, 2024))
    interval = run_decomposition_interval_set(df, pop, row_driver, DECOMP_INTERVAL_YEARS)
    annual_cumulative = make_cumulative_decomposition(annual)
    interval_cumulative = make_cumulative_decomposition(interval)
    annual_vs_endpoint = compare_chain_to_endpoint(annual_cumulative, endpoint_results)
    interval_vs_endpoint = compare_chain_to_endpoint(interval_cumulative, endpoint_results)
    final_summary = annual_vs_endpoint[
        [
            "location_name",
            "sex_name",
            "measure_name",
            "population_growth",
            "chain_population_growth",
            "population_aging",
            "chain_population_aging",
            "rate_change",
            "chain_rate_change",
            "total_change",
            "chain_total_change",
            "total_change_difference",
        ]
    ].rename(
        columns={
            "population_growth": "endpoint_population_growth",
            "chain_population_growth": "annual_chain_population_growth",
            "population_aging": "endpoint_population_aging",
            "chain_population_aging": "annual_chain_population_aging",
            "rate_change": "endpoint_rate_change",
            "chain_rate_change": "annual_chain_rate_change",
            "total_change": "endpoint_total_change",
            "chain_total_change": "annual_chain_total_change",
        }
    )
    return annual, annual_cumulative, interval, interval_cumulative, annual_vs_endpoint, interval_vs_endpoint, final_summary


def build_age_specific_comparison(df: pd.DataFrame) -> pd.DataFrame:
    out = subset(
        df,
        measure_name=TREND_PATTERN_MEASURES,
        metric_name=("Rate", "Number", "Percent"),
        age_name=ADULT,
    )
    out = out[out["year"].between(1990, 2023)].copy()
    out["age_order"] = out["age_name"].map(age_order_index)
    return out.sort_values(["metric_name", "sex_name", "measure_name", "location_name", "year", "age_order"]).drop(columns=["age_order"]).reset_index(drop=True)


def build_age_difference_table(age_specific: pd.DataFrame) -> pd.DataFrame:
    piv = age_specific.pivot_table(
        index=["sex_name", "measure_name", "metric_name", "age_name", "year"],
        columns="location_name",
        values="val",
        aggfunc="first",
    ).reset_index()
    piv["china_us_pct_difference"] = np.where(
        piv["United States of America"] > 0,
        100.0 * (piv["China"] / piv["United States of America"] - 1.0),
        np.nan,
    )
    piv["china_us_absolute_difference"] = piv["China"] - piv["United States of America"]
    piv["age_order"] = piv["age_name"].map(age_order_index)
    return piv.sort_values(["metric_name", "sex_name", "measure_name", "year", "age_order"]).drop(columns=["age_order"]).reset_index(drop=True)


def recent_aligned_apc_periods(years: pd.Series) -> pd.DataFrame:
    years = pd.Series(years, copy=False).astype(int)
    block_index = ((APC_YEAR_END - years) // APC_PERIOD_WIDTH).astype(int)
    period_end = APC_YEAR_END - block_index * APC_PERIOD_WIDTH
    period_start = np.maximum(APC_YEAR_START, period_end - APC_PERIOD_WIDTH + 1)
    return pd.DataFrame(
        {
            "period_start": period_start.astype(int),
            "period_end": period_end.astype(int),
            "period_label": period_start.astype(str) + "-" + period_end.astype(str),
            "period_midpoint": (period_start + period_end) / 2.0,
            "period_years": (period_end - period_start + 1).astype(int),
        },
        index=years.index,
    )


def normw(w):
    w = np.asarray(w, float)
    pos = w[np.isfinite(w) & (w > 0)]
    scale = float(np.nanmedian(pos)) if pos.size else 1.0
    scale = 1.0 if not np.isfinite(scale) or scale <= 0 else scale
    w = w / scale
    w[~np.isfinite(w) | (w <= 0)] = 1.0
    return w


def apc_from_slope(slope, se=None):
    if se is None or np.isnan(se):
        return 100.0 * (math.exp(slope) - 1.0), np.nan, np.nan
    lower = 100.0 * (math.exp(slope - 1.96 * se) - 1.0)
    upper = 100.0 * (math.exp(slope + 1.96 * se) - 1.0)
    return 100.0 * (math.exp(slope) - 1.0), lower, upper


def wlm(X, y, w, a=0.0, m=None):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    w = normw(w)
    keep = np.isfinite(y) & np.all(np.isfinite(X), axis=1) & np.isfinite(w)
    X, y, w = X[keep], y[keep], w[keep]
    root_w = np.sqrt(w)
    Xw = X * root_w[:, None]
    yw = y * root_w
    if m is None:
        m = np.ones(X.shape[1], float)
        m[0] = 0.0
    penalty = a * np.diag(np.asarray(m, float))
    xtx = Xw.T @ Xw + penalty
    beta = np.linalg.pinv(xtx) @ (Xw.T @ yw)
    resid = y - X @ beta
    dof = max(X.shape[0] - X.shape[1], 1)
    sigma2 = float(np.sum(w * resid**2) / dof)
    cov = sigma2 * np.linalg.pinv(xtx)
    return {"beta": beta, "cov": cov, "resid": resid, "n": int(X.shape[0])}


def contrast(beta, cov, vector):
    vector = np.asarray(vector, float)
    estimate = float(vector @ beta)
    variance = float(vector @ cov @ vector)
    se = math.sqrt(max(variance, 0.0))
    return estimate, se, estimate - 1.96 * se, estimate + 1.96 * se


def indicator_matrix(values, levels):
    return np.column_stack([(values == level).astype(float) for level in levels])


def apc_contrast_vector(beta, start, index, reference):
    vector = np.zeros_like(beta)
    vector[start + index] = 1.0
    vector[start + reference] = -1.0
    return vector


def build_apc_inputs(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = subset(df, measure_name=APC_MEASURES, age_name=APCAGES, location_name=LOCS, sex_name=SEXES)
    base = base[base["year"].between(APC_YEAR_START, APC_YEAR_END)].copy()
    keys = ["location_name", "sex_name", "measure_name", "age_name", "year"]
    numbers = renamed_values(subset(base, metric_name="Number"), "count")
    rates = renamed_values(subset(base, metric_name="Rate"), "rate")
    percents = renamed_values(subset(base, metric_name="Percent"), "percent")
    annual = (
        numbers.merge(rates[[*keys, "rate_value", "rate_lower", "rate_upper"]], on=keys, how="inner", validate="one_to_one")
        .merge(percents[[*keys, "percent_value", "percent_lower", "percent_upper"]], on=keys, how="left", validate="one_to_one")
    )
    annual["population_estimate"] = np.where(annual["rate_value"] > 0, annual["count_value"] / annual["rate_value"] * 100000.0, np.nan)
    annual["burden_total_estimate"] = np.where(annual["percent_value"] > 0, annual["count_value"] / annual["percent_value"], np.nan)
    annual["age_midpoint"] = annual["age_name"].map(MID)
    annual = annual.join(recent_aligned_apc_periods(annual["year"]))
    annual["cohort_midpoint"] = annual["period_midpoint"] - annual["age_midpoint"]
    annual["cohort_label"] = annual["cohort_midpoint"].round().astype(int).astype(str)
    grouping = [
        "location_name",
        "sex_name",
        "measure_name",
        "age_name",
        "age_midpoint",
        "period_start",
        "period_end",
        "period_label",
        "period_midpoint",
        "period_years",
        "cohort_label",
        "cohort_midpoint",
    ]
    cells = annual.groupby(grouping, as_index=False).agg(
        count_value=("count_value", "sum"),
        population_sum=("population_estimate", "sum"),
        burden_total_sum=("burden_total_estimate", "sum"),
    )
    cells["rate_value"] = np.where(cells["population_sum"] > 0, cells["count_value"] / cells["population_sum"] * 100000.0, np.nan)
    cells["percent_value"] = np.where(cells["burden_total_sum"] > 0, cells["count_value"] / cells["burden_total_sum"], np.nan)
    cells["age_order"] = cells["age_name"].map(age_order_index)
    cells = cells.sort_values(["sex_name", "measure_name", "location_name", "period_start", "age_order"]).drop(columns=["age_order"]).reset_index(drop=True)
    return annual, cells


def build_apc_outputs(df: pd.DataFrame):
    annual, cells = build_apc_inputs(df)
    summary, age_curve, drift, period_rr, cohort_rr = [], [], [], [], []
    for (location, sex, measure), cell in cells.groupby(["location_name", "sex_name", "measure_name"]):
        cell = cell.sort_values(["period_start", "age_midpoint"]).copy()
        ann = subset(annual, location_name=location, sex_name=sex, measure_name=measure).sort_values(["year", "age_midpoint"])
        ages = [age for age in APCAGES if age in set(cell["age_name"])]
        pmeta = cell[["period_label", "period_start", "period_midpoint"]].drop_duplicates().sort_values("period_start")
        cmeta = cell[["cohort_label", "cohort_midpoint"]].drop_duplicates().sort_values("cohort_midpoint")
        periods = pmeta["period_label"].tolist()
        cohorts = cmeta["cohort_label"].tolist()
        pmap = dict(zip(pmeta["period_label"], pmeta["period_midpoint"]))
        cmap = dict(zip(cmeta["cohort_label"], cmeta["cohort_midpoint"]))

        X = np.column_stack(
            [
                np.ones(len(cell)),
                indicator_matrix(cell["age_name"], ages),
                indicator_matrix(cell["period_label"], periods),
                indicator_matrix(cell["cohort_label"], cohorts),
            ]
        )
        mask = np.ones(X.shape[1], float)
        mask[0] = 0.0
        fit = wlm(
            X,
            np.log(cell["rate_value"].clip(lower=1e-12).to_numpy(float)),
            cell["population_sum"].to_numpy(float),
            5.0,
            mask,
        )

        age_start = 1
        per_start = age_start + len(ages)
        coh_start = per_start + len(periods)
        ref_age = APCREF if APCREF in ages else ages[len(ages) // 2]
        ref_per = periods[len(periods) // 2]
        ref_coh = cohorts[len(cohorts) // 2]
        ia, ip, ic = ages.index(ref_age), periods.index(ref_per), cohorts.index(ref_coh)

        net_x = np.column_stack(
            [
                np.ones(len(ann)),
                ann["year"].to_numpy(float) - ann["year"].mean(),
                indicator_matrix(ann["age_name"], ages),
            ]
        )
        net = wlm(
            net_x,
            np.log(ann["rate_value"].clip(lower=1e-12).to_numpy(float)),
            ann["population_estimate"].to_numpy(float),
        )
        nd, ndl, ndu = apc_from_slope(float(net["beta"][1]), math.sqrt(max(float(net["cov"][1, 1]), 0.0)))
        summary.append(
            {
                "location_name": location,
                "sex_name": sex,
                "measure_name": measure,
                "model_year_start": APC_YEAR_START,
                "model_year_end": APC_YEAR_END,
                "n_age_groups": len(ages),
                "n_period_groups": len(periods),
                "n_cohort_groups": len(cohorts),
                "reference_age": ref_age,
                "reference_period": ref_per,
                "reference_cohort": ref_coh,
                "ridge_alpha": 5.0,
                "net_drift": nd,
                "net_drift_lower": ndl,
                "net_drift_upper": ndu,
            }
        )

        for j, age in enumerate(ages):
            lp_c = np.zeros_like(fit["beta"])
            lp_c[0] = lp_c[age_start + j] = lp_c[per_start + ip] = lp_c[coh_start + ic] = 1.0
            lp, _, ll, lu = contrast(fit["beta"], fit["cov"], lp_c)
            rr, se, rl, ru = contrast(fit["beta"], fit["cov"], apc_contrast_vector(fit["beta"], age_start, j, ia))
            age_curve.append(
                {
                    "location_name": location,
                    "sex_name": sex,
                    "measure_name": measure,
                    "age_name": age,
                    "age_midpoint": MID[age],
                    "reference_age": ref_age,
                    "reference_period": ref_per,
                    "reference_cohort": ref_coh,
                    "fitted_rate": math.exp(lp),
                    "fitted_rate_lower": math.exp(ll),
                    "fitted_rate_upper": math.exp(lu),
                    "age_rr_vs_reference": math.exp(rr),
                    "age_rr_vs_reference_lower": math.exp(rl),
                    "age_rr_vs_reference_upper": math.exp(ru),
                    "age_rr_vs_reference_se": se,
                }
            )
            age_ann = ann[ann["age_name"] == age].sort_values("year")
            local_x = np.column_stack([np.ones(len(age_ann)), age_ann["year"].to_numpy(float) - age_ann["year"].mean()])
            local = wlm(
                local_x,
                np.log(age_ann["rate_value"].clip(lower=1e-12).to_numpy(float)),
                age_ann["population_estimate"].to_numpy(float),
            )
            ld, ll0, lu0 = apc_from_slope(float(local["beta"][1]), math.sqrt(max(float(local["cov"][1, 1]), 0.0)))
            drift.append(
                {
                    "location_name": location,
                    "sex_name": sex,
                    "measure_name": measure,
                    "age_name": age,
                    "age_midpoint": MID[age],
                    "local_drift": ld,
                    "local_drift_lower": ll0,
                    "local_drift_upper": lu0,
                }
            )

        for j, period in enumerate(periods):
            rr, _, rl, ru = contrast(fit["beta"], fit["cov"], apc_contrast_vector(fit["beta"], per_start, j, ip))
            period_rr.append(
                {
                    "location_name": location,
                    "sex_name": sex,
                    "measure_name": measure,
                    "period_label": period,
                    "period_midpoint": pmap[period],
                    "reference_period": ref_per,
                    "period_rr": math.exp(rr),
                    "period_rr_lower": math.exp(rl),
                    "period_rr_upper": math.exp(ru),
                }
            )

        for j, cohort in enumerate(cohorts):
            rr, _, rl, ru = contrast(fit["beta"], fit["cov"], apc_contrast_vector(fit["beta"], coh_start, j, ic))
            cohort_rr.append(
                {
                    "location_name": location,
                    "sex_name": sex,
                    "measure_name": measure,
                    "cohort_label": cohort,
                    "cohort_midpoint": cmap[cohort],
                    "reference_cohort": ref_coh,
                    "cohort_rr": math.exp(rr),
                    "cohort_rr_lower": math.exp(rl),
                    "cohort_rr_upper": math.exp(ru),
                }
            )

    return (
        pd.DataFrame(summary).sort_values(["sex_name", "measure_name", "location_name"]).reset_index(drop=True),
        pd.DataFrame(age_curve).sort_values(["sex_name", "measure_name", "location_name", "age_midpoint"]).reset_index(drop=True),
        pd.DataFrame(drift).sort_values(["sex_name", "measure_name", "location_name", "age_midpoint"]).reset_index(drop=True),
        pd.DataFrame(period_rr).sort_values(["sex_name", "measure_name", "location_name", "period_midpoint"]).reset_index(drop=True),
        pd.DataFrame(cohort_rr).sort_values(["sex_name", "measure_name", "location_name", "cohort_midpoint"]).reset_index(drop=True),
        cells.sort_values(["sex_name", "measure_name", "location_name", "period_start", "age_midpoint"]).reset_index(drop=True),
    )


def line_label(location: str, sex: str) -> str:
    return f"{SLOC[location]} {sex}"


def plot_trend_panel(df: pd.DataFrame, value_label: str, title: str, key: str, measures=TMEAS) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    axes = axes.ravel()
    for ax, measure in zip(axes, measures):
        panel = subset(df, measure_name=measure)
        for (location, sex), s in panel.groupby(["location_name", "sex_name"]):
            s = s.sort_values("year")
            ax.plot(s["year"], s["val"], color=COL[location], linestyle=SEX_STYLE[sex], linewidth=2, label=line_label(location, sex))
        ax.set_title(measure)
        ax.set_ylabel(value_label)
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)
    fig.suptitle(title, y=0.98, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.07, 1, 0.95])
    return save_fig(fig, key)


def plot_segmented_trend_panel(obs: pd.DataFrame, fitted: pd.DataFrame, summary: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    axes = axes.ravel()
    for ax, measure in zip(axes, TREND_PATTERN_MEASURES):
        panel = subset(obs, measure_name=measure)
        fit_panel = subset(fitted, measure_name=measure)
        for (location, sex), s in panel.groupby(["location_name", "sex_name"]):
            s = s.sort_values("year")
            ax.scatter(s["year"], s["val"], color=COL[location], marker=SEX_MARKER[sex], alpha=0.35, s=14)
        for (location, sex), s in fit_panel.groupby(["location_name", "sex_name"]):
            s = s.sort_values("year")
            ax.plot(s["year"], s["fitted"], color=COL[location], linestyle=SEX_STYLE[sex], linewidth=2.1, label=line_label(location, sex))
        ax.set_title(measure)
        ax.set_ylabel("Age-standardized rate per 100,000")
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)
    fig.suptitle("Segmented log-linear trends in schizophrenia age-standardized rates", y=0.98, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.07, 1, 0.95])
    return save_fig(fig, "segmented")


def plot_age_specific_endpoint(age: pd.DataFrame, metric: str, key: str, value_label: str, title: str) -> Path:
    plot_df = subset(age, metric_name=metric)
    if metric == "Percent":
        plot_df = plot_df.copy()
        plot_df[VALUE_COLS] = plot_df[VALUE_COLS] * 100.0
    plot_df = plot_df[plot_df["year"].isin(ENDPOINT_YEARS)].copy()
    fig, axes = plt.subplots(len(TMEAS), len(SEXES), figsize=(15, 12), sharex=True)
    age_positions = np.arange(len(ADULT))
    for r, measure in enumerate(TMEAS):
        for c, sex in enumerate(SEXES):
            ax = axes[r, c]
            panel = subset(plot_df, measure_name=measure, sex_name=sex)
            for (location, year), s in panel.groupby(["location_name", "year"]):
                s = s.set_index("age_name").reindex(ADULT).reset_index()
                ax.plot(
                    age_positions,
                    s["val"],
                    color=COL[location],
                    linestyle="-" if year == 2023 else ":",
                    marker="o" if year == 2023 else "s",
                    linewidth=1.8,
                    markersize=3,
                    label=f"{SLOC[location]} {year}",
                )
            ax.set_title(f"{measure} - {sex}")
            ax.set_ylabel(value_label)
            ax.grid(alpha=0.25)
            if r == len(TMEAS) - 1:
                ax.set_xticks(age_positions)
                ax.set_xticklabels(ADULT, rotation=45, ha="right")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)
    fig.suptitle(title, y=0.995, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.06, 1, 0.97])
    return save_fig(fig, key)


def symmetric_limit(values: pd.Series) -> float:
    finite = values[np.isfinite(values)]
    if finite.empty:
        return 1.0
    lim = np.nanpercentile(np.abs(finite), 95)
    return max(float(lim), 1.0)


def plot_age_difference_heatmap(diff: pd.DataFrame, metric: str, key: str, title: str) -> Path:
    panel_df = subset(diff, metric_name=metric)
    limit = symmetric_limit(panel_df["china_us_pct_difference"])
    fig, axes = plt.subplots(len(TMEAS), len(SEXES), figsize=(14, 12), sharex=True, sharey=True, constrained_layout=True)
    im = None
    for r, measure in enumerate(TMEAS):
        for c, sex in enumerate(SEXES):
            ax = axes[r, c]
            s = subset(panel_df, measure_name=measure, sex_name=sex)
            mat = s.pivot_table(index="age_name", columns="year", values="china_us_pct_difference", aggfunc="first").reindex(ADULT)
            im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
            ax.set_title(f"{measure} - {sex}")
            if c == 0:
                ax.set_yticks(np.arange(len(ADULT)))
                ax.set_yticklabels(ADULT)
            if r == len(TMEAS) - 1:
                years = list(mat.columns)
                tick_idx = [i for i, y in enumerate(years) if y in (1990, 2000, 2010, 2020, 2023)]
                ax.set_xticks(tick_idx)
                ax.set_xticklabels([years[i] for i in tick_idx], rotation=45, ha="right")
    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
        cbar.set_label("China vs United States difference (%)")
    fig.suptitle(title, fontsize=14, fontweight="semibold")
    return save_fig(fig, key)


def plot_age_metric_heatmap(age: pd.DataFrame, metric: str, key: str, title: str, value_label: str) -> Path:
    panel_df = subset(age, metric_name=metric).copy()
    if metric == "Percent":
        panel_df[VALUE_COLS] = panel_df[VALUE_COLS] * 100.0
    combos = [(location, sex) for sex in SEXES for location in LOCS]
    fig, axes = plt.subplots(len(TMEAS), len(combos), figsize=(16, 12), sharex=True, sharey=True, constrained_layout=True)
    for r, measure in enumerate(TMEAS):
        measure_values = subset(panel_df, measure_name=measure)["val"]
        vmax = np.nanpercentile(measure_values, 98) if measure_values.notna().any() else 1.0
        vmax = vmax if vmax > 0 else 1.0
        for c, (location, sex) in enumerate(combos):
            ax = axes[r, c]
            s = subset(panel_df, measure_name=measure, location_name=location, sex_name=sex)
            mat = s.pivot_table(index="age_name", columns="year", values="val", aggfunc="first").reindex(ADULT)
            im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)
            ax.set_title(f"{measure}\n{SLOC[location]} {sex}")
            if c == 0:
                ax.set_yticks(np.arange(len(ADULT)))
                ax.set_yticklabels(ADULT)
            if r == len(TMEAS) - 1:
                years = list(mat.columns)
                tick_idx = [i for i, y in enumerate(years) if y in (1990, 2000, 2010, 2020, 2023)]
                ax.set_xticks(tick_idx)
                ax.set_xticklabels([years[i] for i in tick_idx], rotation=45, ha="right")
        fig.colorbar(im, ax=axes[r, :].ravel().tolist(), shrink=0.65, pad=0.01, label=f"{measure} {value_label}")
    fig.suptitle(title, fontsize=14, fontweight="semibold")
    return save_fig(fig, key)


def plot_decomposition(df: pd.DataFrame) -> Path:
    components = ["population_growth", "population_aging", "rate_change"]
    labels = ["Population growth", "Population aging", "Rate change"]
    colors = ["#4c78a8", "#f58518", "#54a24b"]
    fig, axes = plt.subplots(len(MAIN), len(SEXES), figsize=(13, 9), sharex=True)
    for r, measure in enumerate(MAIN):
        for c, sex in enumerate(SEXES):
            ax = axes[r, c]
            panel = subset(df, measure_name=measure, sex_name=sex).set_index("location_name").reindex(LOCS).reset_index()
            x = np.arange(len(LOCS))
            pos_bottom = np.zeros(len(panel))
            neg_bottom = np.zeros(len(panel))
            for comp, label, color in zip(components, labels, colors):
                vals = panel[comp].to_numpy(dtype=float)
                bottoms = np.where(vals >= 0, pos_bottom, neg_bottom)
                ax.bar(x, vals, bottom=bottoms, color=color, width=0.6, label=label)
                pos_bottom += np.where(vals >= 0, vals, 0)
                neg_bottom += np.where(vals < 0, vals, 0)
            ax.axhline(0, color="black", linewidth=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels([SLOC[loc] for loc in LOCS])
            ax.set_title(f"{measure} - {sex}")
            ax.set_ylabel("Adult count change")
            ax.grid(axis="y", alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3)
    fig.suptitle("Decomposition of adult schizophrenia burden change, 1990-2023", y=0.98, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.07, 1, 0.95])
    return save_fig(fig, "decomposition")


def plot_annual_chain_cumulative(df: pd.DataFrame) -> Path:
    components = [
        ("cumulative_population_growth", "Population growth", "#4c78a8"),
        ("cumulative_population_aging", "Population aging", "#f58518"),
        ("cumulative_rate_change", "Rate change", "#54a24b"),
    ]
    rows = [(measure, sex) for measure in MAIN for sex in SEXES]
    fig, axes = plt.subplots(len(rows), len(LOCS), figsize=(14, 16), sharex=True)
    for r, (measure, sex) in enumerate(rows):
        for c, location in enumerate(LOCS):
            ax = axes[r, c]
            panel = subset(df, measure_name=measure, sex_name=sex, location_name=location).sort_values("end_year")
            for col, label, color in components:
                ax.plot(panel["end_year"], panel[col], color=color, linewidth=1.8, label=label)
            ax.plot(panel["end_year"], panel["cumulative_total_change"], color="#111827", linewidth=1.5, linestyle="--", label="Total change")
            ax.axhline(0, color="black", linewidth=0.7)
            ax.set_title(f"{measure} - {sex} - {SLOC[location]}")
            ax.set_ylabel("Cumulative change since 1990")
            ax.grid(alpha=0.25)
            if r == len(rows) - 1:
                ax.set_xlabel("End year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)
    fig.suptitle("Annual chained decomposition of adult schizophrenia burden change", y=0.99, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.045, 1, 0.97])
    return save_fig(fig, "annual_chain")


def plot_interval_chain_components(df: pd.DataFrame) -> Path:
    components = [
        ("population_growth", "Population growth", "#4c78a8"),
        ("population_aging", "Population aging", "#f58518"),
        ("rate_change", "Rate change", "#54a24b"),
    ]
    rows = [(measure, sex) for measure in MAIN for sex in SEXES]
    fig, axes = plt.subplots(len(rows), len(LOCS), figsize=(15, 16), sharex=True)
    for r, (measure, sex) in enumerate(rows):
        for c, location in enumerate(LOCS):
            ax = axes[r, c]
            panel = subset(df, measure_name=measure, sex_name=sex, location_name=location).sort_values("start_year")
            x = np.arange(len(panel))
            pos_bottom = np.zeros(len(panel))
            neg_bottom = np.zeros(len(panel))
            for comp, label, color in components:
                vals = panel[comp].to_numpy(dtype=float)
                bottoms = np.where(vals >= 0, pos_bottom, neg_bottom)
                ax.bar(x, vals, bottom=bottoms, color=color, width=0.68, label=label)
                pos_bottom += np.where(vals >= 0, vals, 0)
                neg_bottom += np.where(vals < 0, vals, 0)
            ax.scatter(x, panel["total_change"], color="#111827", marker="D", s=18, label="Total change", zorder=4)
            ax.axhline(0, color="black", linewidth=0.7)
            ax.set_title(f"{measure} - {sex} - {SLOC[location]}")
            ax.set_ylabel("Interval change")
            ax.set_xticks(x)
            ax.set_xticklabels(panel["interval_label"], rotation=45, ha="right")
            ax.grid(axis="y", alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)
    fig.suptitle("Interval chained decomposition of adult schizophrenia burden change", y=0.99, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=[0, 0.045, 1, 0.97])
    return save_fig(fig, "interval_chain")


def plot_event_impact_panel(summary: pd.DataFrame) -> Path:
    if summary.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No event-impact rows available", ha="center", va="center")
        ax.axis("off")
        return save_fig(fig, "event_impact")

    plot_df = summary.copy()
    plot_df["event_label"] = (
        plot_df["event_name"]
        + " ("
        + plot_df["event_year"].astype(str)
        + ") | "
        + plot_df["location_name"].map(SLOC)
        + " "
        + plot_df["sex_name"]
    )
    row_order = plot_df[["event_year", "event_label"]].drop_duplicates().sort_values(["event_year", "event_label"])["event_label"].tolist()
    mat = plot_df.pivot_table(index="event_label", columns="measure_name", values="max_excess_pct", aggfunc="first").reindex(row_order)[list(TMEAS)]
    limit = symmetric_limit(mat.stack())
    fig_height = max(6.0, 0.42 * len(mat) + 2.0)
    fig, ax = plt.subplots(figsize=(11, fig_height), constrained_layout=True)
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index)
    ax.set_title("ASR deviation during event windows")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            value = mat.iloc[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.1f}%", ha="center", va="center", fontsize=8, color="#111827")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Max deviation vs pre-event baseline (%)")
    return save_fig(fig, "event_impact")


def plot_apc_panel(cells: pd.DataFrame) -> Path:
    plot_df = subset(cells, measure_name=APC_PLOT_MEASURES)
    rows = [(measure, sex) for measure in APC_PLOT_MEASURES for sex in SEXES]
    fig, axes = plt.subplots(len(rows), len(LOCS), figsize=(12, 11), sharex=True, sharey=True, constrained_layout=True)
    im = None
    for r, (measure, sex) in enumerate(rows):
        measure_values = subset(plot_df, measure_name=measure)["rate_value"]
        vmax = np.nanpercentile(measure_values, 98) if measure_values.notna().any() else 1.0
        vmax = vmax if vmax > 0 else 1.0
        for c, location in enumerate(LOCS):
            ax = axes[r, c]
            s = subset(plot_df, location_name=location, sex_name=sex, measure_name=measure)
            mat = s.pivot_table(index="age_name", columns="period_label", values="rate_value", aggfunc="first").reindex(APCAGES)
            im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="YlOrRd", vmin=0, vmax=vmax)
            ax.set_title(f"{measure} - {SLOC[location]} {sex}")
            if c == 0:
                ax.set_yticks(np.arange(len(APCAGES)))
                ax.set_yticklabels(APCAGES)
            if r == len(rows) - 1:
                ax.set_xticks(np.arange(len(mat.columns)))
                ax.set_xticklabels(mat.columns, rotation=45, ha="right")
    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
        cbar.set_label("Rate per 100,000")
    fig.suptitle("Age-period schizophrenia rate surfaces", fontsize=14, fontweight="semibold")
    return save_fig(fig, "apc")


def plot_apc_marginal_panel(age: pd.DataFrame, drift: pd.DataFrame, period: pd.DataFrame, cohort: pd.DataFrame) -> Path:
    row_specs = [(measure, sex) for measure in APC_MEASURES for sex in SEXES]
    fig, axes = plt.subplots(nrows=len(row_specs), ncols=4, figsize=(22, 22), squeeze=False)
    for i, (measure, sex) in enumerate(row_specs):
        ax1, ax2, ax3, ax4 = axes[i]
        for location in LOCS:
            color = COL[location]
            label = SLOC[location]
            a = subset(age, measure_name=measure, sex_name=sex, location_name=location).sort_values("age_midpoint")
            d = subset(drift, measure_name=measure, sex_name=sex, location_name=location).sort_values("age_midpoint")
            p = subset(period, measure_name=measure, sex_name=sex, location_name=location).sort_values("period_midpoint")
            c = subset(cohort, measure_name=measure, sex_name=sex, location_name=location).sort_values("cohort_midpoint")

            if not a.empty:
                x1 = np.arange(len(a))
                ax1.plot(x1, a["fitted_rate"], color=color, linewidth=1.8, marker="o", markersize=3, label=label)
                ax1.fill_between(x1, a["fitted_rate_lower"], a["fitted_rate_upper"], color=color, alpha=0.14)
            if not d.empty:
                x2 = np.arange(len(d))
                ax2.plot(x2, d["local_drift"], color=color, linewidth=1.8, marker="o", markersize=3, label=label)
                ax2.fill_between(x2, d["local_drift_lower"], d["local_drift_upper"], color=color, alpha=0.14)
            if not p.empty:
                ax3.plot(p["period_midpoint"], p["period_rr"], color=color, linewidth=1.8, marker="o", markersize=3, label=label)
                ax3.fill_between(p["period_midpoint"], p["period_rr_lower"], p["period_rr_upper"], color=color, alpha=0.14)
            if not c.empty:
                ax4.plot(c["cohort_midpoint"], c["cohort_rr"], color=color, linewidth=1.8, marker="o", markersize=3, label=label)
                ax4.fill_between(c["cohort_midpoint"], c["cohort_rr_lower"], c["cohort_rr_upper"], color=color, alpha=0.14)

        ax1.set_title(f"{measure} - {sex} | Age curve")
        ax1.set_ylabel("Fitted rate per 100,000")
        if not subset(age, measure_name=measure, sex_name=sex).empty:
            ref = subset(age, measure_name=measure, sex_name=sex).sort_values("age_midpoint").drop_duplicates("age_name")
            ax1.set_xticks(np.arange(len(ref)))
            ax1.set_xticklabels(ref["age_name"], rotation=65, ha="right")

        ax2.set_title(f"{measure} - {sex} | Local drift")
        ax2.set_ylabel("Annual % change")
        ax2.axhline(0.0, color="#666666", linestyle="--", linewidth=1)
        if not subset(drift, measure_name=measure, sex_name=sex).empty:
            ref = subset(drift, measure_name=measure, sex_name=sex).sort_values("age_midpoint").drop_duplicates("age_name")
            ax2.set_xticks(np.arange(len(ref)))
            ax2.set_xticklabels(ref["age_name"], rotation=65, ha="right")

        ax3.set_title(f"{measure} - {sex} | Period RR")
        ax3.set_ylabel("RR vs reference period")
        ax3.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
        ax3.set_xlabel("Period midpoint")

        ax4.set_title(f"{measure} - {sex} | Cohort RR")
        ax4.set_ylabel("RR vs reference cohort")
        ax4.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
        ax4.set_xlabel("Birth cohort midpoint")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2)
    fig.suptitle("Schizophrenia ridge-regularized APC marginal effects", y=0.995, fontsize=15, fontweight="semibold")
    fig.text(
        0.5,
        0.012,
        "Annual observations are collapsed into recent-aligned periods; cohort contrasts are shown as relative-risk patterns.",
        ha="center",
        va="bottom",
        fontsize=9,
        wrap=True,
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.982))
    return save_fig(fig, "apc_marginal")


def plot_apc_cohort_panel(cohort: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(nrows=len(APC_MEASURES), ncols=len(SEXES), figsize=(15, 10), sharex=True, squeeze=False)
    for r, measure in enumerate(APC_MEASURES):
        for cidx, sex in enumerate(SEXES):
            ax = axes[r, cidx]
            for location in LOCS:
                s = subset(cohort, measure_name=measure, sex_name=sex, location_name=location).sort_values("cohort_midpoint")
                if s.empty:
                    continue
                ax.plot(s["cohort_midpoint"], s["cohort_rr"], color=COL[location], linewidth=2.0, marker="o", markersize=3, label=SLOC[location])
                ax.fill_between(s["cohort_midpoint"], s["cohort_rr_lower"], s["cohort_rr_upper"], color=COL[location], alpha=0.14)
            ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
            ax.set_title(f"{measure} - {sex}")
            ax.set_ylabel("Cohort RR")
            ax.grid(alpha=0.25)
            if r == len(APC_MEASURES) - 1:
                ax.set_xlabel("Birth cohort midpoint")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2)
    fig.suptitle("Schizophrenia APC cohort relative-risk patterns", y=0.99, fontsize=14, fontweight="semibold")
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 0.965))
    return save_fig(fig, "apc_cohort")


def build_plot_explanations(results, config=None) -> dict[str, str]:
    return {
        "counts": "All-age numbers by location, sex, measure, and year.",
        "percent": "Preferred GBD percent-metric series, using age-standardized percent rows when present and all-age percent rows otherwise.",
        "asr": "Age-standardized rates per 100,000 by location, sex, measure, and year.",
        "segmented": "Observed and fitted segmented log-linear age-standardized rate trends selected by BIC.",
        "decomposition": "Das Gupta style average replacement decomposition of adult count changes into population growth, aging, and rate-change components.",
        "annual_chain": "Adjacent-year decomposition components accumulated from 1990; component allocations can differ from endpoint decomposition because the factor replacement is path-dependent.",
        "interval_chain": "Coarse chained interval decomposition, using the same factors as the endpoint decomposition.",
        "apc": "Ridge-regularized APC contrasts for recent-aligned age-period cells.",
        "apc_marginal": "APC marginal plots for fitted age curves, local drift, period RR, and cohort RR.",
        "apc_cohort": "Cohort relative-risk patterns from the ridge-regularized APC contrasts.",
    }


def build_all_outputs(bd: pd.DataFrame, ba: pd.DataFrame, config=None):
    seg, segd, segf = run_segmented_trend_analysis(bd)
    bs = build_burden_summary_table(bd)
    cm = build_comparative_metrics_table(bs)
    de = build_decomposition_table(bd, SEXES)
    dei = build_interval_decomposition_table(bd, SEXES)
    annual_chain, annual_cumulative, interval_chain, interval_cumulative, annual_vs_endpoint, interval_vs_endpoint, annual_final = build_chained_decomposition_tables(bd, de, SEXES)
    age = build_age_specific_comparison(bd)
    age_diff = build_age_difference_table(age)
    apc_summary, apc_age, apc_local, apc_period_rr, apc_cohort_rr, apc_cells = build_apc_outputs(bd)

    counts = subset(bd, metric_name="Number", age_name=ALL, measure_name=TMEAS).sort_values(["measure_name", "location_name", "sex_name", "year"]).copy()
    asr = subset(bd, metric_name="Rate", age_name=ASR, measure_name=TMEAS).sort_values(["measure_name", "location_name", "sex_name", "year"]).copy()
    percent = select_preferred_percent_rows(subset(bd, measure_name=TMEAS)).sort_values(["measure_name", "location_name", "sex_name", "year"]).copy()

    percent_plot = percent.copy()
    percent_plot[VALUE_COLS] = percent_plot[VALUE_COLS] * 100.0

    figures = {
        "counts": plot_trend_panel(counts, "Number", "Schizophrenia all-age counts in China and the United States (1990-2023)", "counts"),
        "percent": plot_trend_panel(percent_plot, "GBD percent metric (%)", "Schizophrenia measure-specific GBD percent metric (1990-2023)", "percent"),
        "asr": plot_trend_panel(asr, "Age-standardized rate per 100,000", "Schizophrenia age-standardized rates (1990-2023)", "asr"),
        "segmented": plot_segmented_trend_panel(asr, segf, seg),
        "age_rate": plot_age_specific_endpoint(age, "Rate", "age_rate", "Rate per 100,000", "Age-specific schizophrenia rates in 1990 and 2023"),
        "age_rate_difference": plot_age_difference_heatmap(age_diff, "Rate", "age_rate_difference", "China-vs-US age-specific schizophrenia rate differences"),
        "age_number": plot_age_specific_endpoint(age, "Number", "age_number", "Number", "Age-specific schizophrenia numbers in 1990 and 2023"),
        "age_number_difference": plot_age_difference_heatmap(age_diff, "Number", "age_number_difference", "China-vs-US age-specific schizophrenia number differences"),
        "age_percent": plot_age_specific_endpoint(age, "Percent", "age_percent", "GBD percent metric (%)", "Age-specific schizophrenia GBD percent metric in 1990 and 2023"),
        "age_percent_difference": plot_age_difference_heatmap(age_diff, "Percent", "age_percent_difference", "China-vs-US age-specific GBD percent-metric differences"),
        "age_number_heatmap": plot_age_metric_heatmap(age, "Number", "age_number_heatmap", "Age-specific schizophrenia number heatmaps, 1990-2023", "number"),
        "age_percent_heatmap": plot_age_metric_heatmap(age, "Percent", "age_percent_heatmap", "Age-specific schizophrenia GBD percent-metric heatmaps, 1990-2023", "GBD percent metric (%)"),
        "age_rate_heatmap": plot_age_metric_heatmap(age, "Rate", "age_rate_heatmap", "Age-specific schizophrenia rate heatmaps, 1990-2023", "rate per 100,000"),
        "decomposition": plot_decomposition(de),
        "annual_chain": plot_annual_chain_cumulative(annual_cumulative),
        "interval_chain": plot_interval_chain_components(interval_chain),
        "apc": plot_apc_panel(apc_cells),
        "apc_marginal": plot_apc_marginal_panel(apc_age, apc_local, apc_period_rr, apc_cohort_rr),
        "apc_cohort": plot_apc_cohort_panel(apc_cohort_rr),
    }

    tables = {
        "burden_audit": ba,
        "burden_summary": bs,
        "comparative_metrics": cm,
        "segmented_summary": seg,
        "segmented_segments": segd,
        "segmented_fitted": segf,
        "decomposition": de,
        "decomposition_interval": dei,
        "annual_chain_decomposition": annual_chain,
        "annual_chain_cumulative": annual_cumulative,
        "interval_chain_decomposition": interval_chain,
        "interval_chain_cumulative": interval_cumulative,
        "annual_chain_vs_endpoint": annual_vs_endpoint,
        "interval_chain_vs_endpoint": interval_vs_endpoint,
        "annual_chain_final_summary": annual_final,
        "age_specific_comparison": age,
        "age_difference": age_diff,
        "trend_counts": counts,
        "trend_asr": asr,
        "trend_percent": percent,
        "apc_summary": apc_summary,
        "apc_age_curve": apc_age,
        "apc_local_drift": apc_local,
        "apc_period_rr": apc_period_rr,
        "apc_cohort_rr": apc_cohort_rr,
        "apc_cells": apc_cells,
    }
    results = {"tables": tables, "figures": figures}
    results["explanations"] = build_plot_explanations(results, config)
    return results


def write_outputs(results) -> None:
    ensure_output_dirs()
    for path in STALE_OUTPUTS:
        if path.exists():
            path.unlink()
    for key, df in results["tables"].items():
        df.to_csv(DER / TABLE[key], index=False)
    EXP.write_text(json.dumps(results["explanations"], indent=2), encoding="utf-8")


def run_full_pipeline(burden_csv=None, config=None, save_outputs=True):
    configure_style()
    ensure_output_dirs()
    bd, bp = load_inputs(burden_csv)
    sexes = tuple(sex for sex in SEXES if sex in set(bd["sex_name"]))
    bd = bd[bd["sex_name"].isin(sexes)].copy()
    cfg = SimpleNamespace(sexes=sexes)
    ba = build_burden_audit_table(bd, bp, sexes)
    results = build_all_outputs(bd, ba, cfg)
    results["burden_source"] = bp
    results["config"] = cfg
    if save_outputs:
        write_outputs(results)
    return results


def print_run_summary(results) -> None:
    print(f"Burden source:        {results['burden_source']}")
    print(f"Active sexes:         {', '.join(results['config'].sexes)}")
    print("\nSaved figures:")
    for key, path in results["figures"].items():
        print(f"  {key}: {path}")
    print("\nSaved tables:")
    for key, filename in TABLE.items():
        print(f"  {key}: {DER / filename}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build China-vs-US schizophrenia GBD 1990-2023 analyses from the prepared cause export."
    )
    parser.add_argument("--burden-csv", type=Path, default=None)
    parser.add_argument("--no-save", action="store_true")
    args = parser.parse_args()
    print_run_summary(run_full_pipeline(args.burden_csv, None, not args.no_save))


if __name__ == "__main__":
    main()
