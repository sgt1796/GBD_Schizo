# Runs the NCI Age Period Cohort Web Tool's published apc2() implementation.
# Source: https://github.com/CBIIT/nci-webtools-dceg-age-period-cohort
# Pinned revision: 9e86d92d3f95b76b610a37f6e3c539ee126b5efd
# Only point estimates are exported: GBD modeled counts do not support the
# implementation's Poisson-count confidence intervals or Wald tests.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript nci_apc_runner.R input.csv output_dir source.R")
source(args[3])
cells <- read.csv(args[1], stringsAsFactors = FALSE)
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)

tables <- list(summary = list(), age_curve = list(), local_drift = list(),
               period_rr = list(), cohort_rr = list())
for (window_name in unique(cells$window)) {
  for (location in unique(cells$location_name)) {
    for (sex in unique(cells$sex_name)) {
      panel <- subset(cells, window == window_name & location_name == location & sex_name == sex)
      panel <- panel[order(panel$period_start, panel$age_lower), ]
      if (nrow(panel) != 72) stop("APC panel must have 12 ages by 6 periods")
      if (any(panel$events <= 0 | panel$population <= 0)) stop("APC cells must be positive")
      starts <- sort(unique(panel$period_start))
      if (length(starts) != 6 || any(diff(starts) != 5)) stop("Invalid five-year periods")
      if (!identical(sort(unique(panel$age_lower)), seq(10L, 65L, 5L))) stop("Invalid age groups")
      model_input <- list(name = paste(location, sex, window_name),
                          description = "GBD 2023 modeled incidence counts and population",
                          events = matrix(panel$events, nrow = 12),
                          offset = matrix(panel$population, nrow = 12),
                          a = seq(10, 70, 5), p = seq(starts[1], starts[1] + 30, 5),
                          ages = as.character(seq(10, 65, 5)),
                          periods = as.character(starts))
      fit <- apc2(model_input)
      common <- data.frame(window = window_name, location_name = location,
                           sex_name = sex, measure_name = "Incidence")
      add <- function(name, coordinate, values, value_name) {
        row <- common[rep(1, length(coordinate)), , drop = FALSE]
        row[[if (name %in% c("age_curve", "local_drift")) "age_midpoint" else
              if (name == "period_rr") "period_midpoint" else "cohort_midpoint"]] <- coordinate
        row[[value_name]] <- values
        tables[[name]][[length(tables[[name]]) + 1]] <<- row
      }
      row <- common
      row$net_drift <- unname(fit$NetDrift[1, 1])
      row$reference_age <- fit$Inputs$RVals[1]
      row$reference_period <- fit$Inputs$RVals[2]
      row$reference_cohort <- fit$Inputs$RVals[3]
      tables$summary[[length(tables$summary) + 1]] <- row
      add("age_curve", fit$LongAge[, 1], fit$LongAge[, 2], "longitudinal_age_rate_per_100000")
      add("local_drift", fit$LocalDrifts[, 1], fit$LocalDrifts[, 2], "local_drift")
      add("period_rr", fit$PeriodRR[, 1], fit$PeriodRR[, 2], "period_rr")
      add("cohort_rr", fit$CohortRR[, 1], fit$CohortRR[, 2], "cohort_rr")
    }
  }
}
for (name in names(tables)) {
  write.csv(do.call(rbind, tables[[name]]), file.path(args[2], paste0("apc_", name, ".csv")), row.names = FALSE)
}
