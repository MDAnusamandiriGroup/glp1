options(stringsAsFactors = FALSE)
Sys.setenv(PUBLICATION_PROFILE = "ejpc")

source("R/parameters.R")
source("R/model_rhs.R")
source("R/simulation.R")
source("R/calibration.R")
source("R/monte_carlo.R")
source("R/profile.R")
source("R/robustness.R")

required_packages <- c("deSolve", "ggplot2", "tidyr", "tidyselect")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    ". Run source('analysis/00_install_packages.R') first."
  )
}

psa_path <- "results/shared/probabilistic_sensitivity.rds"
if (!file.exists(psa_path)) {
  stop(
    "Existing PSA draws were not found at ", psa_path, ". ",
    "This finalization script intentionally does not rerun the PSA."
  )
}

dir.create("results/shared", recursive = TRUE, showWarnings = FALSE)
dir.create("results/ejpc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/ejpc", recursive = TRUE, showWarnings = FALSE)

# Rebuild deterministic EJPC outputs from the calibrated base case.
base_p <- calibrate_variant(list(
  effect_weight_exponent = 1,
  immediate_loss_after_discontinuation = FALSE
))
saveRDS(base_p, "results/shared/calibrated_parameters.rds")
write.csv(
  calibration_table(base_p),
  "results/ejpc/select_calibration_validation.csv",
  row.names = FALSE
)

core <- compare_publication_scenarios(base_p)
feedback <- compare_publication_scenarios(base_p, restart_feedback = TRUE)
saveRDS(core$simulations, "results/shared/core_simulations.rds")

core_summary <- profile_scenario_summary(core, "ejpc")
core_summary$analysis <- "constant restart primary"
feedback_summary <- profile_scenario_summary(feedback, "ejpc")
feedback_summary$analysis <- "weight-regain-dependent restart sensitivity"
write.csv(
  rbind(core_summary, feedback_summary),
  "results/ejpc/scenario_summary.csv",
  row.names = FALSE
)

curve <- persistence_response_curve(
  base_p,
  persistence_grid = seq(0.01, 1, by = 0.01),
  benchmark = "select"
)
thresholds <- refined_persistence_thresholds(base_p)
two_way <- two_way_sensitivity(base_p, benchmark = "select")
write.csv(
  curve, "results/ejpc/persistence_response_curve.csv", row.names = FALSE
)
write.csv(
  thresholds, "results/ejpc/persistence_thresholds.csv", row.names = FALSE
)
write.csv(
  two_way, "results/ejpc/two_way_sensitivity.csv", row.names = FALSE
)

# Structural uncertainty: recalibrate every variant to both SELECT trial arms.
structural_variants <- list(
  list(
    label = "gradual waning; eta=0.5",
    overrides = list(
      effect_weight_exponent = 0.5,
      immediate_loss_after_discontinuation = FALSE
    )
  ),
  list(
    label = "gradual waning; eta=1.0 (base case)",
    overrides = list(
      effect_weight_exponent = 1,
      immediate_loss_after_discontinuation = FALSE
    )
  ),
  list(
    label = "gradual waning; eta=1.5",
    overrides = list(
      effect_weight_exponent = 1.5,
      immediate_loss_after_discontinuation = FALSE
    )
  ),
  list(
    label = "immediate loss after discontinuation",
    overrides = list(
      effect_weight_exponent = 1,
      immediate_loss_after_discontinuation = TRUE
    )
  )
)
structural <- do.call(rbind, lapply(structural_variants, function(x) {
  scenario_metrics(
    calibrate_variant(x$overrides),
    analysis = "structural sensitivity",
    setting = x$label
  )
}))
write.csv(
  structural,
  "results/ejpc/structural_sensitivity.csv",
  row.names = FALSE
)

# Phase-chain convergence at a fixed two-year rebound-memory window.
rebound_memory_years <-
  base_p$n_rebound_stages * base_p$stage_duration_years
convergence <- do.call(rbind, lapply(c(4L, 8L, 12L), function(n_stages) {
  scenario_metrics(
    calibrate_variant(list(
      n_rebound_stages = n_stages,
      stage_duration_years = rebound_memory_years / n_stages,
      effect_weight_exponent = 1,
      immediate_loss_after_discontinuation = FALSE
    )),
    analysis = "numerical convergence",
    setting = paste(n_stages, "rebound stages")
  )
}))
write.csv(
  convergence,
  "results/ejpc/numerical_convergence.csv",
  row.names = FALSE
)

# Trial, primary-policy, and long-horizon projections.
horizons <- c(base_p$select_followup_years, 5, 10)
horizon_results <- do.call(rbind, lapply(horizons, function(horizon) {
  scenario_metrics(
    base_p,
    analysis = "horizon sensitivity",
    setting = sprintf("%.4g-year horizon", horizon),
    horizon = horizon
  )
}))
write.csv(
  horizon_results,
  "results/ejpc/horizon_sensitivity.csv",
  row.names = FALSE
)

# Re-summarize the stored PSA draws; no Monte Carlo or ODE draws are repeated.
psa <- normalize_legacy_psa(readRDS(psa_path), base_p = base_p)
legacy_conversion <- attr(psa, "ejpc_legacy_conversion")
compatible_psa_path <-
  "results/shared/probabilistic_sensitivity_ejpc_compatible.rds"
saveRDS(psa, compatible_psa_path)
psa_summary <- summarize_psa(psa, profile = "ejpc")
write.csv(
  psa_summary,
  "results/ejpc/probabilistic_sensitivity_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    profile = "ejpc",
    n_psa_reused = nrow(psa),
    psa_rerun = FALSE,
    legacy_benchmark_conversion = if (is.null(legacy_conversion)) {
      "not required"
    } else {
      legacy_conversion
    },
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  ),
  "results/ejpc/robustness_run_metadata.csv",
  row.names = FALSE
)

if (nrow(psa) != 10000) {
  warning(
    "The stored PSA object contains ", nrow(psa),
    " draws rather than the 10,000 planned for the final manuscript."
  )
}

if (!is.null(legacy_conversion)) {
  message(legacy_conversion)
}

# Recreate the five EJPC figures from deterministic outputs and stored PSA draws.
source("analysis/05_make_figures.R", local = new.env(parent = globalenv()))

run_tests <- tolower(Sys.getenv("RUN_MODEL_TESTS", unset = "false")) %in%
  c("1", "true", "yes")
if (run_tests) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  Sys.setenv(GLP1_PROJECT_ROOT = project_root)
  testthat::test_dir("tests/testthat", reporter = "summary")
} else {
  message("Optional developer tests skipped; analysis outputs are complete.")
}

message(
  "EJPC robustness analyses and figures are complete. Existing PSA draws were ",
  "reused; analysis/04_probabilistic_sensitivity.R was not run."
)
