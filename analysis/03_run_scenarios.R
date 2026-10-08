source("R/parameters.R")
source("R/model_rhs.R")
source("R/simulation.R")
source("R/profile.R")

p <- readRDS("results/shared/calibrated_parameters.rds")
profiles <- requested_profiles()
core <- compare_publication_scenarios(p)
feedback_core <- compare_publication_scenarios(p, restart_feedback = TRUE)

saveRDS(core$simulations, "results/shared/core_simulations.rds")
saveRDS(feedback_core$simulations, "results/shared/feedback_simulations.rds")

for (profile in profiles) {
  output_dir <- file.path("results", profile)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  benchmark <- profile_benchmark(profile)

  core_summary <- profile_scenario_summary(core, profile)
  core_summary$analysis <- "constant restart primary"
  feedback_summary <- profile_scenario_summary(feedback_core, profile)
  feedback_summary$analysis <-
    "weight-regain-dependent restart sensitivity"

  curve <- persistence_response_curve(p, benchmark = benchmark)
  thresholds <- persistence_thresholds(curve)
  two_way <- two_way_sensitivity(p, benchmark = benchmark)

  write.csv(
    rbind(core_summary, feedback_summary),
    file.path(output_dir, "scenario_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    curve, file.path(output_dir, "persistence_response_curve.csv"),
    row.names = FALSE
  )
  write.csv(
    thresholds, file.path(output_dir, "persistence_thresholds.csv"),
    row.names = FALSE
  )
  write.csv(
    two_way, file.path(output_dir, "two_way_sensitivity.csv"),
    row.names = FALSE
  )
}
