project_root <- Sys.getenv("GLP1_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  project_root <- normalizePath(file.path(getwd(), "..", ".."))
}
source(file.path(project_root, "R", "parameters.R"))
source(file.path(project_root, "R", "model_rhs.R"))
source(file.path(project_root, "R", "simulation.R"))
source(file.path(project_root, "R", "calibration.R"))
source(file.path(project_root, "R", "profile.R"))
source(file.path(project_root, "R", "monte_carlo.R"))

testthat::test_that("population mass is conserved", {
  testthat::skip_if_not_installed("deSolve")
  p <- calibrate_select_model(read_model_parameters())
  sim <- simulate_model(p, "real_world_initiators", horizon = 2)
  state_columns <- state_names(p$n_rebound_stages)
  totals <- rowSums(sim[state_columns])
  testthat::expect_lt(max(abs(totals - p$cohort_size)), 1e-3)
})

testthat::test_that("full persistence prevents first MACE", {
  testthat::skip_if_not_installed("deSolve")
  p <- calibrate_select_model(read_model_parameters())
  comparison <- compare_core_scenarios(p, horizon = 3)$summary
  no_tx <- comparison$first_mace[comparison$scenario == "no_treatment"]
  perfect <- comparison$first_mace[
    comparison$scenario == "perfect_persistence"
  ]
  testthat::expect_lt(perfect, no_tx)
})

testthat::test_that("higher persistence does not reduce retained benefit", {
  testthat::skip_if_not_installed("deSolve")
  p <- calibrate_select_model(read_model_parameters())
  curve <- persistence_response_curve(
    p, horizon = 3, persistence_grid = seq(0.1, 1, by = 0.1)
  )
  testthat::expect_true(all(diff(curve$retained_benefit_fraction) >= -1e-7))
})

testthat::test_that("immediate-loss sensitivity removes residual CV benefit", {
  p <- read_model_parameters()
  p$immediate_loss_after_discontinuation <- TRUE
  profile <- rebound_profile(p)
  testthat::expect_equal(profile$remaining_cv_benefit, rep(0, nrow(profile)))
  testthat::expect_equal(profile$hr_mace, rep(1, nrow(profile)))
})

testthat::test_that("legacy PSA outcomes are derived without new draws", {
  legacy <- data.frame(
    fatal_first_mace_real_world = c(20, 30),
    first_mace_real_world = c(100, 120),
    mace_prevented_real_world = c(10, 12)
  )
  derived <- add_derived_psa_outcomes(legacy)
  testthat::expect_equal(
    derived$fatal_or_undetermined_first_mace_prevented,
    c(2, 3)
  )
})

testthat::test_that("SELECT calibration reproduces both trial arms", {
  testthat::skip_if_not_installed("deSolve")
  p <- calibrate_select_model(read_model_parameters())
  validation <- calibration_table(p)
  trial_rows <- validation$metric %in% c(
    "placebo_first_mace_risk", "semaglutide_first_mace_risk"
  )
  testthat::expect_lt(max(abs(
    validation$model_value[trial_rows] -
      validation$target_or_reference[trial_rows]
  )), 1e-7)
  testthat::expect_lt(p$hr_mace_on, p$select_hr_itt)
})

testthat::test_that("publication profiles use different valid benchmarks", {
  testthat::skip_if_not_installed("deSolve")
  p <- calibrate_select_model(read_model_parameters())
  comparison <- compare_publication_scenarios(p, horizon = 5)
  ejpc <- profile_scenario_summary(comparison, "ejpc")
  dom <- profile_scenario_summary(comparison, "dom")
  real_row <- comparison$summary$scenario == "real_world_initiators"
  testthat::expect_true(
    ejpc$retained_benefit_fraction[real_row] >=
      dom$retained_benefit_fraction[real_row]
  )
})

testthat::test_that("deterministic base-case checksum is stable", {
  testthat::skip_if_not_installed("deSolve")
  p <- calibrate_select_model(read_model_parameters())
  comparison <- compare_publication_scenarios(p, horizon = 5)$summary
  expected <- c(
    no_treatment = 11700.0353,
    perfect_persistence = 9356.5953,
    select_like = 9634.4828,
    real_world_initiators = 10106.9242
  )
  observed <- setNames(comparison$first_mace, comparison$scenario)
  testthat::expect_equal(
    unname(observed[names(expected)]), unname(expected), tolerance = 0.05
  )
})
