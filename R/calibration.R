select_trial_discontinuation_rate <- function(p) {
  -log1p(-p$select_treatment_discontinuation_total) /
    p$select_followup_years
}

calibrate_select_model <- function(p) {
  p$lambda_mace_direct <- p$lambda_mace
  p$select_delta <- select_trial_discontinuation_rate(p)
  p$select_discontinue_probability_1y <- 1 - exp(-p$select_delta)

  placebo_target <- p$select_mace_events_placebo / p$select_n_placebo
  semaglutide_target <-
    p$select_mace_events_semaglutide / p$select_n_semaglutide

  placebo_objective <- function(lambda) {
    p_try <- p
    p_try$lambda_mace <- lambda
    sim <- simulate_model(
      p_try, scenario = "no_treatment",
      horizon = p_try$select_followup_years
    )
    summarize_simulation(sim)$first_mace / p_try$cohort_size -
      placebo_target
  }

  lambda_bounds <- c(
    max(1e-6, p$lambda_mace_direct * 0.25),
    p$lambda_mace_direct * 4
  )
  p$lambda_mace <- stats::uniroot(
    placebo_objective,
    interval = lambda_bounds,
    tol = 1e-10
  )$root
  p$lambda_mace_calibration_factor <-
    p$lambda_mace / p$lambda_mace_direct

  treatment_objective <- function(hr_on) {
    p_try <- p
    p_try$hr_mace_on <- hr_on
    sim <- simulate_model(
      p_try, scenario = "select_like",
      horizon = p_try$select_followup_years
    )
    summarize_simulation(sim)$first_mace / p_try$cohort_size -
      semaglutide_target
  }

  p$hr_mace_on <- stats::uniroot(
    treatment_objective,
    interval = c(0.20, 0.999),
    tol = 1e-10
  )$root
  p$hr_mace_on_calibrated <- p$hr_mace_on
  p
}

calibration_table <- function(p) {
  placebo <- summarize_simulation(simulate_model(
    p, "no_treatment", horizon = p$select_followup_years
  ))
  treatment <- summarize_simulation(simulate_model(
    p, "select_like", horizon = p$select_followup_years
  ))

  placebo_model <- placebo$first_mace / p$cohort_size
  treatment_model <- treatment$first_mace / p$cohort_size
  placebo_observed <- p$select_mace_events_placebo / p$select_n_placebo
  treatment_observed <-
    p$select_mace_events_semaglutide / p$select_n_semaglutide

  data.frame(
    metric = c(
      "placebo_first_mace_risk",
      "semaglutide_first_mace_risk",
      "crude_risk_ratio",
      "calibrated_on_treatment_hr",
      "trial_discontinuation_probability_total",
      "trial_discontinuation_probability_1y",
      "calibrated_placebo_mace_hazard"
    ),
    model_value = c(
      placebo_model,
      treatment_model,
      treatment_model / placebo_model,
      p$hr_mace_on,
      p$select_treatment_discontinuation_total,
      p$select_discontinue_probability_1y,
      p$lambda_mace
    ),
    target_or_reference = c(
      placebo_observed,
      treatment_observed,
      treatment_observed / placebo_observed,
      p$select_hr_itt,
      p$select_treatment_discontinuation_total,
      NA_real_,
      p$lambda_mace_direct
    ),
    interpretation = c(
      "calibration target",
      "calibration target",
      "derived validation",
      "latent state effect; compare cautiously with ITT HR",
      "calibration input over the mean trial follow-up",
      "annualized from total trial discontinuation",
      "calibrated value; reference is direct event/person-year rate"
    ),
    stringsAsFactors = FALSE
  )
}
