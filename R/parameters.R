annual_probability_to_rate <- function(probability) {
  if (any(probability < 0 | probability >= 1)) {
    stop("Annual probabilities must be in [0, 1).")
  }
  -log1p(-probability)
}

read_model_parameters <- function(path = "config/parameters.csv") {
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  values <- setNames(as.list(tab$value), tab$parameter)

  values$lambda_mace <- values$select_mace_events_placebo /
    values$select_mace_py_placebo
  values$lambda_ncv_death <- values$select_ncv_deaths_placebo /
    values$select_mace_py_placebo
  values$fatal_mace_fraction <-
    values$select_fatal_or_undetermined_first_mace /
    values$select_mace_events_placebo
  values$delta <- annual_probability_to_rate(
    values$discontinue_probability_1y
  )
  values$rho <- annual_probability_to_rate(values$restart_probability_1y)
  values$alpha <- annual_probability_to_rate(values$uptake_probability_1y)
  values$rebound_k_year <- values$rebound_k_week * 52.1775
  values$stage_rate <- 1 / values$stage_duration_years
  values$n_rebound_stages <- as.integer(values$n_rebound_stages)
  values$hr_mace_on <- values$select_hr_itt
  values$select_delta <- -log1p(
    -values$select_treatment_discontinuation_total
  ) / values$select_followup_years
  values$select_discontinue_probability_1y <-
    1 - exp(-values$select_delta)
  values
}

rebound_profile <- function(p, restart_feedback = FALSE) {
  stage <- seq_len(p$n_rebound_stages)
  time_since_stop <- (stage - 0.5) / p$stage_rate
  regain_fraction <- p$rebound_plateau *
    (1 - exp(-p$rebound_k_year * time_since_stop))
  remaining_weight_benefit <- pmax(0, 1 - regain_fraction)
  remaining_cv_benefit <- if (
    isTRUE(p$immediate_loss_after_discontinuation)
  ) {
    rep(0, length(stage))
  } else {
    remaining_weight_benefit ^ p$effect_weight_exponent
  }
  hr_mace <- 1 - (1 - p$hr_mace_on) * remaining_cv_benefit
  hr_ncv <- 1 - (1 - p$hr_ncv_death_on) * remaining_cv_benefit

  restart_rate <- rep(p$rho, length(stage))
  if (restart_feedback) {
    regain_at_one_year <- p$rebound_plateau *
      (1 - exp(-p$rebound_k_year))
    regain_pct_body_weight <- p$on_treatment_weight_loss_pct *
      (regain_fraction - regain_at_one_year)
    restart_rate <- p$rho *
      p$restart_hr_per_1pct_regain ^ regain_pct_body_weight
  }

  data.frame(
    stage = stage,
    time_since_stop_years = time_since_stop,
    regain_fraction = regain_fraction,
    remaining_weight_benefit = remaining_weight_benefit,
    remaining_cv_benefit = remaining_cv_benefit,
    hr_mace = hr_mace,
    hr_ncv_death = hr_ncv,
    restart_rate = restart_rate
  )
}
