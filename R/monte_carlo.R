sample_lognormal_ci <- function(point, lower, upper) {
  sd_log <- (log(upper) - log(lower)) / (2 * 1.96)
  rlnorm(1, meanlog = log(point), sdlog = sd_log)
}

sample_logitnormal_ci <- function(point, lower, upper) {
  sd_logit <- (qlogis(upper) - qlogis(lower)) / (2 * 1.96)
  plogis(rnorm(1, mean = qlogis(point), sd = sd_logit))
}

parameter_row <- function(parameter_table, parameter) {
  parameter_table[parameter_table$parameter == parameter, , drop = FALSE]
}

sample_model_parameters <- function(base_p, parameter_table) {
  p <- base_p

  hr <- parameter_row(parameter_table, "select_hr_itt")
  p$select_hr_draw <- sample_lognormal_ci(
    hr$value, hr$lower_95, hr$upper_95
  )
  benefit_scaler <-
    (1 - p$select_hr_draw) / (1 - base_p$select_hr_itt)
  p$hr_mace_on <- 1 -
    (1 - base_p$hr_mace_on_calibrated) * benefit_scaler
  p$hr_mace_on <- min(max(p$hr_mace_on, 0.05), 0.999)

  disc <- parameter_row(parameter_table, "discontinue_probability_1y")
  p$discontinue_probability_1y <- sample_logitnormal_ci(
    disc$value, disc$lower_95, disc$upper_95
  )
  p$delta <- annual_probability_to_rate(p$discontinue_probability_1y)

  restart <- parameter_row(parameter_table, "restart_probability_1y")
  p$restart_probability_1y <- sample_logitnormal_ci(
    restart$value, restart$lower_95, restart$upper_95
  )
  p$rho <- annual_probability_to_rate(p$restart_probability_1y)

  plateau <- parameter_row(parameter_table, "rebound_plateau")
  p$rebound_plateau <- sample_logitnormal_ci(
    plateau$value, plateau$lower_95, plateau$upper_95
  )

  rebound_k <- parameter_row(parameter_table, "rebound_k_week")
  p$rebound_k_week <- sample_lognormal_ci(
    rebound_k$value, rebound_k$lower_95, rebound_k$upper_95
  )
  p$rebound_k_year <- p$rebound_k_week * 52.1775

  p$lambda_mace <- rgamma(
    1,
    shape = p$select_mace_events_placebo + 0.5,
    rate = p$select_mace_py_placebo
  ) * base_p$lambda_mace_calibration_factor
  p$lambda_ncv_death <- rgamma(
    1,
    shape = p$select_ncv_deaths_placebo + 0.5,
    rate = p$select_mace_py_placebo
  )
  p$fatal_mace_fraction <- rbeta(
    1,
    p$select_fatal_or_undetermined_first_mace + 0.5,
    p$select_mace_events_placebo -
      p$select_fatal_or_undetermined_first_mace + 0.5
  )

  # Structural rather than sampling uncertainty.
  p$effect_weight_exponent <- runif(1, 0.5, 1.5)
  p
}

run_probabilistic_sensitivity <- function(n, base_p, parameter_table,
                                          seed = 20260928,
                                          restart_feedback = FALSE) {
  set.seed(seed)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    p_i <- sample_model_parameters(base_p, parameter_table)
    comparison <- compare_publication_scenarios(
      p_i,
      horizon = p_i$horizon_years,
      restart_feedback = restart_feedback
    )$summary
    real <- comparison[comparison$scenario == "real_world_initiators", ]
    full <- comparison[comparison$scenario == "perfect_persistence", ]
    trial <- comparison[comparison$scenario == "select_like", ]
    no_tx <- comparison[comparison$scenario == "no_treatment", ]

    rows[[i]] <- data.frame(
      draw = i,
      first_mace_real_world = real$first_mace,
      first_mace_full_persistence = full$first_mace,
      first_mace_select_like = trial$first_mace,
      first_mace_no_treatment = no_tx$first_mace,
      mace_prevented_real_world = no_tx$first_mace - real$first_mace,
      retained_vs_full_persistence =
        (no_tx$first_mace - real$first_mace) /
        (no_tx$first_mace - full$first_mace),
      retained_vs_select_like =
        (no_tx$first_mace - real$first_mace) /
        (no_tx$first_mace - trial$first_mace),
      fatal_or_undetermined_first_mace_no_treatment =
        no_tx$fatal_or_undetermined_first_mace,
      fatal_or_undetermined_first_mace_real_world =
        real$fatal_or_undetermined_first_mace,
      fatal_or_undetermined_first_mace_prevented =
        no_tx$fatal_or_undetermined_first_mace -
        real$fatal_or_undetermined_first_mace,
      # Retained for compatibility with PSA objects generated before v0.2.
      fatal_first_mace_real_world = real$fatal_or_undetermined_first_mace,
      hr_mace_on = p_i$hr_mace_on,
      select_hr_draw = p_i$select_hr_draw,
      discontinuation_1y = p_i$discontinue_probability_1y,
      restart_1y = p_i$restart_probability_1y,
      rebound_plateau = p_i$rebound_plateau,
      rebound_k_week = p_i$rebound_k_week,
      effect_weight_exponent = p_i$effect_weight_exponent
    )
  }
  do.call(rbind, rows)
}

add_derived_psa_outcomes <- function(psa) {
  target <- "fatal_or_undetermined_first_mace_prevented"
  if (target %in% names(psa)) {
    return(psa)
  }

  if (all(c(
    "fatal_or_undetermined_first_mace_no_treatment",
    "fatal_or_undetermined_first_mace_real_world"
  ) %in% names(psa))) {
    psa[[target]] <-
      psa$fatal_or_undetermined_first_mace_no_treatment -
      psa$fatal_or_undetermined_first_mace_real_world
    return(psa)
  }

  legacy_fields <- c(
    "fatal_first_mace_real_world", "first_mace_real_world",
    "mace_prevented_real_world"
  )
  if (!all(legacy_fields %in% names(psa))) {
    stop(
      "PSA object lacks the fields required to derive fatal-or-undetermined ",
      "first-MACE events prevented."
    )
  }

  fatal_fraction <- ifelse(
    psa$first_mace_real_world > 0,
    psa$fatal_first_mace_real_world / psa$first_mace_real_world,
    NA_real_
  )
  psa[[target]] <- psa$mace_prevented_real_world * fatal_fraction
  psa
}

normalize_legacy_psa <- function(psa, base_p = NULL) {
  if (
    !"first_mace_full_persistence" %in% names(psa) &&
      "first_mace_perfect" %in% names(psa)
  ) {
    psa$first_mace_full_persistence <- psa$first_mace_perfect
  }

  if (
    !"retained_vs_full_persistence" %in% names(psa) &&
      "retained_benefit_fraction" %in% names(psa)
  ) {
    psa$retained_vs_full_persistence <- psa$retained_benefit_fraction
  }

  if (!"retained_vs_select_like" %in% names(psa)) {
    if (is.null(base_p)) {
      stop(
        "Legacy PSA requires calibrated parameters to convert the retained-",
        "benefit benchmark to SELECT-like implementation."
      )
    }
    comparison <- compare_publication_scenarios(base_p)$summary
    no_tx <- comparison$first_mace[comparison$scenario == "no_treatment"]
    full <- comparison$first_mace[
      comparison$scenario == "perfect_persistence"
    ]
    trial <- comparison$first_mace[comparison$scenario == "select_like"]
    full_to_select_scale <- (no_tx - full) / (no_tx - trial)
    psa$retained_vs_select_like <-
      psa$retained_vs_full_persistence * full_to_select_scale
    attr(psa, "ejpc_legacy_conversion") <- paste0(
      "Post-hoc conversion from the stored full-persistence retention ratio ",
      "using deterministic full-to-SELECT benefit scale=",
      signif(full_to_select_scale, 8), ". No PSA draws rerun."
    )
  }
  psa
}

summarize_psa <- function(psa, profile = c("ejpc", "dom")) {
  profile <- match.arg(profile)
  psa <- add_derived_psa_outcomes(psa)
  retention_metric <- if (profile == "ejpc") {
    "retained_vs_select_like"
  } else {
    "retained_vs_full_persistence"
  }
  metrics <- c(
    "mace_prevented_real_world", retention_metric,
    "fatal_or_undetermined_first_mace_prevented"
  )
  missing_metrics <- setdiff(metrics, names(psa))
  if (length(missing_metrics)) {
    stop("PSA object is missing metric(s): ", paste(missing_metrics, collapse = ", "))
  }
  out <- do.call(rbind, lapply(metrics, function(metric) {
    x <- psa[[metric]]
    data.frame(
      metric = metric,
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      lower_95 = unname(quantile(x, 0.025, na.rm = TRUE)),
      upper_95 = unname(quantile(x, 0.975, na.rm = TRUE))
    )
  }))
  out$profile <- profile
  out
}
