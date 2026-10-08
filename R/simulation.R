simulate_model <- function(p,
                           scenario = c(
                             "real_world_initiators", "perfect_persistence",
                             "select_like", "no_treatment", "population_uptake"
                           ),
                           horizon = p$horizon_years,
                           dt = p$time_step_years,
                           restart_feedback = FALSE) {
  scenario <- match.arg(scenario)
  p_run <- p

  if (scenario == "perfect_persistence") {
    p_run$delta <- 0
    p_run$alpha <- 0
  } else if (scenario == "select_like") {
    p_run$delta <- p_run$select_delta
    p_run$rho <- 0
    p_run$alpha <- 0
  } else if (scenario == "no_treatment") {
    p_run$delta <- 0
    p_run$alpha <- 0
  } else if (scenario == "real_world_initiators") {
    p_run$alpha <- 0
  }

  profile <- rebound_profile(p_run, restart_feedback = restart_feedback)
  times <- unique(c(seq(0, horizon, by = dt), horizon))
  y0 <- initial_state(p_run, scenario)

  out <- deSolve::ode(
    y = y0,
    times = times,
    func = glp1_rhs,
    parms = list(p = p_run, profile = profile),
    method = "lsoda",
    atol = 1e-8,
    rtol = 1e-8
  )
  out <- as.data.frame(out)
  out$scenario <- scenario
  out
}

compare_publication_scenarios <- function(
    p, horizon = p$horizon_years, restart_feedback = FALSE) {
  names_to_run <- c(
    "no_treatment", "perfect_persistence", "select_like",
    "real_world_initiators"
  )
  simulations <- lapply(names_to_run, function(s) {
    simulate_model(
      p, scenario = s, horizon = horizon,
      restart_feedback = restart_feedback
    )
  })
  names(simulations) <- names_to_run
  summary <- do.call(rbind, lapply(simulations, summarize_simulation))
  rownames(summary) <- NULL

  no_tx <- summary$first_mace[summary$scenario == "no_treatment"]
  full <- summary$first_mace[summary$scenario == "perfect_persistence"]
  trial <- summary$first_mace[summary$scenario == "select_like"]

  summary$mace_prevented_vs_no_treatment <- no_tx - summary$first_mace
  summary$retained_vs_full_persistence <-
    (no_tx - summary$first_mace) / (no_tx - full)
  summary$retained_vs_select_like <-
    (no_tx - summary$first_mace) / (no_tx - trial)

  list(summary = summary, simulations = simulations)
}

profile_scenario_summary <- function(comparison, profile = c("ejpc", "dom")) {
  profile <- match.arg(profile)
  out <- comparison$summary
  out$retained_benefit_fraction <- if (profile == "ejpc") {
    out$retained_vs_select_like
  } else {
    out$retained_vs_full_persistence
  }
  out$benchmark <- if (profile == "ejpc") {
    "SELECT-like trial implementation"
  } else {
    "Full treatment persistence"
  }
  out
}

summarize_simulation <- function(simulation) {
  last <- simulation[nrow(simulation), , drop = FALSE]
  data.frame(
    scenario = last$scenario,
    horizon_years = last$time,
    first_mace = last$M_nonfatal + last$D_fatal_mace,
    nonfatal_first_mace = last$M_nonfatal,
    fatal_or_undetermined_first_mace = last$D_fatal_mace,
    non_cv_death_before_first_mace = last$D_ncv,
    event_free_alive = sum(unlist(
      last[c("U", "T", grep("^R[0-9]+$", names(last), value = TRUE))],
      use.names = FALSE
    )),
    total_population = sum(unlist(
      last[state_names(
        length(grep("^R[0-9]+$", names(last), value = TRUE))
      )],
      use.names = FALSE
    ))
  )
}

two_way_sensitivity <- function(p, horizon = p$horizon_years,
                                discontinuation_grid = seq(0.1, 0.9, by = 0.1),
                                restart_grid = seq(0, 0.8, by = 0.1),
                                benchmark = c("full", "select")) {
  benchmark <- match.arg(benchmark)
  baseline <- compare_publication_scenarios(p, horizon = horizon)$summary
  no_tx <- baseline$first_mace[baseline$scenario == "no_treatment"]
  comparator_name <- if (benchmark == "select") {
    "select_like"
  } else {
    "perfect_persistence"
  }
  comparator <- baseline$first_mace[baseline$scenario == comparator_name]
  maximum_benefit <- no_tx - comparator

  grid <- expand.grid(
    discontinuation_1y = discontinuation_grid,
    restart_1y = restart_grid
  )
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    p_i <- p
    p_i$delta <- annual_probability_to_rate(grid$discontinuation_1y[i])
    p_i$rho <- annual_probability_to_rate(grid$restart_1y[i])
    sim <- simulate_model(
      p_i, scenario = "real_world_initiators", horizon = horizon
    )
    outcome <- summarize_simulation(sim)
    data.frame(
      discontinuation_1y = grid$discontinuation_1y[i],
      restart_1y = grid$restart_1y[i],
      first_mace = outcome$first_mace,
      retained_benefit_fraction =
        (no_tx - outcome$first_mace) / maximum_benefit
    )
  })
  do.call(rbind, rows)
}

compare_core_scenarios <- function(p, horizon = p$horizon_years,
                                   restart_feedback = FALSE) {
  names_to_run <- c(
    "no_treatment", "perfect_persistence", "real_world_initiators"
  )
  simulations <- lapply(names_to_run, function(s) {
    simulate_model(
      p, scenario = s, horizon = horizon,
      restart_feedback = restart_feedback
    )
  })
  names(simulations) <- names_to_run
  summary <- do.call(rbind, lapply(simulations, summarize_simulation))
  rownames(summary) <- NULL

  no_tx <- summary$first_mace[summary$scenario == "no_treatment"]
  perfect <- summary$first_mace[summary$scenario == "perfect_persistence"]
  real <- summary$first_mace[summary$scenario == "real_world_initiators"]
  denominator <- no_tx - perfect
  retained <- if (denominator > 0) (no_tx - real) / denominator else NA_real_

  summary$mace_prevented_vs_no_treatment <- no_tx - summary$first_mace
  summary$retained_benefit_fraction <- ifelse(
    summary$scenario == "real_world_initiators", retained,
    ifelse(summary$scenario == "perfect_persistence", 1, 0)
  )

  list(summary = summary, simulations = simulations)
}

persistence_response_curve <- function(p, horizon = p$horizon_years,
                                       persistence_grid = seq(0.05, 1, by = 0.01),
                                       restart_feedback = FALSE,
                                       benchmark = c("full", "select")) {
  benchmark <- match.arg(benchmark)
  baseline <- compare_publication_scenarios(
    p, horizon = horizon, restart_feedback = restart_feedback
  )$summary
  no_tx <- baseline$first_mace[baseline$scenario == "no_treatment"]
  comparator_name <- if (benchmark == "select") {
    "select_like"
  } else {
    "perfect_persistence"
  }
  comparator <- baseline$first_mace[baseline$scenario == comparator_name]
  maximum_benefit <- no_tx - comparator

  rows <- lapply(persistence_grid, function(persist_1y) {
    p_i <- p
    p_i$delta <- -log(persist_1y)
    sim <- simulate_model(
      p_i,
      scenario = "real_world_initiators",
      horizon = horizon,
      restart_feedback = restart_feedback
    )
    s <- summarize_simulation(sim)
    data.frame(
      persistence_1y = persist_1y,
      discontinuation_1y = 1 - persist_1y,
      first_mace = s$first_mace,
      mace_prevented = no_tx - s$first_mace,
      retained_benefit_fraction = (no_tx - s$first_mace) / maximum_benefit
    )
  })
  do.call(rbind, rows)
}

persistence_thresholds <- function(curve, targets = c(0.50, 0.75, 0.90)) {
  do.call(rbind, lapply(targets, function(target) {
    eligible <- curve[curve$retained_benefit_fraction >= target, , drop = FALSE]
    data.frame(
      target_retained_benefit = target,
      minimum_persistence_1y = if (nrow(eligible)) {
        min(eligible$persistence_1y)
      } else {
        NA_real_
      }
    )
  }))
}
