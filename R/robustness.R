calibrate_variant <- function(overrides = list()) {
  p_i <- read_model_parameters()
  for (name in names(overrides)) {
    p_i[[name]] <- overrides[[name]]
  }
  p_i$n_rebound_stages <- as.integer(p_i$n_rebound_stages)
  p_i$stage_rate <- 1 / p_i$stage_duration_years
  calibrate_select_model(p_i)
}

scenario_metrics <- function(p_i, analysis, setting,
                             horizon = p_i$horizon_years) {
  comparison <- compare_publication_scenarios(p_i, horizon = horizon)$summary
  no_tx <- comparison[comparison$scenario == "no_treatment", ]
  full <- comparison[comparison$scenario == "perfect_persistence", ]
  trial <- comparison[comparison$scenario == "select_like", ]
  real <- comparison[comparison$scenario == "real_world_initiators", ]

  data.frame(
    analysis = analysis,
    setting = setting,
    horizon_years = horizon,
    n_rebound_stages = p_i$n_rebound_stages,
    stage_duration_years = p_i$stage_duration_years,
    effect_weight_exponent = p_i$effect_weight_exponent,
    immediate_loss_after_discontinuation =
      isTRUE(p_i$immediate_loss_after_discontinuation),
    first_mace_no_treatment = no_tx$first_mace,
    first_mace_select_like = trial$first_mace,
    first_mace_full_persistence = full$first_mace,
    first_mace_real_world = real$first_mace,
    mace_prevented_real_world = no_tx$first_mace - real$first_mace,
    excess_mace_vs_select_like = real$first_mace - trial$first_mace,
    retained_vs_select_like =
      (no_tx$first_mace - real$first_mace) /
      (no_tx$first_mace - trial$first_mace),
    retained_vs_full_persistence =
      (no_tx$first_mace - real$first_mace) /
      (no_tx$first_mace - full$first_mace),
    fatal_or_undetermined_first_mace_prevented =
      no_tx$fatal_or_undetermined_first_mace -
      real$fatal_or_undetermined_first_mace,
    stringsAsFactors = FALSE
  )
}

retention_function <- function(p_i, horizon = p_i$horizon_years) {
  baseline <- compare_publication_scenarios(p_i, horizon = horizon)$summary
  no_tx <- baseline$first_mace[baseline$scenario == "no_treatment"]
  trial <- baseline$first_mace[baseline$scenario == "select_like"]
  benchmark_benefit <- no_tx - trial

  function(persistence_1y) {
    p_run <- p_i
    persistence_1y <- min(max(persistence_1y, 1e-8), 1)
    p_run$delta <- -log(persistence_1y)
    real <- summarize_simulation(simulate_model(
      p_run,
      scenario = "real_world_initiators",
      horizon = horizon
    ))
    (no_tx - real$first_mace) / benchmark_benefit
  }
}

refined_persistence_thresholds <- function(
    p_i, targets = c(0.50, 0.75, 0.90), lower_bound = 1e-4) {
  retained <- retention_function(p_i)
  lower_retained <- retained(lower_bound)
  upper_retained <- retained(1)

  do.call(rbind, lapply(targets, function(target) {
    if (lower_retained >= target) {
      threshold <- lower_bound
      status <- "at_or_below_lower_bound"
    } else if (upper_retained < target) {
      threshold <- NA_real_
      status <- "not_reached"
    } else {
      threshold <- stats::uniroot(
        function(x) retained(x) - target,
        interval = c(lower_bound, 1),
        tol = 1e-6
      )$root
      status <- "estimated"
    }
    data.frame(
      target_retained_benefit = target,
      minimum_persistence_1y = threshold,
      threshold_status = status,
      lower_search_bound = lower_bound,
      retained_benefit_at_lower_bound = lower_retained,
      stringsAsFactors = FALSE
    )
  }))
}
