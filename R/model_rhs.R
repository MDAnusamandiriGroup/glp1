state_names <- function(n_rebound_stages) {
  c(
    "U", "T", paste0("R", seq_len(n_rebound_stages)),
    "M_nonfatal", "D_fatal_mace", "D_ncv"
  )
}

initial_state <- function(p, scenario = c(
                            "real_world_initiators", "perfect_persistence",
                            "select_like", "no_treatment", "population_uptake"
                          )) {
  scenario <- match.arg(scenario)
  y <- setNames(rep(0, length(state_names(p$n_rebound_stages))),
                state_names(p$n_rebound_stages))
  if (scenario %in% c(
    "real_world_initiators", "perfect_persistence", "select_like"
  )) {
    y["T"] <- p$cohort_size
  } else {
    y["U"] <- p$cohort_size
  }
  y
}

glp1_rhs <- function(time, state, pars) {
  p <- pars$p
  profile <- pars$profile
  n_r <- p$n_rebound_stages
  r_names <- paste0("R", seq_len(n_r))
  alive_names <- c("U", "T", r_names)
  d <- setNames(rep(0, length(state)), names(state))

  # Treatment dynamics
  uptake <- p$alpha * state["U"]
  d["U"] <- d["U"] - uptake
  d["T"] <- d["T"] + uptake

  discontinuation <- p$delta * state["T"]
  d["T"] <- d["T"] - discontinuation
  d["R1"] <- d["R1"] + discontinuation

  for (j in seq_len(n_r)) {
    rj <- r_names[j]
    restart <- profile$restart_rate[j] * state[rj]
    d[rj] <- d[rj] - restart
    d["T"] <- d["T"] + restart

    if (j < n_r) {
      progression <- p$stage_rate * state[rj]
      d[rj] <- d[rj] - progression
      d[r_names[j + 1]] <- d[r_names[j + 1]] + progression
    }
  }

  # First-event and competing-risk dynamics
  hr_mace <- c(U = 1, T = p$hr_mace_on,
               setNames(profile$hr_mace, r_names))
  hr_ncv <- c(U = 1, T = p$hr_ncv_death_on,
              setNames(profile$hr_ncv_death, r_names))

  for (compartment in alive_names) {
    n_people <- state[compartment]
    mace_flow <- p$lambda_mace * hr_mace[compartment] * n_people
    fatal_flow <- p$fatal_mace_fraction * mace_flow
    nonfatal_flow <- (1 - p$fatal_mace_fraction) * mace_flow
    ncv_flow <- p$lambda_ncv_death * hr_ncv[compartment] * n_people

    d[compartment] <- d[compartment] -
      fatal_flow - nonfatal_flow - ncv_flow
    d["M_nonfatal"] <- d["M_nonfatal"] + nonfatal_flow
    d["D_fatal_mace"] <- d["D_fatal_mace"] + fatal_flow
    d["D_ncv"] <- d["D_ncv"] + ncv_flow
  }

  list(d)
}
