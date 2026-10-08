source("R/parameters.R")
source("R/plotting.R")
source("R/profile.R")

p <- readRDS("results/shared/calibrated_parameters.rds")
profiles <- requested_profiles()
simulations <- readRDS("results/shared/core_simulations.rds")
if (!file.exists("results/shared/probabilistic_sensitivity.rds")) {
  stop(
    "Existing PSA draws were not found at ",
    "results/shared/probabilistic_sensitivity.rds. ",
    "Run analysis/04_probabilistic_sensitivity.R once before making figures."
  )
}
compatible_psa_path <-
  "results/shared/probabilistic_sensitivity_ejpc_compatible.rds"
psa <- readRDS(if (
  identical(profiles, "ejpc") && file.exists(compatible_psa_path)
) {
  compatible_psa_path
} else {
  "results/shared/probabilistic_sensitivity.rds"
})

for (profile in profiles) {
  result_dir <- file.path("results", profile)
  figure_dir <- file.path("figures", profile)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  curve <- read.csv(file.path(result_dir, "persistence_response_curve.csv"))
  thresholds_path <- file.path(result_dir, "persistence_thresholds.csv")
  thresholds <- if (file.exists(thresholds_path)) {
    read.csv(thresholds_path)
  } else {
    NULL
  }
  two_way <- read.csv(file.path(result_dir, "two_way_sensitivity.csv"))

  p1 <- plot_treatment_trajectories(
    simulations$real_world_initiators,
    p$n_rebound_stages
  )
  p2 <- plot_mace_trajectories(simulations, profile = profile)
  p3 <- plot_persistence_curve(
    curve,
    profile = profile,
    real_world_persistence_1y = 1 - p$discontinue_probability_1y,
    thresholds = thresholds
  )
  p4 <- plot_psa(psa, profile = profile)
  p5 <- plot_two_way_sensitivity(
    two_way,
    profile = profile,
    observed_discontinuation_1y = p$discontinue_probability_1y,
    observed_restart_1y = p$restart_probability_1y
  )

  ggplot2::ggsave(
    file.path(figure_dir, "figure_1_treatment_states.png"), p1,
    width = 7.5, height = 5.0, dpi = 600, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, "figure_2_mace_trajectories.png"), p2,
    width = 7.5, height = 5.0, dpi = 600, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, "figure_3_persistence_response.png"), p3,
    width = 7.5, height = 5.0, dpi = 600, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, "figure_4_psa.png"), p4,
    width = 7.5, height = 5.0, dpi = 600, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, "figure_5_two_way_sensitivity.png"), p5,
    width = 7.5, height = 5.2, dpi = 600, bg = "white"
  )
}
