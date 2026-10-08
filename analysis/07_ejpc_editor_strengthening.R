options(stringsAsFactors = FALSE)

# This script is intentionally post-processing only:
# - it does NOT rerun the 10,000-draw PSA;
# - it does NOT overwrite any existing result;
# - it converts completed robustness outputs into editor-facing tables/figure.

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing package: ", paste(missing_packages, collapse = ", "),
    ". Run source('analysis/00_install_packages.R') first."
  )
}

required_inputs <- c(
  "results/ejpc/structural_sensitivity.csv",
  "config/parameters.csv",
  "data/derived/nhanes_select_like_population.csv"
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Required completed outputs are missing:\n- ",
    paste(missing_inputs, collapse = "\n- "),
    "\nThis public release needs the listed inputs. Do not rerun the PSA merely to regenerate these tables."
  )
}

dir.create("results/ejpc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/ejpc", recursive = TRUE, showWarnings = FALSE)

structural <- read.csv(
  "results/ejpc/structural_sensitivity.csv",
  check.names = FALSE
)

wanted <- c(
  "gradual waning; eta=0.5",
  "gradual waning; eta=1.0 (base case)",
  "gradual waning; eta=1.5",
  "immediate loss after discontinuation"
)
if (!all(wanted %in% structural$setting)) {
  stop(
    "structural_sensitivity.csv does not contain all four prespecified ",
    "waning scenarios. Rerun analysis/06_ejpc_robustness.R; PSA is reused."
  )
}

matrix <- structural[match(wanted, structural$setting), ]
matrix$scenario_label <- c(
  "Slower waning (eta=0.5)",
  "Weight-linked waning (base case; eta=1.0)",
  "Faster waning (eta=1.5)",
  "Immediate loss of protection"
)
matrix$structural_interpretation <- c(
  "More cardiovascular protection retained for a given weight-regain fraction",
  "Cardiovascular protection wanes proportionally with remaining weight benefit",
  "Cardiovascular protection wanes faster than remaining weight benefit",
  "No cardiovascular protection remains after discontinuation"
)
matrix$retained_vs_select_like_percent <-
  100 * matrix$retained_vs_select_like
matrix$incremental_loss_vs_select_like_percent <-
  100 * (1 - matrix$retained_vs_select_like)

editor_matrix <- matrix[, c(
  "scenario_label",
  "structural_interpretation",
  "first_mace_real_world",
  "mace_prevented_real_world",
  "excess_mace_vs_select_like",
  "retained_vs_select_like_percent",
  "incremental_loss_vs_select_like_percent"
)]

write.csv(
  editor_matrix,
  "results/ejpc/editor_waning_scenario_matrix.csv",
  row.names = FALSE
)

editor_matrix$scenario_label <- factor(
  editor_matrix$scenario_label,
  levels = rev(editor_matrix$scenario_label)
)
editor_matrix$is_base <- grepl("base case", editor_matrix$scenario_label)

p <- ggplot2::ggplot(
  editor_matrix,
  ggplot2::aes(
    x = scenario_label,
    y = retained_vs_select_like_percent,
    fill = is_base
  )
) +
  ggplot2::geom_col(width = 0.67, show.legend = FALSE) +
  ggplot2::geom_hline(
    yintercept = 100,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.45
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.1f%%", retained_vs_select_like_percent)
    ),
    hjust = -0.12,
    size = 3.7
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "#9CB8CB", `TRUE` = "#005B82")) +
  ggplot2::scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "SELECT-like cardiovascular benefit retained (%)",
    title = "Structural uncertainty in post-discontinuation waning",
    subtitle = "Five-year outcomes; all scenarios recalibrated to both SELECT trial arms"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(8, 18, 8, 8)
  )

ggplot2::ggsave(
  "figures/ejpc/supplementary_figure_s2_waning_scenarios.png",
  p,
  width = 8.0,
  height = 4.8,
  dpi = 600,
  bg = "white"
)
ggplot2::ggsave(
  "figures/ejpc/supplementary_figure_s2_waning_scenarios.tiff",
  p,
  width = 8.0,
  height = 4.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

parameter_table <- read.csv(
  "config/parameters.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
population <- read.csv(
  "data/derived/nhanes_select_like_population.csv",
  check.names = FALSE
)

get_parameter <- function(parameter) {
  value <- parameter_table$value[parameter_table$parameter == parameter]
  if (length(value) != 1L) {
    stop("Model parameter is missing or duplicated: ", parameter)
  }
  as.numeric(value)
}

placebo_risk <-
  get_parameter("select_mace_events_placebo") /
  get_parameter("select_n_placebo")
semaglutide_risk <-
  get_parameter("select_mace_events_semaglutide") /
  get_parameter("select_n_semaglutide")
crude_risk_ratio <- semaglutide_risk / placebo_risk

nhanes_total <- population$estimate[
  population$metric == "eligible_population_total"
]
nhanes_se <- population$standard_error[
  population$metric == "eligible_population_total"
]
if (length(nhanes_total) != 1L || length(nhanes_se) != 1L) {
  stop("NHANES population output is missing or duplicated.")
}

validation <- data.frame(
  validation_type = c(
    "Calibration target",
    "Calibration target",
    "Derived calibration check",
    "Independent population face-validity",
    "Independent population face-validity"
  ),
  quantity = c(
    "First-MACE risk: SELECT placebo arm (39.8 months)",
    "First-MACE risk: SELECT semaglutide arm (39.8 months)",
    "Crude SELECT risk ratio (not the time-to-event HR)",
    "US SELECT-eligible population: Lusk et al.",
    "US SELECT-eligible population: Yan et al."
  ),
  model_estimate = c(
    placebo_risk,
    semaglutide_risk,
    crude_risk_ratio,
    nhanes_total,
    nhanes_total
  ),
  published_benchmark = c(
    placebo_risk,
    semaglutide_risk,
    crude_risk_ratio,
    4473681,
    3700000
  ),
  benchmark_source = c(
    "SELECT",
    "SELECT",
    "SELECT",
    "Lusk et al., Am Heart J 2024; PMID 39182939",
    "Yan et al., JACC Adv 2025; PMID 40402121"
  ),
  interpretation = c(
    "Exact agreement expected because this is a calibration target",
    "Exact agreement expected because this is a calibration target",
    "Derived from the two calibrated arm risks; not external validation",
    "Independent estimate using different NHANES cycles and eligibility coding",
    "Independent estimate using different NHANES cycles and eligibility coding"
  ),
  stringsAsFactors = FALSE
)
validation$absolute_difference <-
  validation$model_estimate - validation$published_benchmark
validation$relative_difference_percent <-
  100 * validation$absolute_difference / validation$published_benchmark
validation$internal_nhanes_standard_error <- c(NA, NA, NA, nhanes_se, nhanes_se)

write.csv(
  validation,
  "results/ejpc/editor_calibration_and_face_validity.csv",
  row.names = FALSE
)

key_messages <- data.frame(
  item = c(
    "Primary interpretation",
    "Implementation gap",
    "Structural range",
    "Persistence-response interpretation",
    "Boundary of inference"
  ),
  manuscript_language = c(
    "Real-world treatment trajectories retained approximately three-quarters of SELECT-like cardiovascular benefit over five years.",
    "Discontinuation and reinitiation produced a clinically meaningful gap relative to a SELECT-like implementation.",
    sprintf(
      "Across prespecified post-discontinuation waning assumptions, retained benefit ranged from %.1f%% to %.1f%%.",
      min(editor_matrix$retained_vs_select_like_percent),
      max(editor_matrix$retained_vs_select_like_percent)
    ),
    "Preserving the final increment from 75% to 90% of SELECT-like benefit required more than doubling one-year persistence.",
    "Results are model projections informed by trial and observational inputs, not causal estimates from linked individual-level outcome data."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  key_messages,
  "results/ejpc/editor_key_messages.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    script = "analysis/07_ejpc_editor_strengthening.R",
    psa_rerun = FALSE,
    ode_rerun = FALSE,
    existing_structural_results_reused = TRUE,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  ),
  "results/ejpc/editor_strengthening_metadata.csv",
  row.names = FALSE
)

message(
  "EJPC editor-strengthening outputs completed. No PSA draws or ODE ",
  "simulations were rerun."
)
