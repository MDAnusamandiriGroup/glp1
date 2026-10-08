publication_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = base_size + 3, face = "bold", margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1, color = "grey30",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(
        size = base_size - 2, color = "grey35", hjust = 0
      ),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 18, 12, 18)
    )
}

plot_treatment_trajectories <- function(simulation, n_rebound_stages) {
  r_names <- paste0("R", seq_len(n_rebound_stages))
  plot_data <- data.frame(
    time = simulation$time,
    Untreated = simulation$U,
    Treated = simulation$T,
    Discontinued = rowSums(simulation[r_names])
  )
  state_columns <- c("Treated", "Discontinued")
  if (max(abs(plot_data$Untreated), na.rm = TRUE) > 1e-8) {
    state_columns <- c("Untreated", state_columns)
  }
  plot_long <- tidyr::pivot_longer(
    plot_data,
    cols = tidyselect::all_of(state_columns),
    names_to = "state",
    values_to = "people"
  )
  initial_total <- sum(unlist(
    plot_data[1, c("Untreated", "Treated", "Discontinued")],
    use.names = FALSE
  ))
  plot_long$percent <- 100 * plot_long$people / initial_total
  plot_long$state <- factor(
    plot_long$state,
    levels = c("Treated", "Discontinued", "Untreated")
  )

  ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = time, y = percent, color = state)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = c(
      Treated = "#009E73", Discontinued = "#D55E00",
      Untreated = "#0072B2"
    )) +
    ggplot2::labs(
      x = "Years",
      y = "Initial cohort alive and first-MACE-free (%)",
      color = NULL,
      title = "Treatment states after initiation",
      caption = paste(
        "All participants initiate treatment at time zero; absorbing first-MACE",
        "and competing-death states are not shown."
      )
    ) +
    publication_theme()
}

plot_two_way_sensitivity <- function(
    sensitivity, profile = c("ejpc", "dom"),
    observed_discontinuation_1y = NULL,
    observed_restart_1y = NULL) {
  profile <- match.arg(profile)
  p <- ggplot2::ggplot(
    sensitivity,
    ggplot2::aes(
      x = 100 * discontinuation_1y,
      y = 100 * restart_1y,
      fill = 100 * retained_benefit_fraction
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_contour(
      data = sensitivity,
      ggplot2::aes(
        x = 100 * discontinuation_1y,
        y = 100 * restart_1y,
        z = 100 * retained_benefit_fraction
      ),
      inherit.aes = FALSE,
      breaks = c(75, 90, 100), color = "white", linewidth = 0.35,
      alpha = 0.85
    ) +
    ggplot2::scale_fill_viridis_c(
      name = if (profile == "ejpc") {
        "SELECT-like benefit\nretained (%)"
      } else {
        "Full-persistence\nbenefit retained (%)"
      }
    ) +
    ggplot2::labs(
      x = "Discontinuing within one year (%)",
      y = "Re-initiating within one year (%)",
      title = "Discontinuation and re-initiation",
      subtitle = "White contours mark 75%, 90%, and 100% benefit retained",
      caption = if (
        !is.null(observed_discontinuation_1y) &&
          !is.null(observed_restart_1y)
      ) {
        sprintf(
          "Marker: observed estimates (%.1f%% discontinuation; %.1f%% re-initiation).",
          100 * observed_discontinuation_1y,
          100 * observed_restart_1y
        )
      } else {
        NULL
      }
    ) +
    ggplot2::scale_x_continuous(breaks = seq(10, 90, by = 20)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 80, by = 20)) +
    publication_theme() +
    ggplot2::theme(legend.position = "right")

  if (
    !is.null(observed_discontinuation_1y) &&
      !is.null(observed_restart_1y)
  ) {
    observed <- data.frame(
      x = 100 * observed_discontinuation_1y,
      y = 100 * observed_restart_1y
    )
    p <- p + ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21, size = 3.4, stroke = 0.9,
      fill = "white", color = "black"
    )
  }
  p
}

plot_mace_trajectories <- function(simulations, profile = c("ejpc", "dom")) {
  profile <- match.arg(profile)
  keep <- if (profile == "ejpc") {
    c("no_treatment", "real_world_initiators", "select_like")
  } else {
    c("no_treatment", "real_world_initiators", "perfect_persistence")
  }
  simulations <- simulations[keep]
  dat <- do.call(rbind, lapply(simulations, function(x) {
    data.frame(
      time = x$time,
      first_mace = x$M_nonfatal + x$D_fatal_mace,
      scenario = x$scenario
    )
  }))
  dat$scenario <- factor(
    dat$scenario,
    levels = keep,
    labels = if (profile == "ejpc") {
      c("No treatment", "Real-world", "SELECT-like")
    } else {
      c("No treatment", "Real-world persistence", "Full persistence")
    }
  )
  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = time, y = first_mace, color = scenario)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = c(
      "No treatment" = "#D55E00",
      "Real-world" = "#009E73",
      "SELECT-like" = "#0072B2",
      "Real-world persistence" = "#009E73",
      "Full persistence" = "#0072B2"
    )) +
    ggplot2::scale_x_continuous(
      breaks = 0:ceiling(max(dat$time)),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Years", y = "Cumulative first MACE per 100,000",
      color = NULL, title = "Projected first MACE"
    ) +
    publication_theme() +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    )
}

plot_persistence_curve <- function(
    curve, profile = c("ejpc", "dom"),
    real_world_persistence_1y = NULL,
    thresholds = NULL) {
  profile <- match.arg(profile)
  p <- ggplot2::ggplot(
    curve,
    ggplot2::aes(x = 100 * persistence_1y,
                 y = 100 * retained_benefit_fraction)
  ) +
    ggplot2::geom_hline(
      yintercept = c(75, 90, 100), linetype = "dashed", color = "grey65"
    ) +
    ggplot2::geom_line(linewidth = 1.1, color = "#6A3D9A") +
    ggplot2::labs(
      x = "Persistent at one year (%)",
      y = if (profile == "ejpc") {
        "SELECT-like trial benefit retained (%)"
      } else {
        "Full-persistence benefit retained (%)"
      },
      title = "Persistence-response curve",
      caption = if (profile == "ejpc") {
        "Values above 100% exceed the SELECT-like implementation benchmark."
      } else {
        NULL
      }
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 100, by = 20), limits = c(0, 100)
    ) +
    publication_theme()

  if (!is.null(real_world_persistence_1y)) {
    observed_x <- 100 * real_world_persistence_1y
    p <- p +
      ggplot2::geom_vline(
        xintercept = observed_x, linetype = "dotted",
        linewidth = 0.7, color = "#0072B2"
      ) +
      ggplot2::annotate(
        "text", x = observed_x + 1.5,
        y = max(100 * curve$retained_benefit_fraction, na.rm = TRUE),
        label = sprintf("Observed persistence %.1f%%", observed_x),
        hjust = 0, vjust = 1.1, size = 3.2, color = "#0072B2"
      )
  }

  if (!is.null(thresholds)) {
    points <- thresholds[
      is.finite(thresholds$minimum_persistence_1y) &
        thresholds$target_retained_benefit %in% c(0.75, 0.90),
      , drop = FALSE
    ]
    if (nrow(points)) {
      p <- p + ggplot2::geom_point(
        data = points,
        ggplot2::aes(
          x = 100 * minimum_persistence_1y,
          y = 100 * target_retained_benefit
        ),
        inherit.aes = FALSE,
        size = 2.5, color = "#6A3D9A"
      )
    }
  }
  p
}

plot_psa <- function(psa, profile = c("ejpc", "dom")) {
  profile <- match.arg(profile)
  metric <- if (profile == "ejpc") {
    "retained_vs_select_like"
  } else {
    "retained_vs_full_persistence"
  }
  x <- 100 * psa[[metric]]
  interval <- stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
  ggplot2::ggplot(
    psa,
    ggplot2::aes(x = 100 * .data[[metric]])
  ) +
    ggplot2::geom_histogram(
      bins = 40, fill = "#0072B2", color = "white", linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = interval[c(1, 3)], linetype = "dashed",
      linewidth = 0.55, color = "grey35"
    ) +
    ggplot2::geom_vline(
      xintercept = interval[2], linewidth = 0.75, color = "#D55E00"
    ) +
    ggplot2::labs(
      x = if (profile == "ejpc") {
        "SELECT-like benefit retained (%)"
      } else {
        "Full-persistence benefit retained (%)"
      },
      y = "Simulation draws",
      title = "Uncertainty in retained benefit",
      subtitle = sprintf(
        "Median %.1f%% (95%% uncertainty interval %.1f-%.1f%%)",
        interval[2], interval[1], interval[3]
      )
    ) +
    publication_theme()
}
