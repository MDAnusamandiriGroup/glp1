options(stringsAsFactors = FALSE)

profiles <- tolower(Sys.getenv("PUBLICATION_PROFILE", unset = "both"))
if (!profiles %in% c("ejpc", "dom", "both")) {
  stop("PUBLICATION_PROFILE must be 'ejpc', 'dom', or 'both'.")
}

required_directories <- c(
  "data/raw", "data/derived", "results/shared", "figures"
)
if (profiles %in% c("ejpc", "both")) {
  required_directories <- c(
    required_directories, "results/ejpc", "figures/ejpc"
  )
}
if (profiles %in% c("dom", "both")) {
  required_directories <- c(
    required_directories, "results/dom", "figures/dom"
  )
}
invisible(lapply(
  required_directories, dir.create,
  recursive = TRUE, showWarnings = FALSE
))

scripts <- c(
  "analysis/00_install_packages.R",
  "analysis/01_build_nhanes_population.R",
  "analysis/02_validate_select.R",
  "analysis/03_run_scenarios.R",
  "analysis/04_probabilistic_sensitivity.R",
  "analysis/05_make_figures.R"
)

for (script in scripts) {
  message("Running ", script)
  source(script, local = new.env(parent = globalenv()))
}

message(
  "Pipeline complete for profile: ", profiles,
  ". See results/ and figures/."
)
