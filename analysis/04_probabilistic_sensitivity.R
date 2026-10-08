source("R/parameters.R")
source("R/model_rhs.R")
source("R/simulation.R")
source("R/monte_carlo.R")
source("R/profile.R")

p <- readRDS("results/shared/calibrated_parameters.rds")
profiles <- requested_profiles()
parameter_table <- read.csv(
  "config/parameters.csv", stringsAsFactors = FALSE, check.names = FALSE
)
n_psa <- as.integer(Sys.getenv("N_PSA", unset = "1000"))
if (is.na(n_psa) || n_psa < 1) {
  stop("N_PSA must be a positive integer.")
}

psa <- run_probabilistic_sensitivity(
  n = n_psa,
  base_p = p,
  parameter_table = parameter_table
)

saveRDS(psa, "results/shared/probabilistic_sensitivity.rds")
write.csv(
  data.frame(
    n_psa = n_psa,
    seed = 20260928,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ),
  "results/shared/run_metadata.csv",
  row.names = FALSE
)

for (profile in profiles) {
  output_dir <- file.path("results", profile)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  psa_summary <- summarize_psa(psa, profile = profile)
  write.csv(
    psa_summary,
    file.path(output_dir, "probabilistic_sensitivity_summary.csv"),
    row.names = FALSE
  )
}
