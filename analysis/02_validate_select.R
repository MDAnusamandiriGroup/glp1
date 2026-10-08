source("R/parameters.R")
source("R/model_rhs.R")
source("R/simulation.R")
source("R/calibration.R")

p <- calibrate_select_model(read_model_parameters())
validation <- calibration_table(p)

dir.create("results/shared", recursive = TRUE, showWarnings = FALSE)
write.csv(
  validation,
  "results/shared/select_calibration_validation.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    parameter = c(
      "lambda_mace_direct", "lambda_mace_calibrated",
      "lambda_mace_calibration_factor", "hr_mace_on_calibrated",
      "select_discontinuation_probability_1y"
    ),
    value = c(
      p$lambda_mace_direct, p$lambda_mace,
      p$lambda_mace_calibration_factor, p$hr_mace_on,
      p$select_discontinue_probability_1y
    )
  ),
  "results/shared/calibrated_parameters.csv",
  row.names = FALSE
)
saveRDS(p, "results/shared/calibrated_parameters.rds")
