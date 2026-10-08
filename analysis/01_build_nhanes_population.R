source("data-raw/download_nhanes.R")

read_xpt_upper <- function(path) {
  x <- haven::read_xpt(path)
  names(x) <- toupper(names(x))
  x
}

demo <- read_xpt_upper("data/raw/DEMO_L.XPT")
bmx <- read_xpt_upper("data/raw/BMX_L.XPT")
diq <- read_xpt_upper("data/raw/DIQ_L.XPT")
ghb <- read_xpt_upper("data/raw/GHB_L.XPT")
mcq <- read_xpt_upper("data/raw/MCQ_L.XPT")

keep_existing <- function(data, columns) {
  data[, intersect(columns, names(data)), drop = FALSE]
}

analytic <- demo |>
  keep_existing(c(
    "SEQN", "RIDAGEYR", "RIAGENDR", "SDMVSTRA", "SDMVPSU",
    "WTPH2YR", "WTMEC2YR"
  )) |>
  dplyr::left_join(keep_existing(bmx, c("SEQN", "BMXBMI")), by = "SEQN") |>
  dplyr::left_join(keep_existing(diq, c("SEQN", "DIQ010")), by = "SEQN") |>
  dplyr::left_join(keep_existing(ghb, c("SEQN", "LBXGH", "WTPH2YR")),
                   by = "SEQN", suffix = c("", ".GHB")) |>
  dplyr::left_join(
    keep_existing(mcq, c("SEQN", "MCQ160E", "MCQ160F")), by = "SEQN"
  )

if (!"WTPH2YR" %in% names(analytic) && "WTPH2YR.GHB" %in% names(analytic)) {
  analytic$WTPH2YR <- analytic$WTPH2YR.GHB
}

required <- c(
  "SEQN", "RIDAGEYR", "SDMVSTRA", "SDMVPSU", "WTPH2YR",
  "BMXBMI", "DIQ010", "LBXGH", "MCQ160E", "MCQ160F"
)
missing_required <- setdiff(required, names(analytic))
if (length(missing_required)) {
  stop("Missing expected NHANES variables: ", paste(missing_required, collapse = ", "))
}

analytic <- analytic |>
  dplyr::mutate(
    prior_mi_or_stroke = dplyr::case_when(
      MCQ160E == 1 | MCQ160F == 1 ~ 1,
      MCQ160E == 2 & MCQ160F == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    diabetes_excluded = dplyr::case_when(
      DIQ010 == 1 | LBXGH >= 6.5 ~ 1,
      DIQ010 %in% c(2, 3) & LBXGH < 6.5 ~ 0,
      TRUE ~ NA_real_
    ),
    eligible = dplyr::case_when(
      is.na(RIDAGEYR) | is.na(BMXBMI) | is.na(prior_mi_or_stroke) |
        is.na(diabetes_excluded) ~ NA_real_,
      RIDAGEYR >= 45 & BMXBMI >= 27 & prior_mi_or_stroke == 1 &
        diabetes_excluded == 0 ~ 1,
      TRUE ~ 0
    ),
    age_group = cut(
      RIDAGEYR,
      breaks = c(45, 55, 65, 75, Inf),
      right = FALSE,
      labels = c("45-54", "55-64", "65-74", "75+")
    )
  )

design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTPH2YR,
  nest = TRUE,
  data = analytic[!is.na(analytic$eligible) & analytic$WTPH2YR > 0, ]
)

eligible_prevalence <- survey::svymean(~eligible, design, na.rm = TRUE)
eligible_total <- survey::svytotal(~eligible, design, na.rm = TRUE)

population_summary <- data.frame(
  metric = c("eligible_prevalence", "eligible_population_total"),
  estimate = c(
    unname(stats::coef(eligible_prevalence)[1]),
    unname(stats::coef(eligible_total)[1])
  ),
  standard_error = c(
    unname(survey::SE(eligible_prevalence)[1]),
    unname(survey::SE(eligible_total)[1])
  ),
  unweighted_complete_n = sum(!is.na(analytic$eligible) & analytic$WTPH2YR > 0),
  unweighted_eligible_n = sum(analytic$eligible == 1, na.rm = TRUE)
)

eligible_design <- subset(design, eligible == 1 & !is.na(age_group))
age_distribution <- survey::svymean(~age_group, eligible_design, na.rm = TRUE)
age_summary <- data.frame(
  age_group = sub("age_group", "", names(stats::coef(age_distribution))),
  proportion = unname(stats::coef(age_distribution)),
  standard_error = unname(survey::SE(age_distribution))
)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
write.csv(
  population_summary,
  "data/derived/nhanes_select_like_population.csv",
  row.names = FALSE
)
write.csv(
  age_summary,
  "data/derived/nhanes_select_like_age_distribution.csv",
  row.names = FALSE
)

