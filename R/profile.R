requested_profiles <- function() {
  requested <- tolower(Sys.getenv("PUBLICATION_PROFILE", unset = "both"))
  if (!requested %in% c("ejpc", "dom", "both")) {
    stop("PUBLICATION_PROFILE must be 'ejpc', 'dom', or 'both'.")
  }
  if (requested == "both") c("ejpc", "dom") else requested
}

profile_benchmark <- function(profile) {
  if (profile == "ejpc") "select" else "full"
}
