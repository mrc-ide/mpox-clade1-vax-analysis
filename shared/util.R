version_check <- function(package, version) {
  if (packageVersion(package) < version) {
    stop(sprintf(
      paste("Please upgrade %s to at least %s"),
      package, version))
  }
}

check_region <- function(region) {
  if (!(region %in% c("equateur", "sudkivu", "both"))) {
    stop("Region must be one of 'equateur', 'sudkivu' or 'both'")
  }
}

check_mixing_matrix <- function(mixing_matrix) {
  if (!(mixing_matrix %in% c("Zimbabwe", "synthetic_home", "synthetic_all"))) {
    stop("mixing_matrix must be one of 'Zimbabwe', 'synthetic_home' or
         'synthetic_all'")
  }
}

check_assumptions <- function(assumptions) {
  possible_assumptions <- c("standard")
  if (!(assumptions %in% possible_assumptions)) {
    stop(paste("assumptions must be one of:",
               paste(possible_assumptions, collapse = ", ")))
  }
}
