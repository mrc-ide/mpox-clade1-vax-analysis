update_inputs <- function(from_deterministic, to_deterministic) {
  ## Set working directory to project directory
  setwd(orderly2:::orderly_src_root(NULL, TRUE))
  
  ## Find the latest report for the given region
  latest <- 
    orderly2::orderly_search(quote(latest(parameter:short_run == FALSE &&
                                            parameter:deterministic == this:deterministic)),
                             name = "pmcmc_burundi",
                             parameters = list(deterministic = from_deterministic))
  
  ## Folder location in severity_parameters
  type <- if (to_deterministic) "deterministic" else "stochastic"
  pars_folder <- paste("src/pmcmc_burundi/parameters/bujumbura", type, sep = "/")
  
  ## Copy the latest info and proposal csvs into that folder location
  orderly2::orderly_copy_files(latest,
                               files = c("proposal.csv" = "outputs/proposal.csv",
                                         "info.csv" = "outputs/info.csv"),
                               dest = pars_folder)
}

## example usage: update the inputs for Equateur from the latest determinstic
## fit to use for deterministic fitting
