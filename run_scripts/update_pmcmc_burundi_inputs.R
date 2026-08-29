update_inputs <- function(assumptions, from_deterministic, to_deterministic) {
  ## Set working directory to project directory
  ##setwd(orderly2:::orderly_src_root(NULL, TRUE))
  
  ## Find the latest report for the given region
  latest <- 
    orderly2::orderly_search(quote(latest(parameter:short_run == FALSE &&
                                            parameter:assumptions == this:assumptions &&
                                            parameter:deterministic == this:deterministic)),
                             name = "pmcmc_burundi",
                             parameters = list(deterministic = from_deterministic,
                                               assumptions = assumptions))
  
  ## Folder location in severity_parameters
  type <- if (to_deterministic) "deterministic" else "stochastic"
  pars_folder <- paste("src/pmcmc_burundi/parameters/bujumbura", type, assumptions, sep = "/")
  
  ## Copy the latest info and proposal csvs into that folder location
  orderly2::orderly_copy_files(latest,
                               files = c("proposal.csv" = "outputs/proposal.csv",
                                         "info.csv" = "outputs/info.csv"),
                               dest = pars_folder)
}

## example usage: update the inputs for Equateur from the latest determinstic
## fit to use for deterministic fitting
# update_inputs(from_deterministic = TRUE, to_deterministic = TRUE)

remove_parameter <- function(assumptions, deterministic, par_name) {
  type <- if (deterministic) "deterministic" else "stochastic"
  pars_folder <- paste("src/pmcmc_burundi/parameters/bujumbura", type, assumptions, 
                       sep = "/")
  
  info_file <- paste(pars_folder, "info.csv", sep = "/")
  info <- read.csv(info_file)
  
  if (region == "both") {
    shared <- any(info$region[info$name == par_name] == "both")
  }
  
  info <- info[info$name != par_name, ]
  write.csv(info, info_file, row.names = FALSE)
  
  proposal_file <- paste(pars_folder, "proposal.csv", sep = "/")
  proposal <- read.csv(proposal_file, row.names = 1)
  if (region == "both") {
    colnames(proposal) <- rownames(proposal)
    if (shared) {
      proposal <- proposal[rownames(proposal) != par_name, ]
      proposal <- proposal[, colnames(proposal) != par_name]
    } else {
      proposal <- 
        proposal[!grepl(paste0("^", par_name, "<"), rownames(proposal)), ]
      proposal <- 
        proposal[, !grepl(paste0("^", par_name, "<"), colnames(proposal))]
    }
  } else {
    proposal <- proposal[rownames(proposal) != par_name, ]
    proposal <- proposal[, colnames(proposal) != par_name]
  }
  
  write.csv(proposal, proposal_file)
}
