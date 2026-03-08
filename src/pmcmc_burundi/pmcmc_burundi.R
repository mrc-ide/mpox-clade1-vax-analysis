library(monty)
library(tidyverse)
library(mpoxseir)
library(posterior)
library(bayesplot)
library(abind)
library(eigen1)

orderly2::orderly_resource("support.R")
orderly2::orderly_resource("packer.R")
source("support.R")
source("packer.R")

orderly2::orderly_shared_resource("util.R")
orderly2::orderly_shared_resource("index.R")
orderly2::orderly_shared_resource("calc-Rt.R")
source("util.R")
source("index.R")
source("calc-Rt.R")

orderly_pars <- orderly2::orderly_parameters(deterministic = TRUE,
                                             mixing_matrix = "Zimbabwe",
                                             short_run = TRUE)
list2env(orderly_pars, environment())

check_mixing_matrix(mixing_matrix)

orderly2::orderly_dependency(name = "data_burundi", "latest()",
                             files =  c("inputs/data.rds" = "outputs/pmcmc_data.rds"))

r <- "bujumbura"
for (type in c("deterministic", "stochastic")) {
  orderly2::orderly_resource(
    c(paste("parameters", r, type, "info.csv", sep = "/"),
      paste("parameters", r, type, "proposal.csv", sep = "/")))
  
}


orderly2::orderly_artefact(description = "Samples object",
                           files = "outputs/samples.rds")

orderly2::orderly_artefact(description = "Formatted data for fitting",
                           files = "outputs/fitting_data.rds")

orderly2::orderly_artefact(description = "Tuned inputs for future runs",
                           files = c("outputs/info.csv",
                                     "outputs/proposal.csv"))

version_check("mpoxseir", "0.2.28")
version_check("monty", "0.3.28")
version_check("dust2", "0.3.22")


start_date <- "2024-05-01"

## end at last complete week
end_date <- "2025-01-05"

## Save a snapshot at first time point in data for counterfactual
snapshots <- "2024-08-18"

raw_data <- readRDS("inputs/data.rds")

data <- get_data(raw_data, start_date, end_date)

mcmc_pars <- create_mcmc_pars(deterministic, region = r)

n_chains <- 4
if (short_run) {
  n_burnin <- 5
  n_steps <- 20
  n_sample <- 20
  n_particles <- 8
} else {
  n_burnin <- 5000
  n_steps <- if (deterministic) 15000 else 105000
  n_sample <- 1000
  n_particles <- 400 ## ideally a multiple of 8 for efficient parallelisation
}

if (deterministic) {
  n_particles <- 1
}


control <- fit_control(deterministic, n_steps = n_steps, n_burnin = n_burnin,
                       n_sample = n_sample, n_chains = n_chains,
                       n_particles = n_particles)

filter <- create_filter(data, deterministic, start_date, control)
packer <- create_packer(mixing_matrix, region = r)

pmcmc_results <- run_pmcmc(filter, packer, mcmc_pars, control$pmcmc,
                           deterministic, snapshots)

dir.create("outputs", FALSE, TRUE)

pmcmc_results <- thin_samples(pmcmc_results, control)
pmcmc_results <- calc_Rt_fits(pmcmc_results, r)

saveRDS(pmcmc_results, "outputs/samples.rds")
saveRDS(data, "outputs/fitting_data.rds")
save_pars_inputs(pmcmc_results, mcmc_pars, r)
