# pmcmc plots to assess fits run in pmcmc task

orderly_pars <- orderly2::orderly_parameters(short_run = FALSE,
                                             deterministic = FALSE,
                                             mixing_matrix = "Zimbabwe")
list2env(orderly_pars, environment())
orderly2::orderly_shared_resource("color_palette.R")
source("color_palette.R")

# inputs
orderly2::orderly_dependency(
  name = "pmcmc_burundi",
  "latest(parameter:deterministic == this:deterministic && parameter:mixing_matrix == this:mixing_matrix && parameter:short_run == this:short_run)",
  files = c("inputs/samples.rds" = "outputs/samples.rds",
            "inputs/fitting_data.rds" = "outputs/fitting_data.rds"))

# outputs
orderly2::orderly_artefact(description = "Traceplot",
                           files = c("outputs/traceplots.png"))
orderly2::orderly_artefact(description = "Rankplots",
                           files = c("outputs/rankplots.png"))
orderly2::orderly_artefact(description = "Pairs plots",
                           files = c("outputs/pairs_plots.png"))
orderly2::orderly_artefact(description = "Convergence diagnostics",
                           files = c("outputs/convergence_diagnostics.rds"))

orderly2::orderly_artefact(description = "Case fit plot",
                           files = "outputs/cases.png")

orderly2::orderly_artefact(description = "Case group fit plot",
                           files = "outputs/cases_groups.png")

orderly2::orderly_artefact(description = "Rt plot",
                           files = "outputs/Rt.png")

orderly2::orderly_artefact(description = "Posterior parameters plot",
                           files = "outputs/posteriors.png")



orderly2::orderly_resource("support.R")
source("support.R")
orderly2::orderly_shared_resource("mcmc_diagnostics.R")
source("mcmc_diagnostics.R")
orderly2::orderly_shared_resource("diagnostic_plots.R")
source("diagnostic_plots.R")
library(tidyverse)
library(GGally)
library(patchwork)
library(posterior)


# read in data dependencies
pmcmc_results <- readRDS("inputs/samples.RDS")
fitting_data <- readRDS("inputs/fitting_data.RDS")

dir.create("outputs", FALSE, TRUE)
ggsave(plot_trajectories(pmcmc_results, fitting_data, "cases"),
       filename = "outputs/cases.png", width = 15, height = 10, dpi = 600)
ggsave(plot_fitted_by_group(pmcmc_results, fitting_data, "cases"),
       filename = "outputs/cases_groups.png",
       width = 15, height = 10, dpi = 600)  
ggsave(plot_Rt(pmcmc_results, fitting_data),
       filename = paste0("outputs/Rt.png"),
       width = 7, height = 5)

ggsave(plot_posteriors(pmcmc_results, "burundi"),
       filename = "outputs/posteriors.png",
       width = 10, height = 5)
ggsave(traceplots(pmcmc_results), 
       filename = "outputs/traceplots.png",
       width = 15, height = 9, dpi = 600, bg = "white")
ggsave(rankplots(pmcmc_results), 
       filename = "outputs/rankplots.png",
       width = 15, height = 9, dpi = 600, bg = "white")
ggsave(pairs_plot(pmcmc_results),
       filename = "outputs/pairs_plots.png",
       width =  10, height = 10, dpi = 600)

ggsave(plot_prop_cases_by_age(pmcmc_results, fitting_data),
       filename = "outputs/proportion_cases_by_age.png",
       width = 10, height = 5)

convergence_diagnostics <- calc_diagnostics(pmcmc_results)
saveRDS(convergence_diagnostics, "outputs/convergence_diagnostics.rds")
