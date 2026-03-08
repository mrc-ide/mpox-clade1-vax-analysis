# pmcmc plots to assess fits run in pmcmc task

orderly_pars <- orderly2::orderly_parameters(region = 'sudkivu',
                                             mixing_matrix = 'Zimbabwe',
                                             short_run = FALSE,
                                             fit_by_age = TRUE,
                                             deterministic = FALSE,
                                             fit_KPs = TRUE,
                                             assumptions = "standard")
list2env(orderly_pars, environment())
orderly2::orderly_shared_resource("color_palette.R")
source("color_palette.R")

# inputs
orderly2::orderly_dependency(
  name = "pmcmc",
  "latest(parameter:region == this:region && parameter:deterministic == this:deterministic && parameter:fit_by_age == this:fit_by_age && parameter:fit_KPs == this:fit_KPs && parameter:mixing_matrix == this:mixing_matrix && parameter:short_run == this:short_run  && parameter:assumptions == this:assumptions)",
  files = c("inputs/samples.rds" = "outputs/samples.rds",
            "inputs/fitting_data.rds" = "outputs/fitting_data.rds"))

# outputs
orderly2::orderly_artefact(description = "Traceplot",
                           files = c("outputs/traceplots.png"))
orderly2::orderly_artefact(description = "Pairs plots",
                           files = c("outputs/pairs_plots.png"))
orderly2::orderly_artefact(description = "Rankplots",
                           files = c("outputs/rankplots.png"))
orderly2::orderly_artefact(description = "Convergence diagnostics",
                           files = c("outputs/convergence_diagnostics.rds"))

if (region == "both") {
  orderly2::orderly_artefact(description = "Case fit plot",
                             files = c("outputs/cases_equateur.png",
                                       "outputs/cases_sudkivu.png"))
  
  orderly2::orderly_artefact(description = "Death fit plot",
                             files = c("outputs/deaths_equateur.png",
                                       "outputs/deaths_sudkivu.png"))
  
  orderly2::orderly_artefact(description = "SW cases and CFR fit plot",
                             files = c("outputs/sw_cases_and_cfr.png",
                                       "outputs/sw_cases_and_cfr.png"))
  
  orderly2::orderly_artefact(description = "Case groups fit plot",
                             files = c("outputs/cases_groups_equateur.png",
                                       "outputs/cases_groups_sudkivu.png"))
  
  orderly2::orderly_artefact(description = "Death groups fit plot",
                             files = c("outputs/deaths_groups_equateur.png",
                                       "outputs/deaths_groups_sudkivu.png"))
  
  orderly2::orderly_artefact(description = "Rt plot",
                             files = c("outputs/Rt_equateur.png",
                                       "outputs/Rt_sudkivu.png"))
} else {
  orderly2::orderly_artefact(description = "Case fit plot",
                             files = "outputs/cases.png")
  
  orderly2::orderly_artefact(description = "Death fit plot",
                             files = "outputs/deaths.png")
  
  orderly2::orderly_artefact(description = "SW cases and CFR fit plot",
                             files = "outputs/sw_cases_and_cfr.png")
  
  orderly2::orderly_artefact(description = "Case group fit plot",
                             files = "outputs/cases_groups.png")
  
  orderly2::orderly_artefact(description = "Death group fit plot",
                             files = "outputs/deaths_groups.png")
  
  orderly2::orderly_artefact(description = "Rt plot",
                             files = "outputs/Rt.png")

  orderly2::orderly_artefact(description = "Posterior parameters plot",
                             files = "outputs/posteriors.png")
}



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
if (region == "both") {
  for (r in c("equateur", "sudkivu")) {
    ggsave(plot_trajectories(pmcmc_results, fitting_data, "cases", r),
           filename = paste0("outputs/cases_", r, ".png"),
           width = 15, height = 10, dpi = 600)
    ggsave(plot_trajectories(pmcmc_results, fitting_data, "deaths", r),
           filename = paste0("outputs/deaths_", r, ".png"),
           width = 15, height = 10, dpi = 600)
    ggsave(plot_fitted_by_group(pmcmc_results,fitting_data, "cases", r),
           filename = paste0("outputs/cases_groups_", r, ".png"),
           width = 15, height = 10)
    ggsave(plot_fitted_by_group(pmcmc_results,fitting_data, "deaths", r),
           filename = paste0("outputs/deaths_groups_", r, ".png"),
           width = 15, height = 10)
    ggsave(plot_Rt(pmcmc_results, fitting_data, r),
           filename = paste0("outputs/Rt_", r, ".png"),
           width = 7, height = 5)
    ggsave(plot_SW_cases_and_cfr(pmcmc_results, fitting_data, r),
           filename = paste0("outputs/sw_cases_and_cfr_", r, ".png"),
           width = 15, height = 5, dpi = 600)
  }
} else {
  ggsave(plot_trajectories(pmcmc_results, fitting_data, "cases"),
         filename = "outputs/cases.png", width = 15, height = 10, dpi = 600)
  ggsave(plot_trajectories(pmcmc_results, fitting_data, "deaths"),
         filename = "outputs/deaths.png", width = 15, height = 10, dpi = 600)
  ggsave(plot_fitted_by_group(pmcmc_results, fitting_data, "cases"),
         filename = "outputs/cases_groups.png", width = 15, height = 10)  
  ggsave(plot_fitted_by_group(pmcmc_results,fitting_data, "deaths"),
         filename = "outputs/deaths_groups.png", width = 15, height = 10)
  ggsave(plot_Rt(pmcmc_results, fitting_data),
         filename = paste0("outputs/Rt.png"),
         width = 7, height = 5)
  ggsave(plot_SW_cases_and_cfr(pmcmc_results, fitting_data),
         filename = "outputs/sw_cases_and_cfr.png",
         width = 15, height = 5, dpi = 600)

}
ggsave(plot_posteriors(pmcmc_results, region),
       filename = "outputs/posteriors.png",
       width = 10, height = 5)
ggsave(traceplots(pmcmc_results), 
       filename = "outputs/traceplots.png",
       width = 15, height = 9, dpi = 600, bg = "white")
ggsave(rankplots(pmcmc_results), 
       filename = "outputs/rankplots.png",
       width = 15, height = 9, dpi = 600,  bg = "white")
ggsave(pairs_plot(pmcmc_results),
       filename = "outputs/pairs_plots.png",
       width = 10, height = 10, dpi = 600)


convergence_diagnostics <- calc_diagnostics(pmcmc_results)
saveRDS(convergence_diagnostics, "outputs/convergence_diagnostics.rds")
