for (r in c("equateur", "sudkivu")) {
  orderly2::orderly_dependency(
    name  = "pmcmc_plots",
    quote(latest(parameter:region == environment:r && parameter:deterministic == FALSE && parameter:fit_by_age == TRUE && parameter:fit_KPs == TRUE && parameter:mixing_matrix == "Zimbabwe" && parameter:short_run == FALSE && parameter:assumptions == "standard")),
    files = c("outputs/fits_plots/${r}_cases.png" = "outputs/cases.png",
              "outputs/fits_plots/${r}_deaths.png" = "outputs/deaths.png",
              "outputs/fits_plots/${r}_sw_cases_and_cfr.png" = "outputs/sw_cases_and_cfr.png",
              "outputs/diagnostic_plots/${r}_pairs_plots.png" = "outputs/pairs_plots.png",
              "outputs/diagnostic_plots/${r}_rankplots.png" = "outputs/rankplots.png",
              "outputs/diagnostic_plots/${r}_traceplots.png" = "outputs/traceplots.png",
              "inputs/${r}_convergence_diagnostics.rds" = "outputs/convergence_diagnostics.rds"))
}


orderly2::orderly_dependency(
  name = "pmcmc_plots_burundi",
  "latest(parameter:deterministic == FALSE && parameter:mixing_matrix == 'Zimbabwe' && parameter:short_run == FALSE)",
  files = c("outputs/fits_plots/bujumbura_cases.png" = "outputs/cases.png",
            "outputs/fits_plots/bujumbura_proportion_cases_by_age.png" = "outputs/proportion_cases_by_age.png",
            "outputs/diagnostic_plots/bujumbura_pairs_plots.png" = "outputs/pairs_plots.png",
            "outputs/diagnostic_plots/bujumbura_rankplots.png" = "outputs/rankplots.png",
            "outputs/diagnostic_plots/bujumbura_traceplots.png" = "outputs/traceplots.png",
            "inputs/bujumbura_convergence_diagnostics.rds" = "outputs/convergence_diagnostics.rds"))

library(gt)
library(tidyverse)

orderly2::orderly_resource("diagnostics.Rmd")
orderly2::orderly_artefact(description = "diagnostics",
                           files = "outputs/diagnostics.docx")

rmarkdown::render("diagnostics.Rmd", output_file = "outputs/diagnostics.docx")
