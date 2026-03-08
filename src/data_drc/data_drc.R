### read in and save the aggregated DRC pmcmc input data

data <- readRDS("pmcmc_inputs/data.rds")
data_linelist <- readRDS("pmcmc_inputs/data_linelist.rds")
data_who_shiny <- readRDS("pmcmc_inputs/data_who_shiny.rds")


orderly2::orderly_artefact(description = "Formatted data for fitting",
                           files = "outputs/fitting_data.rds")

orderly2::orderly_artefact(description = "Formatted data from surveillance",
                           files = "outputs/cases.rds")

orderly2::orderly_artefact(description = "Formatted data from WHO",
                           files = "outputs/who_data.rds")

dir.create("outputs", FALSE, TRUE)

saveRDS(data, "outputs/fitting_data.rds")
saveRDS(data_linelist, "outputs/cases.rds")
saveRDS(data_who_shiny, "outputs/who_data.rds")
