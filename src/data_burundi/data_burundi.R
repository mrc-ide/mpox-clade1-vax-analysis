### read in and save the aggregated Burundi pmcmc input data

data <- readRDS("pmcmc_inputs_burundi/data.rds")


orderly2::orderly_artefact(description = "Formatted data for fitting",
                           files = "outputs/pmcmc_data.rds")


dir.create("outputs", FALSE, TRUE)

saveRDS(data, "outputs/pmcmc_data.rds")
