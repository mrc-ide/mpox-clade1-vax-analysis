#----------------------------------------------#
# Run full vaccination pipeline FOR BURUNDI
#----------------------------------------------#

#----------####
# set up
#----------####
library(ggplot2)
library(tidyverse)
library(hipercow)


requiredPackages_CRAN <- c("posterior","bayesplot", "Hmisc")
for (package in requiredPackages_CRAN) { 
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(package)
}
requiredPackages_ide <- c("monty", "dust2", "odin2", "orderly2")
for (package in requiredPackages_ide) { 
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(
      package,
      repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
}
requiredPackages_GH <- c("mpoxseir")
for (package in requiredPackages_GH) { 
  if (!requireNamespace(package, quietly = TRUE))
    devtools::install_github(paste0("mrc-ide/", package))
}


#---------------------####
# Step 1 : pmcmc
#---------------------####

# set the parameters
fit_by_age <- TRUE
fit_KPs <- TRUE
region <- "bujumbura" 
use_both_fit <- FALSE
deterministic <- FALSE
short_run <- FALSE
mixing_matrix <- "Zimbabwe"
vaccines_onset <- "start"


# run locally ... full run
orderly2::orderly_run("pmcmc_burundi",
                      parameters = list(mixing_matrix = mixing_matrix,
                                        deterministic = deterministic,
                                        short_run = FALSE))

#----------------------####
# Step 2: pmcmc plots 
#----------------------####

# check you have an output from step 3 that matches your parameters
orderly2::orderly_search(orderly2::orderly_query(quote(parameter:short_run == FALSE &&
                                                         parameter:mixing_matrix == "Zimbabwe" &&
                                                         parameter:deterministic == FALSE),
                                                 name = "pmcmc_burundi"))
orderly2::orderly_run("pmcmc_plots_burundi", 
                      parameters = list(deterministic = deterministic,
                                        short_run = short_run,
                                        mixing_matrix = mixing_matrix))



#--------------------------------####
# Step 3: set up scenarios to run
#--------------------------------####

# check you have an output from previous steps that matches your parameters
orderly2::orderly_search(orderly2::orderly_query(quote(parameter:mixing_matrix == "Zimbabwe"  &&
                                                         parameter:short_run == FALSE &&
                                                         parameter:deterministic == FALSE),
                                                 name = "pmcmc_burundi"))

# set up scenarios to run
scenario_grid_burundi <-  read.csv("./shared/Burundi_scenario_grid.csv") 

# convert daily doses to have a minimum week long vaccine roll out
scenario_grid_burundi$doses_per_day_total<-sapply(
  1:nrow(scenario_grid_burundi),
  function(x) floor(min(scenario_grid_burundi$doses_per_day_total[x], 
                        (scenario_grid_burundi$total_doses_adults[x]+scenario_grid_burundi$total_doses_children[x])/7 )))

scenario_grid_burundi |>
  select(-region) |>
  mutate(total_doses_children=as.numeric(total_doses_children),
         total_doses_adults=as.numeric(total_doses_adults),
         doses_per_day_total=as.numeric(doses_per_day_total),
         daily_vax_split_children=as.numeric(daily_vax_split_children),
         days_between_doses=as.numeric(days_between_doses),
         prioritisation_children = as.character(prioritisation_children),
         uptake_realised=as.numeric(uptake_realised))-> scenario_grid_burundi

scenario_grid_burundi$scenario_num <- 1:nrow(scenario_grid_burundi)
scenario_grid_burundi <- scenario_grid_burundi |> mutate(scenario_new = scenario_name)

# save this scenario grid subset to the shared folder for later tasks to read in
saveRDS(scenario_grid_burundi, file = "./shared/scenario_grid_subset_burundi.RDS") 


#--------------------------------####
# Step 4: run
#--------------------------------####

# set working directory to run_scenario_burundi
scenario_plots <- scenario_grid_burundi

for (i in 1:nrow(scenario_plots)) {
orderly2::orderly_run(
  "run_scenario_burundi",
  parameters = with(scenario_plots,
                    list(deterministic = deterministic,
                    short_run = short_run,
                    mixing_matrix=mixing_matrix,
                    vaccines_onset = vaccines_onset,
                    vaccine_used = vaccine_used,
                    t_ve=t_ve[i],
                    total_doses_children=total_doses_children[i],
                    total_doses_adults=total_doses_adults[i],
                    doses_per_day_total=doses_per_day_total[i],
                    daily_vax_split_children=daily_vax_split_children[i],
                    vaccine_dose_scenario =vaccine_dose_scenario[i],
                    days_between_doses = days_between_doses[i],
                    prioritisation_children =prioritisation_children[i],
                    prioritisation_adults = prioritisation_adults[i],
                    uptake_realised = uptake_realised[i],
                    scenario_num = scenario_num[i]
                    )))
}


#-------------------------------------####
# Step 5: Plot multiple scenarios runs 
#-------------------------------------####

orderly2::orderly_run("run_postprocess_burundi",
                      parameters= list(deterministic = deterministic,
                                       vaccines_onset = vaccines_onset,
                                       short_run = short_run))

