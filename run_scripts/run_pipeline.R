#----------------------------------------------#
# Run full vaccination pipeline 
#----------------------------------------------#
#
#--------------------------#
# Steps in the pipeline 
#--------------------------#

# 1. pmcmc
# 2. pmcmc_plots 
# 3. running scenarios 
# 4. running post process for scenarios

#----------####
# set up
#----------####
beeps <- TRUE
if (beeps) library(beepr)

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
deterministic <- FALSE
short_run <- FALSE
fit_by_age <- TRUE
fit_KPs <- TRUE
region <- "sudkivu" # sudkivu/equateur/both
use_both_fit <- FALSE
R0_SW_reduction <- 0
vaccines_onset <- "start"
assumptions <-"fix_prop_SW"

# run locally ...
# Note this is not recommended and should be done on a HPCU if available.
orderly2::orderly_run("pmcmc",
                      parameters = list(region = region,
                                        deterministic = deterministic,
                                        short_run = short_run,
                                        fit_by_age = fit_by_age,
                                        fit_KPs = fit_KPs,
                                        assumptions = assumptions))

# ... or pull an example from share drive 
# note - check the parameters for what you are pulling
orderly2::orderly_location_pull(
  quote(parameter:region == 'sudkivu' &&
          parameter:short_run == FALSE &&
          parameter:fit_by_age == TRUE &&
          parameter:deterministic == FALSE &&
          parameter:fit_KPs == TRUE &&
          parameter:assumptions == "fix_prop_SW"),
  name = "pmcmc")

orderly2::orderly_location_pull(
  quote(parameter:region == 'equateur' &&
          parameter:short_run == FALSE &&
          parameter:fit_by_age == TRUE &&
          parameter:deterministic == FALSE &&
          parameter:fit_KPs == TRUE &&
          parameter:assumptions == "fix_prop_SW"),
  name = "pmcmc")

orderly2::orderly_location_pull(
  quote(parameter:region == 'both' &&
          parameter:short_run == FALSE &&
          parameter:fit_by_age == TRUE &&
          parameter:deterministic == FALSE &&
          parameter:fit_KPs == TRUE &&
          assumptions == "fix_prop_SW"),
  name = "pmcmc")


#----------------------####
# Step 2: pmcmc plots 
#----------------------####

# run locally ...
# check you have an output from step 3 that matches your parameters
orderly2::orderly_search(orderly2::orderly_query(quote(parameter:region == "sudkivu" &&
                                                         parameter:short_run == FALSE &&
                                                         parameter:fit_by_age == TRUE &&
                                                         parameter:deterministic == FALSE &&
                                                         parameter:fit_KPs == TRUE),
                                                 name = "pmcmc"))
orderly2::orderly_run("pmcmc_plots", 
                      parameters = list(region = region, 
                                        deterministic = deterministic,
                                        short_run = short_run,
                                        fit_by_age = fit_by_age,
                                        fit_KPs = fit_KPs))

# ... or pull example from share drive
# note - check the parameters for what you are pulling
orderly2::orderly_location_pull(
  quote(parameter:region == 'sudkivu' &&
          parameter:short_run == FALSE &&
          parameter:fit_by_age == TRUE &&
          parameter:deterministic == FALSE &&
          parameter:fit_KPs == TRUE &&
          assumptions == assumptions),
  name = "pmcmc_plots") 

orderly2::orderly_location_pull(
  quote(parameter:region == 'equateur' &&
          parameter:short_run == FALSE &&
          parameter:fit_by_age == TRUE &&
          parameter:deterministic == FALSE &&
          parameter:fit_KPs == TRUE &&
          assumptions == assumptions),
  name = "pmcmc_plots") 

orderly2::orderly_location_pull(
  quote(parameter:region == 'both' &&
          parameter:short_run == FALSE &&
          parameter:fit_by_age == TRUE &&
          parameter:deterministic == FALSE &&
          parameter:fit_KPs == TRUE &&
          assumptions == assumptions),
  name = "pmcmc_plots") 


#--------------------------------####
# Step 3: Run scenarios 
#--------------------------------####

# check you have an output from previous steps that matches your parameters
orderly2::orderly_search(orderly2::orderly_query(quote(parameter:region == "equateur" &&
                                                         parameter:short_run == FALSE &&
                                                         parameter:fit_by_age == TRUE &&
                                                         parameter:deterministic == FALSE &&
                                                         parameter:fit_KPs == TRUE),
                                                 name = "pmcmc"))

#scenario_grid <- readRDS(file = "./shared/scenario_grid.RDS") 
scenario_grid<-read_csv("shared/DRC_dailydose.csv")

# Convert daily doses to have a minimum week long vaccine roll out
scenario_grid<-scenario_grid|>mutate(scenario_doses_per_day_total=doses_per_day_total)
scenario_grid$doses_per_day_total<-sapply(
  1:nrow(scenario_grid),
  function(x) floor(min(scenario_grid$doses_per_day_total[x], 
                        (scenario_grid$total_doses_adults[x]+scenario_grid$total_doses_children[x])/7 )))

region_fit <- if (use_both_fit) "both" else region
region_text <- region #if (use_both_fit) c("equateur", "sudkivu") else region

scenario_grid |>
  filter(region %in% region_text) |>
  select(-region) |>
  mutate(total_doses_children=as.numeric(total_doses_children),
         total_doses_adults=as.numeric(total_doses_adults),
         doses_per_day_total=as.numeric(doses_per_day_total),
         daily_vax_split_children=as.numeric(daily_vax_split_children),
         days_between_doses=as.numeric(days_between_doses),
         prioritisation_children = as.character(prioritisation_children),
         prioritisation_adults = as.character(prioritisation_adults),
         vaccine_dose_scenario = as.character(vaccine_dose_scenario),
         uptake_realised=as.numeric(uptake_realised))-> scenario_grid
scenario_grid$scenario_num <- 1:nrow(scenario_grid)
scenario_plots<- scenario_grid|>mutate(scenario_new=scenario_name)

# save this scenario grid subset to the share drive for later tasks to read in
saveRDS(scenario_plots, file = "./shared/scenario_grid_subset.RDS") 


# Run scenario grid locally
for (j in 1:nrow(scenario_plots)) {
  orderly2::orderly_run(
    "run_scenario",
    parameters = with(scenario_plots,
                      list(region = region,
                           deterministic = deterministic,
                           short_run = short_run,
                           fit_by_age = fit_by_age,
                           use_both_fit = use_both_fit,
                           vaccine_used = vaccine_used[j],
                           t_ve=t_ve[j],
                           R0_SW_reduction = R0_SW_reduction,
                           vaccines_onset = vaccines_onset,
                           total_doses_children=total_doses_children[j],
                           total_doses_adults=total_doses_adults[j],
                           doses_per_day_total=doses_per_day_total[j],
                           daily_vax_split_children=daily_vax_split_children[j],
                           vaccine_dose_scenario =vaccine_dose_scenario[j],
                           days_between_doses = days_between_doses[j],
                           prioritisation_children =prioritisation_children[j],
                           prioritisation_adults = prioritisation_adults[j],
                           uptake_realised = uptake_realised[j],
                           scenario_num = scenario_num[j]
                      )))
}
if (beeps) beep()

#-------------------------------------####
# Step 5: Run post process for plots 
#------------------------------------#####

# check you have an output from previous steps that matches your parameters
orderly2::orderly_search(orderly2::orderly_query(quote(parameter:region == "equateur" &&
                                                         parameter:short_run == FALSE &&
                                                         parameter:fit_by_age == TRUE &&
                                                         parameter:deterministic == FALSE &&
                                                         parameter:fit_KPs == TRUE &&
                                                         parameter:R0_SW_reduction == 0 &&
                                                         parameter:vaccines_onset == "end"&&
                                                         parameter:use_both_fit == FALSE),
                                                 name = "run_scenario"))

# Post process all the orderly outputs run on the cluster to create one big data frame with outputs of interest 

# run locally .. 
# Post process for no R0 reduction
orderly2::orderly_run("run_postprocess",
                      parameters= list(region = region, 
                                       deterministic = deterministic,
                                       short_run = short_run,
                                       fit_by_age = fit_by_age,
                                       fit_KPs = fit_KPs,
                                       use_both_fit = use_both_fit,
                                       R0_SW_reduction = R0_SW_reduction,
                                       vaccines_onset = vaccines_onset))


