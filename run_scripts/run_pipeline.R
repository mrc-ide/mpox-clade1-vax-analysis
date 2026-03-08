#----------------------------------------------#
# Run full vaccination pipeline FOR DRC
#----------------------------------------------#
#
#--------------------------#
# Steps in the pipeline 
#--------------------------#
# 1. pmcmc
# 2. pmcmc_plots 
# 3. running scenarios 
# 4. simulation_plots 

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
deterministic <- FALSE
short_run <- FALSE
fit_by_age <- TRUE
fit_KPs <- TRUE
region <- "sudkivu" 
use_both_fit <- FALSE
R0_SW_reduction <- 0
vaccines_onset <- "start"

# run locally 
orderly2::orderly_run("pmcmc",
                      parameters = list(region = region,
                                        deterministic = deterministic,
                                        short_run = short_run,
                                        fit_by_age = fit_by_age,
                                        fit_KPs = fit_KPs))


#----------------------####
# Step 2: pmcmc plots 
#----------------------####

# run locally 
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

scenario_grid<-read_csv("shared/DRC_scenario_grid.csv")

# convert daily doses to have a minimum week long vaccine roll out
scenario_grid$doses_per_day_total<-sapply(
  1:nrow(scenario_grid),
  function(x) floor(min(scenario_grid$doses_per_day_total[x], 
                        (scenario_grid$total_doses_adults[x]+scenario_grid$total_doses_children[x])/7 )))

region_fit <- if (use_both_fit) "both" else region
region_text <- region 

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

# save this scenario grid subset to the shared folder for later tasks to read in
saveRDS(scenario_plots, file = "./shared/scenario_grid_subset.RDS") 


#--------------------------------####
# Step 4: Run scenario  
#--------------------------------####

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

# plotting dynamics of one locally-run simulation outputs
output <- readRDS(paste0("./archive/run_scenario/",scenario_run_example[j],"/outputs/model_run_output_1.rds")) %>%
  group_by(Category, TimeStep) %>%
  summarise(median_value = median(Value, na.rm = TRUE),
            lower_95_value = quantile(Value, 0.025, na.rm = TRUE),
            upper_95_value = quantile(Value, 0.975, na.rm = TRUE))

dyn_to_plot <- output %>% filter(Category %in% c("cases_inc","deaths_inc"))
ggplot(dyn_to_plot) +
  geom_line(aes(x = TimeStep, y = median_value),color = "brown",linewidth=0.7) +
  geom_ribbon(aes(x = TimeStep,ymin = lower_95_value, ymax = upper_95_value), fill = "brown", alpha = 0.3) +
  facet_wrap(~Category, scales = "free_y") +
  labs(x = "Time (weeks)",
       y = "Total events" ) +
  theme_minimal(base_size = 13)

#-------------------------------------####
# Step 5: Plot multiple scenarios runs 
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

# post process for no R0 reduction
orderly2::orderly_run("run_postprocess",
                      parameters= list(region = region, 
                                       deterministic = deterministic,
                                       short_run = short_run,
                                       fit_by_age = fit_by_age,
                                       fit_KPs = fit_KPs,
                                       use_both_fit = use_both_fit,
                                       R0_SW_reduction = R0_SW_reduction,
                                       vaccines_onset = vaccines_onset))
