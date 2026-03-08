library(dplyr)
library(ggplot2)
library(tidyr)

orderly_pars <- orderly2::orderly_parameters(deterministic = FALSE,
                                             short_run = FALSE,
                                             fit_by_age = TRUE,
                                             fit_KPs=TRUE,
                                             vaccines_onset= "end")

# Declare files that you promise to produce, and describe them
orderly2::orderly_artefact(description = "Big df with all scenarios outputs",
                           files ="outputs/runs_postprocess.RDS")

orderly2::orderly_artefact(description = "inputs folder",
                           files = "inputs/")

# read in scenario grid 
orderly2::orderly_shared_resource(scenario_grid_subset_burundi.RDS = "scenario_grid_subset_burundi.RDS") 
scenario_grid_burundi <- readRDS("scenario_grid_subset_burundi.RDS")

scenario_grid_burundi |>
  mutate(total_doses_children=as.numeric(total_doses_children),
         total_doses_adults=as.numeric(total_doses_adults),
         doses_per_day_total=as.numeric(doses_per_day_total),
         daily_vax_split_children=as.numeric(daily_vax_split_children),
         days_between_doses=as.numeric(days_between_doses),
         prioritisation_children = as.character(prioritisation_children),
         uptake_realised=as.numeric(uptake_realised))-> scenario_grid_burundi

scenario_grid_burundi$scenario_num <- 1:nrow(scenario_grid_burundi)
scenario_grid_burundi <- scenario_grid_burundi |> mutate(scenario_new = scenario_name)



# Importing all the outputs from the different scenarios 
for (i in 1:nrow(scenario_grid_burundi)) {
  
  # Extracting all parameters from each scenario  
  n_weeks <- 72
  total_doses_children <- scenario_grid_burundi$total_doses_children[i]
  total_doses_adults<- scenario_grid_burundi$total_doses_adults[i]
  doses_per_day_total <- scenario_grid_burundi$doses_per_day_total[i] 
  daily_vax_split_children <- scenario_grid_burundi$daily_vax_split_children[i] 
  vaccine_dose_scenario <- scenario_grid_burundi$vaccine_dose_scenario[i]
  days_between_doses <- scenario_grid_burundi$days_between_doses[i]
  prioritisation_children <- scenario_grid_burundi$prioritisation_children[i]
  prioritisation_adults <- scenario_grid_burundi$prioritisation_adults[i]
  uptake_realised <- scenario_grid_burundi$uptake_realised[i]
  scenario_num <- scenario_grid_burundi$scenario_num[i]            
  
  file_name <- paste0("model_run_output_",scenario_num,".rds")
  
  orderly2::orderly_dependency(
    "run_scenario_burundi",
    quote(latest( parameter:short_run == this:short_run &&
                     parameter:fit_by_age == this:fit_by_age &&
                     parameter:deterministic == this:deterministic &&
                     parameter:fit_KPs== this:fit_KPs &&
                     parameter:total_doses_adults == environment:total_doses_adults && 
                     parameter:total_doses_children == environment:total_doses_children && 
                     parameter:doses_per_day_total == environment:doses_per_day_total &&
                     parameter:daily_vax_split_children == environment:daily_vax_split_children &&
                     parameter:vaccine_dose_scenario == environment:vaccine_dose_scenario &&
                     parameter:days_between_doses == environment:days_between_doses &&
                     parameter:prioritisation_children == environment:prioritisation_children &&
                     parameter:prioritisation_adults == environment:prioritisation_adults &&
                     parameter:uptake_realised == environment:uptake_realised &&
                     parameter:vaccines_onset == environment:vaccines_onset &&
                     parameter:scenario_num == environment:scenario_num)),
    files=c("inputs/${file_name}" = paste0("outputs/model_run_output_",scenario_num,".rds")))
  
}


# -------------------------------------------------------------------- ####
# ----               Output data import                           ---- ####
# -------------------------------------------------------------------- ####


big_df <- data.frame(matrix(0, ncol=24, nrow=0))

# Extracting the counterfactual scenario
# Estimating averted and cumulative values 
# Estimating median and quantiles of the dynamics outputs and creating big 
# data frame with all the results 

counter_scenario <- scenario_grid_burundi %>% filter(vaccine_dose_scenario=="no_vaccine")
counter_id <- counter_scenario$scenario_num

counter_output <- readRDS(paste0("inputs/model_run_output_",counter_id,".rds")) %>%
  rename(counter =Value)

scenarios <- scenario_grid_burundi  %>% filter(vaccine_dose_scenario!="no_vaccine") %>% 
  mutate(scenario_name = paste0(total_doses_children, "_", 
                                total_doses_adults, "_",
                                doses_per_day_total,"_", 
                                vaccine_dose_scenario, "_", 
                                prioritisation_children, "_",prioritisation_adults  ))

for (i in 1:nrow(scenarios)){
  
  # Extracting scenario information and simulation output 
  scenario_id <- scenarios$scenario_num[i]
  output <- readRDS(paste0("inputs/model_run_output_",scenario_id,".rds")) %>% 
    left_join(counter_output, by =c("Category","Particle","TimeStep")) %>%
    mutate(averted= counter-Value)  # Calculating event averted 
  
  ## Estimating the cumulative values: 
  
  # List of categories to calculate cumulative values for
  categories_to_sum <- c("cases_inc", "deaths_inc",
                         "cases_inc_00_04", "cases_inc_05_14", "cases_inc_15_plus",
                         "cases_inc_PBS", "cases_inc_SW","cases_inc_ASW","cases_inc_CSW", "cases_inc_HCW",           
                         "deaths_inc_00_04", "deaths_inc_05_14", "deaths_inc_15_plus",
                         "deaths_inc_PBS","deaths_inc_SW", "deaths_inc_ASW","deaths_inc_CSW", "deaths_inc_HCW")
  
  # Calculate cumulative values over TimeStep for specified categories 
  # rename columns so df can be combined with original output dataframe
  output_cumulative <- output %>%
    filter(Category %in% categories_to_sum) %>%  # Filter for specified categories
    arrange(Particle, Category, TimeStep) %>%    # Arrange by Particle, Category, and TimeStep
    group_by(Particle, Category) %>%             # Group by Particle and Category
    mutate(cumulative_value = cumsum(Value),          # Cumulative sum of Value over TimeStep
           cumulative_counter = cumsum(counter),      # Cumulative sum of counter over TimeStep
           cumulative_averted = cumsum(averted)) %>%    # Cumulative sum of averted over TimeStep
    select ("Category","Particle","TimeStep", "cumulative_value", "cumulative_counter" ,"cumulative_averted")  %>%
    mutate(Category = paste0(Category, "_cum")) %>%        # Append '_cum' to all Category values
    rename (Value = cumulative_value,
            counter = cumulative_counter,
            averted =cumulative_averted) %>%
    ungroup()  #Ungroup after calculation
  
  # Big Data frame with cumulative incidence attached 
  output <- rbind(output,output_cumulative)
  
  Total_doses <-output |>
    filter(Category %in% c("dose1_cumulative", "dose2_cumulative")) |>
    group_by(Particle, TimeStep) |>
    summarise(
      cumulative_dose = sum(Value),
      .groups = "drop")|>
    full_join(output, by=c("TimeStep","Particle"))
  scenario_param <- scenarios %>% filter(scenario_num==scenario_id)
  
  # Estimating median and 95% quantiles for all over the different particles
  summary_output <- Total_doses %>% 
    group_by(Category, TimeStep) %>%
    summarise(
      # Metrics for infections  
      median_value = median(Value, na.rm = TRUE),
      mean_value = mean(Value, na.rm = TRUE),
      lower_95_value = quantile(Value, 0.025, na.rm = TRUE),
      upper_95_value = quantile(Value, 0.975, na.rm = TRUE),
      # Metrics for counter infections (no vaccine)
      median_counter = median(counter, na.rm = TRUE),
      mean_counter = mean(counter, na.rm = TRUE),
      lower_95_counter= quantile(counter, 0.025, na.rm = TRUE),
      upper_95_counter = quantile(counter, 0.975, na.rm = TRUE),
      # Metrics for counter infections averted
      median_averted = median(averted, na.rm = TRUE),
      mean_averted = mean(averted, na.rm = TRUE),
      lower_95_averted = quantile(averted, 0.025, na.rm = TRUE),
      upper_95_averted = quantile(averted, 0.975, na.rm = TRUE),
      # Metrics for percentage of infections averted
      median_percent = median(averted/counter, na.rm = TRUE),
      mean_percent = mean(averted/counter, na.rm = TRUE),
      lower_95_percent= quantile(averted/counter, 0.025, na.rm = TRUE),
      upper_95_percent = quantile(averted/counter, 0.975, na.rm = TRUE),
      # Metrics for outcomes averted per dose
      median_averted_per_dose = median(averted/cumulative_dose, na.rm = TRUE),
      mean_averted_per_dose = mean(averted/cumulative_dose, na.rm = TRUE),
      lower_averted_95_per_dose= quantile(averted/cumulative_dose, 0.025, na.rm = TRUE),
      upper_averted_95_per_dose = quantile(averted/cumulative_dose, 0.975, na.rm = TRUE),
      # Metrics for total dose
      median_cumulative_doses = median(cumulative_dose, na.rm = TRUE),
      mean_cumulative_doses = mean(cumulative_dose, na.rm = TRUE),
      lower_95_cumulative_doses = quantile(cumulative_dose, 0.025, na.rm = TRUE),
      upper_95_cumulative_doses = quantile(cumulative_dose, 0.975, na.rm = TRUE)) %>%
    mutate(scenario_num = scenario_id) %>%
    left_join(scenario_param, by ="scenario_num")
  
  # Creating big data frame with all results 
  big_df <- rbind(big_df,summary_output)
  
}

# if projecting from the start of the data use fits as the baseline up to the end of the data
if(vaccines_onset=="start"){
  mean_ci <- function(x, alpha = 0.05) {
    p <- c(alpha / 2, 1 - alpha / 2)
    q <- quantile(x, c(0.025, 0.975))
    names(q) <- paste0("q", p * 100)

    c(mean = mean(x), q)
  }
  ## load a samples.rds output from
  orderly2::orderly_dependency(name = "pmcmc_burundi",
                               "latest(parameter:deterministic == this:deterministic && parameter:short_run == this:short_run)",
                               files =  c("inputs/samples.rds" = "outputs/samples.rds",
                                          "inputs/fitting_data.rds" = "outputs/fitting_data.rds"))



  ## Fitting outputs
  samples <- readRDS("inputs/samples.rds")
  data <- readRDS("inputs/fitting_data.rds")
  region <- "burundi"
  trajectories <- apply(samples$observations$trajectories, c(1, 2), mean_ci)
  grid <- dimnames(trajectories)
  grid <- append(grid, values = region, after = 2)


  names(grid) <- c("name", "state", "region", "date")
  grid$date <- mpoxseir::mpoxseir_date_as_date(sort(unique(data$data$date)))
  fits <- expand.grid(grid)
  fits$value <- c(trajectories)
  fits_wide<-fits%>%pivot_wider(names_from = name, values_from = value)|>
    rename(Category=state, mean_counter=mean, lower_95_counter=q2.5, upper_95_counter=q97.5)

  big_df<-full_join(big_df%>%mutate(date=mpoxseir::mpoxseir_date_as_date(595+(TimeStep-1)*7)), # start weekly dates from 595
                    fits_wide%>%select(-region), by = c("Category", "date")) %>%  # join data frames by date
    mutate(mean_counter.y = if_else(is.na(mean_counter.y), mean_counter.x, mean_counter.y),
           lower_95_counter.y = if_else(is.na(lower_95_counter.y), lower_95_counter.x, lower_95_counter.y),
           upper_95_counter.y = if_else(is.na(upper_95_counter.y), upper_95_counter.x, upper_95_counter.y))%>%
    mutate(mean_counter=mean_counter.y, lower_95_counter=lower_95_counter.y, upper_95_counter=upper_95_counter.y)%>%
    select(-mean_counter.y,-mean_counter.x,-lower_95_counter.y,-lower_95_counter.x,-upper_95_counter.y,-upper_95_counter.x)

}


### Saving the output in the Outputs folder
dir.create("outputs", FALSE, TRUE)
saveRDS(big_df,"outputs/runs_postprocess.RDS")



