# -------------------------------------------------------------------- #
# ---- Setting up the dependencies and parameters from orderly ---- ###
# -------------------------------------------------------------------- #

orderly_pars <-
  orderly2::orderly_parameters(region = "burundi",
                               use_both_fit = FALSE,
                               mixing_matrix = "Zimbabwe",
                               deterministic = FALSE,
                               short_run = FALSE,
                               fit_by_age = TRUE,
                               fit_KPs = TRUE,
                               n_weeks=104,
                               scenario_num=1,
                               total_doses_children=1000, ## number of doses of LC16m8 (children only)
                               total_doses_adults=2000,  ## number of doses of MVA-BN (over 18 only)
                               doses_per_day_total=200, ## total number of vaccines that can be given daily (for adults and children)
                               daily_vax_split_children=0.5, ## what proportion of doses_per_day will be allocated to vaccinating children over adults
                               vaccine_dose_scenario ="2_doses_prioritise_delay", # "1st_dose_only" or "2_doses_prioritise_delay" or "2_doses_prioritise_1st" or no_vaccine ## are we giving 2nd doses or only first doses
                               days_between_doses = 28, ## for adults, how many days do we wait between giving 1st and 2nd doses (default is 28; if only in a 1st dose scenario then this parameter isn't used),
                               prioritisation_children ="children_all_equal",  # "children_all_equal", "children_age_gradient_all","children_age_younger_older")
                               prioritisation_adults ="adults_all_equal",# "adults_all_equal", "adults_all_key_pops", "adults_HCW_first_KP"
                               uptake_realised = 1,
                               t_ve = 28,
                               vaccine_used = "mix",
                               vaccines_onset = "end")
list2env(orderly_pars, environment())


# load a samples.rds output from
orderly2::orderly_dependency(name = "pmcmc_burundi",
                             "latest(parameter:mixing_matrix == this:mixing_matrix && parameter:deterministic == this:deterministic && parameter:short_run == this:short_run)",
                             files =  c("inputs/samples.rds" = "outputs/samples.rds",
                                        "inputs/fitting_data.rds" = "outputs/fitting_data.rds"))



# Declare files that you promise to produce, and describe them
orderly2::orderly_artefact(description = "Output from mpoxseir runs",
                           files =paste0("outputs/model_run_output_",
                                         scenario_num,".rds"))


# Loading needed packages and files
library(mpoxseir)
library(dust2)
library(abind)
library(dplyr)
library(eigen1)

## Fitting outputs
samples <- readRDS("inputs/samples.rds")
data <- readRDS("inputs/fitting_data.rds")

# Parameters setting function
orderly2::orderly_shared_resource("util.R")
orderly2::orderly_shared_resource("index.R")
orderly2::orderly_shared_resource("calc-Rt.R")
orderly2::orderly_shared_resource("interpolate_daily_doses.R")
orderly2::orderly_shared_resource("get_vaccine_scenario_parameters.R")
orderly2::orderly_shared_resource("vaccine_prioritisation_scenarios.xlsx")
source("util.R")
source("index.R")
source("calc-Rt.R")
source("interpolate_daily_doses.R")
source("get_vaccine_scenario_parameters.R")

version_check("mpoxseir", "0.2.28")
version_check("monty", "0.3.28")
version_check("dust2", "0.3.22")

check_mixing_matrix(mixing_matrix)


# ---------------------------------------------------------------------------
#                       SETTING UP THE PARAMETERS FOR VACCINATION
# ----------------------------------------------------------------------------

## take you sampled parameters and transform each parameter set,
## creating a list of lists
## here's where we will probably want to edit the transformed parameters in
## some way (update vaccine schedules)
pars <- lapply(seq_len(ncol(samples$pars)),
               function (i) samples$packer$unpack(samples$pars[, i]))
if (use_both_fit) {
  pars <- lapply(pars, function(x) x[[region]])
}


# ## get the end time point of the fits
t1 <- ifelse(vaccines_onset=="start", min(data$data$date), max(data$data$date)) ## this means this day is already covered so we care about running forward from t1 + 1


# Run scenario for 72 weeks to take into account delay in Bujumbura onset compared to DRC
# - scenarios will run until end of DRC scenarios 371+104*7 = 1099 from start of Bujumbura epidemic 595; (1099-595)/7 = 72 weeks
if(vaccines_onset=="start"){
  n_weeks <- 72 
}

## Get the vaccine scenario parameters
vax_pars <- get_vaccine_scenario_parameters(
  run_pars = pars[[1]], ## information needed will be consistent across parameter sets
  total_doses_children = as.integer(total_doses_children),
  total_doses_adults = as.integer(total_doses_adults),
  doses_per_day_total= as.integer(doses_per_day_total),
  daily_vax_split_children= daily_vax_split_children,
  vaccine_dose_scenario = vaccine_dose_scenario,
  days_between_doses = days_between_doses,
  prioritisation_children = prioritisation_children,
  prioritisation_adults = prioritisation_adults,
  uptake_realised=uptake_realised,
  t_ve = t_ve,
  t1=t1,
  vaccine_used = vaccine_used
)

for (i in 1:length(pars)) {
  # Update the parameters from the fittings with the vaccine scenario parameters
  for (nm in names(vax_pars$model_vax_inputs)) {
    pars[[i]][[nm]] <- vax_pars$model_vax_inputs[[nm]]
  }
}



# ---------------------------------------------------------------------------
#                               RUNNING THE MODEL
# ----------------------------------------------------------------------------

## setting up the model
## the 1 represents one "particle" per parameter set
sys <- dust2::dust_system_create(mpoxseir::model_targeted_vax, pars, time = t1,
                                 n_particles = 1, n_groups = length(pars),
                                 seed = 1, dt = 1, deterministic = deterministic,
                                 preserve_particle_dimension = FALSE)

## update the state with the end states from the fits (one per parameter set)
if(vaccines_onset == "end") {
  state <- samples$observations$state
  if (use_both_fit) {
    state <- state[, region, ]
  }
} else {
  ## Take the first snapshot (though there only is one currently!)
  if (use_both_fit) {
    state <- samples$observations$snapshots[, region, 1, ]
  } else {
    state <- samples$observations$snapshots[, 1, ]
  }
}
dust2::dust_system_set_state(sys, state)

## if we want to just save certain trajectories (like in the fits), we
## set the index to save here

full_index <- dust2::dust_unpack_index(sys)
index <- save_index(vax = TRUE, rt = TRUE)

## Running the model n_weeks forward from t1
days_per_week <- 7

## length of time set to two years from beginning of vaccination
output_times <- c(0, seq_len(n_weeks) * days_per_week) +  t1
output <- dust2::dust_system_simulate(sys, output_times, index_state = unlist(full_index[index$idx]))
rownames(output) <- index$names

output <- calc_Rt_scenarios(output, pars, region)


### Post process the outputs so it is in a nice data frame

# Get the dimensions of the array
dims <- dim(output)
# Generate a grid of indices for each element in the array
indices <- as.data.frame(do.call(expand.grid, lapply(dims, seq_len)))
# Rename the columns for clarity
colnames(indices) <- c("Category", "Particle", "TimeStep")
# Flatten the array into a vector and bind it to the indices
df <- cbind(indices, Value = as.vector(output))

# Add category names
category_names <- names(output[, 1, 1])

df$Category <- factor(df$Category, labels = category_names)

### Saving the output in the Outputs folder
dir.create("outputs", FALSE, TRUE)
saveRDS(df, paste0("outputs/model_run_output_",scenario_num,".rds"))

