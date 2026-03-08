
## have to call the interpolate_daily_doses function from the shared support file

get_vaccine_scenario_parameters <- function(
    run_pars, ## imagining this as the object which was input into the pmcmc fit that we are taking as the basis for the scenario run (e.g. this will tell us things like n_vax, children_idx_raw, adult_idx_raw, location)
    ### now have a series of parameters which are not explicitly model parameters but are used to help us derive our inputs
    total_doses_children, ## number of doses of LC16m8 (children only)
    total_doses_adults, ## number of doses of MVA-BN (over 18 only)
    doses_per_day_total, ## total number of vaccines that can be given daily (for adults and children)
    daily_vax_split_children, ## what proportion of doses_per_day will be allocated to vaccinating children over adults
    vaccine_dose_scenario, ## are we giving 2nd doses or only first doses or no vaccines
    days_between_doses = 28, ## for adults, how many days do we wait between giving 1st and 2nd doses (default is 28; if only in a 1st dose scenario then this parameter isn't used),
    ## these now are formal model parameters, we use the scenario name to read in from the relevant excel tab
    ## this will depend on region as the hesitancy comes from different regions
    prioritisation_children,
    prioritisation_adults,
    uptake_realised = 1,## what proportion of those indicated to take the vaccine will take the vaccine they're offered, default to 1 (just for sens analysis if needed)
    vaccine_used = "mix", # which vaccine to use 'Lc16m8'/'MVA-BN' (or a 'mix' of both)
    t_ve = 28, ## time to reach peak vaccine efficacy (https://www.gov.uk/government/publications/monkeypox-vaccination-resources/protecting-you-from-monkeypox-information-on-the-smallpox-vaccination)
    R0_SW_reduction = 0, # reduce the R0 in SW by given percent
    vaccines_onset="end", # start projections at start or end of data 
    t1 ## last time of the model
){
  
  ### check the inputs are valid
  if(total_doses_children%%1!=0 | total_doses_adults%%1!=0 | doses_per_day_total%%1!=0 | days_between_doses%%1!=0){
    stop("parameters related to doses (total_doses_children; total_doses_adults; doses_per_day_total; days_between_doses) must be integers")
  }
  
  if(daily_vax_split_children<0|daily_vax_split_children>1){
    stop("daily_vax_split_children must be between 0 and 1")
  }
  
  if(!(vaccine_dose_scenario %in% c("1st_dose_only","2_doses_prioritise_1st",
                                    "2_doses_prioritise_delay","no_vaccine"))){
    stop("vaccine_dose_scenario must be either 1st_dose_only, 2_doses_prioritise_1st, 2_doses_prioritise_delay or no_vaccine")
  }
  
  if(!(prioritisation_children %in% c(
    "children_all_equal",
    "children_age_gradient_all"))){
    stop("prioritisation_children is not valid")
  }
  
  if(!(prioritisation_adults %in% c(
    "adults_all_equal",
    "adults_key_pops",
    "adults_key_pops_incl_PBS",
    "adults_HCW_first_KP",
    "adults_SW_first_KP", 
    "adults_children_first"))){
    stop("prioritisation_adults is not valid")
  }
  
  if(uptake_realised>1|uptake_realised<0){
    stop("uptake_realised must be between 0 and 1")
  }
  
  if((t_ve%%1!=0)){
    stop("must provide integer values for t_ve")
  }
  
  if(total_doses_children==0&total_doses_adults==0&vaccine_dose_scenario!="no_vaccine"){
    stop("if no children or adult doses then vaccine_dose_scenario must be no_vaccine")
  }
  
  idx <- mpoxseir::get_compartment_indices()
  
  doses_per_day_children <- 
    ceiling(daily_vax_split_children * doses_per_day_total)
  doses_per_day_adults <- doses_per_day_total - doses_per_day_children
  
  ## no vaccines
  if (vaccine_dose_scenario=="no_vaccine"){
    
    total_1st_doses_children <- 0
    total_1st_doses_adults <- 0
    total_2nd_doses_adults <- 0
    adults_prioritise_delay <- FALSE ## should have no effect
    
  } else if (vaccine_dose_scenario == "1st_dose_only") {
    
    total_1st_doses_children <- total_doses_children
    total_1st_doses_adults <- total_doses_adults
    total_2nd_doses_adults <- 0
    adults_prioritise_delay <- FALSE ## should have no effect
    
  } else if (vaccine_dose_scenario == "2_doses_prioritise_delay") {
    
    total_1st_doses_children <- total_doses_children
    total_1st_doses_adults <- ceiling(0.5 * total_doses_adults)
    total_2nd_doses_adults <- total_doses_adults - total_1st_doses_adults
    adults_prioritise_delay <- TRUE
    
  } else if (vaccine_dose_scenario == "2_doses_prioritise_1st") {
    
    total_1st_doses_children <- total_doses_children
    total_1st_doses_adults <- ceiling(0.5 * total_doses_adults)
    total_2nd_doses_adults <- total_doses_adults - total_1st_doses_adults
    adults_prioritise_delay <- FALSE
    
  } else {
    stop("vaccine_dose_scenario not supported")
  }
  
  dose_schedule <- build_dose_schedule(total_1st_doses_adults,
                                       total_2nd_doses_adults,
                                       total_1st_doses_children,
                                       days_between_doses,
                                       doses_per_day_children,
                                       doses_per_day_adults,
                                       adults_prioritise_delay)
  
  daily_doses_children_time <- dose_schedule$daily_doses_children_time
  daily_doses_children_value <- dose_schedule$daily_doses_children_value
  daily_doses_adults_time <- dose_schedule$daily_doses_adults_time
  daily_doses_adults_value <- dose_schedule$daily_doses_adults_value
  
  
  children <- any(daily_doses_children_value) > 0
  adults <- any(daily_doses_adults_value) > 0
  
  ## now we add in the extra padding to allow time to full efficacy 
  if (children) {
    daily_doses_children_time <- c(1,daily_doses_children_time+t_ve) + t1 - 1
    daily_doses_children_value <- cbind(rep(0,idx$dim$vax),
                                        daily_doses_children_value)
  }
  
  if (adults) {
    daily_doses_adults_time <- c(1,daily_doses_adults_time+t_ve) + t1 - 1
    daily_doses_adults_value <- cbind(rep(0,idx$dim$vax),
                                      daily_doses_adults_value)     
  }
  
  if(children){
    int_doses_children <- interpolate_daily_doses(daily_doses_children_time,
                                                  daily_doses_children_value)
    
    if(sum(int_doses_children)!= total_doses_children){
      warning("total doses allocated to children does not equal total doses planned to give to children")
      
    }
  }
  
  if(adults){
    int_doses_adults <- interpolate_daily_doses(daily_doses_adults_time,
                                                daily_doses_adults_value)
    if(!(rowSums(int_doses_adults)[idx$vax$unvaccinated]== total_1st_doses_adults)){
      warning("total 1st doses allocated to adults does not equal total 1st doses planned to give to adults")
    }
    
    if(!(rowSums(int_doses_adults)[idx$vax$one_dose]==total_2nd_doses_adults)){
      warning("total 2nd doses allocated to adults does not equal total 2nd doses planned to give to adults")
    }
    
    if(any(cumsum(int_doses_adults[idx$vax$one_dose,])>cumsum(int_doses_adults[idx$vax$unvaccinated,]))){
      stop("2nd doses being given out without enough 1st doses being given")
    }
    
  }    
  
  
  ## needed for the time_run output
  
  if(children){
    vaccination_campaign_length_children <- ncol(int_doses_children)
  } else{
    vaccination_campaign_length_children <- 0
  }
  
  if(adults){
    vaccination_campaign_length_adults <- ncol(int_doses_adults)
  } else{
    vaccination_campaign_length_adults <- 0
  }  
  
  
  ## prioritisation scenario
  if(tolower(vaccine_used)=="mix"){
    ## children
    prioritisation_strategy_children <- as.matrix(
      readxl::read_excel(
        "vaccine_prioritisation_scenarios.xlsx",
        sheet=prioritisation_children,
        col_names=FALSE
      ),
      nrow=run_pars$n_group,
      byrow=FALSE)
    
    if(prioritisation_adults=="adults_children_first"){
      stop("prioritisation_adults 'adults_children_first' is only valid when vaccine_used='MVA-BN'")
    }
    
    # adults
    prioritisation_strategy_adults <- as.matrix(
      readxl::read_excel(
        "vaccine_prioritisation_scenarios.xlsx",
        sheet=prioritisation_adults,
        col_names=FALSE
      ),
      nrow=run_pars$n_group,
      byrow=FALSE)
    
  }else if (tolower(vaccine_used)=="lc16m8"){
    
    ## children
    prioritisation_strategy_children <- as.matrix(
      readxl::read_excel(
        "vaccine_prioritisation_scenarios.xlsx",
        sheet=prioritisation_children,
        col_names=FALSE
      ),
      nrow=run_pars$n_group,
      byrow=FALSE)
    
    prioritisation_strategy_children<-cbind(rep(1,nrow(prioritisation_strategy_children)))
    #Don't vaccinate SW
    prioritisation_strategy_children[c(17,18)] <- 0
    prioritisation_strategy_adults<-cbind(rep(0, nrow(prioritisation_strategy_children)))
    
  }else if(tolower(vaccine_used)=="mva-bn"){
    
  
    if(prioritisation_adults=="adults_all_equal"){
      # adults
      prioritisation_strategy_adults <- as.matrix(
        readxl::read_excel(
          "vaccine_prioritisation_scenarios.xlsx",
          sheet=prioritisation_adults,
          col_names=FALSE
        ),
        nrow=run_pars$n_group,
        byrow=FALSE)
      # children and adults equal priority
      prioritisation_strategy_adults<-cbind(rep(1,nrow(prioritisation_strategy_adults)))
      
    } else if (prioritisation_adults=="adults_children_first"){
      # adults
      prioritisation_strategy_adults <- as.matrix(
        readxl::read_excel(
          "vaccine_prioritisation_scenarios.xlsx",
          sheet="adults_all_equal",
          col_names=FALSE
        ),
        nrow=run_pars$n_group,
        byrow=FALSE)
      # prioritise children
      prioritisation_strategy_adults<-cbind(c(rep(1,2),rep(0,nrow(prioritisation_strategy_adults)-2)),rep(1,nrow(prioritisation_strategy_adults)))
      
    }else{
      # adults
      prioritisation_strategy_adults <- as.matrix(
        readxl::read_excel(
          "vaccine_prioritisation_scenarios.xlsx",
          sheet=prioritisation_adults,
          col_names=FALSE
        ),
        nrow=run_pars$n_group,
        byrow=FALSE)
      # prioritise adults specified
      prioritisation_strategy_adults<-cbind(prioritisation_strategy_adults[,1],rep(1,nrow(prioritisation_strategy_adults)))
    }
    prioritisation_strategy_children<-cbind(rep(0,nrow(prioritisation_strategy_adults)))
    
  }
  
  
  # hesistancy parameters
  # read in the hesistancy sheet (now a shared resource)
  vax_hes_master <- readxl::read_xlsx("vaccine_prioritisation_scenarios.xlsx",
                                      sheet="hesistancy_params")
  
  vax_hes <- as.numeric(unlist(vax_hes_master[,run_pars$region])) * uptake_realised
  
  
  ## set the targets using vax hesistancy
  ## children
  prioritisation_strategy_children <- prioritisation_strategy_children*vax_hes
  
  ## adults
  prioritisation_strategy_adults <- prioritisation_strategy_adults*vax_hes
  ## logic here to check we are on the same version of adults and children
  ## if there are any non-zero values in any of the rows here which doesn't correspond to is_child then have an issue 
  
  ## Update vaccine efficiency based on vaccine_used
  ve_I <- run_pars$ve_I
  
  ## logic here to check we are on the same version of adults and children
  ## if there are any non-zero values in any of the rows here which doesn't correspond to is_child then have an issue 
  if(tolower(vaccine_used)=="mva-bn"){
    # if using MVA-BN is_child is all FALSE
    is_child <- rep(0,length(mpoxseir:::get_group_bins()$children))
    
  }else if (tolower(vaccine_used)=="lc16m8"){
    # if using Lc16m8 is_child is all TRUE
    is_child <- rep(1,length(mpoxseir:::get_group_bins()$children))
    ve_I[as.logical(is_child),3] <- 0.95
    
  }else if (tolower(vaccine_used)=="mix"){
    # if using Lc16m8 AND MVA-BN  
    is_child <- mpoxseir:::get_group_bins()$children 
    ve_I[as.logical(is_child),3] <- 0.95
    
  }else{
    stop("parameter vaccine used not recognised, please use 'MVA-BN', 'Lc16m8' or 'mix'.")
  }
  
  
  if(sum(is_child)!=sum(ceiling(prioritisation_strategy_children[,ncol(prioritisation_strategy_children)]))){
    warning("check groups that have been assigned as children/adults, there might have been a model update")
  }
  
  if(sum((1-is_child))!=sum(ceiling(prioritisation_strategy_adults[,ncol(prioritisation_strategy_adults)]))){
    warning("check groups that have been assigned as children/adults, there might have been a model update")
  }
  
  # number of prioritisation steps
  N_prioritisation_steps_children <- ncol(prioritisation_strategy_children)
  N_prioritisation_steps_adults <- ncol(prioritisation_strategy_adults)
  
  ## output the new model inputs related to vaccination
  
  record_inputs <- list(
    run_pars = run_pars,
    total_doses_children = total_doses_children,
    total_doses_adults = total_doses_adults,
    doses_per_day_total = doses_per_day_total,
    daily_vax_split_children = daily_vax_split_children,
    vaccine_dose_scenario = vaccine_dose_scenario,
    days_between_doses = days_between_doses,
    prioritisation_children = prioritisation_children,
    prioritisation_adults = prioritisation_adults,
    uptake_realised = uptake_realised,
    t_ve = t_ve,
    time_to_run = max(c(daily_doses_children_time,daily_doses_adults_time))
  )
  
  model_vax_inputs <- list(
    # doses
    daily_doses_children_value = daily_doses_children_value,
    daily_doses_children_time = daily_doses_children_time,
    daily_doses_adults_value = daily_doses_adults_value,
    daily_doses_adults_time = daily_doses_adults_time,
    # prioritisation
    prioritisation_strategy_children = prioritisation_strategy_children,
    prioritisation_strategy_adults = prioritisation_strategy_adults,
    # number of prioritisation steps
    N_prioritisation_steps_children = N_prioritisation_steps_children,
    N_prioritisation_steps_adults = N_prioritisation_steps_adults,
    # vaccine specs
    is_child = is_child,
    ve_I = ve_I
  )
  
  return(list(record_inputs = record_inputs,
              model_vax_inputs = model_vax_inputs))
  
}


build_dose_schedule <- function(total_1st_doses_adults,
                                total_2nd_doses_adults,
                                total_1st_doses_children,
                                days_between_doses,
                                doses_per_day_children,
                                doses_per_day_adults,
                                adults_prioritise_delay) {
  
  remaining_doses_adults <- c(total_1st_doses_adults, total_2nd_doses_adults)
  remaining_doses_children <- total_1st_doses_children
  remaining_doses <- remaining_doses_children + sum(remaining_doses_adults)
  
  t <- 0
  daily_doses_adults <- c()
  daily_doses_children <- c()
  
  adults_eligible_for_2nd_dose <- 0
  
  while (remaining_doses > 0) {
    t <- t + 1
    
    if (t > days_between_doses) {
      ## How many adults newly eligible for 2nd dose?
      adults_eligible_for_2nd_dose <- adults_eligible_for_2nd_dose + 
        daily_doses_adults[1, t - days_between_doses]
    }
    
    ## Allocate excess doses from children to adults
    excess_doses_children <- 
      max(0, doses_per_day_children - remaining_doses_children)
    if (excess_doses_children > 0) {
      doses_per_day_adults_t <- doses_per_day_adults + excess_doses_children
      doses_per_day_children_t <- doses_per_day_children - excess_doses_children
    } else {
      doses_per_day_adults_t <- doses_per_day_adults
      doses_per_day_children_t <- doses_per_day_children
    }
    
    ## Allocate excess doses from adults to children
    excess_doses_adults <- 
      max(0, doses_per_day_adults_t - remaining_doses_adults[1] - 
            min(remaining_doses_adults[2], adults_eligible_for_2nd_dose))
    if (excess_doses_adults > 0) {
      doses_per_day_adults_t <- doses_per_day_adults_t - excess_doses_adults
      doses_per_day_children_t <- doses_per_day_children_t + excess_doses_adults
    } 
    
    ## Allocate doses for children
    doses_children_t <- rep(0, 2)
    doses_children_t[1] <- 
      min(doses_per_day_children_t, remaining_doses_children)
    daily_doses_children <- cbind(daily_doses_children, doses_children_t)
    remaining_doses_children <- remaining_doses_children - sum(doses_children_t)
    
    
    ## Allocate doses for adults
    doses_adults_t <- rep(0, 2) 
    
    if (adults_prioritise_delay) {
      ## allocate 2nd doses
      doses_adults_t[2] <- min(doses_per_day_adults_t,
                               remaining_doses_adults[2],
                               adults_eligible_for_2nd_dose)
      ## allocate 1st doses
      doses_adults_t[1] <- min(doses_per_day_adults_t - doses_adults_t[2],
                               remaining_doses_adults[1])
    } else {
      ## allocate 1st doses
      doses_adults_t[1] <- min(doses_per_day_adults_t, 
                               remaining_doses_adults[1])
      ## allocate 2nd doses
      doses_adults_t[2] <- min(doses_per_day_adults_t - doses_adults_t[1],
                               remaining_doses_adults[2],
                               adults_eligible_for_2nd_dose)
    }
    
    daily_doses_adults <- cbind(daily_doses_adults, doses_adults_t)
    remaining_doses_adults <- remaining_doses_adults - doses_adults_t
    
    ## Remove 2nd dose allocation from 2nd dose eligible
    adults_eligible_for_2nd_dose <- 
      adults_eligible_for_2nd_dose - doses_adults_t[2]
    
    remaining_doses <- remaining_doses_children + sum(remaining_doses_adults)
  }
  
  ## Add on last day with 0s to indicate end of vaccine campaign
  daily_doses_adults <- cbind(daily_doses_adults, rep(0, 2))
  daily_doses_children <- cbind(daily_doses_children, rep(0, 2))
  
  colnames(daily_doses_children) <- NULL
  colnames(daily_doses_adults) <- NULL
  
  
  idx <- mpoxseir::get_compartment_indices()
  
  daily_doses_children_full <- array(0, c(4, ncol(daily_doses_children)))
  daily_doses_children_full[c(idx$vax$unvaccinated, idx$vax$one_dose), ] <-
    daily_doses_children
  daily_doses_adults_full <- array(0, c(4, ncol(daily_doses_adults)))
  daily_doses_adults_full[c(idx$vax$unvaccinated, idx$vax$one_dose), ] <-
    daily_doses_adults
  
  ## Convert to format for interpolate
  daily_doses_adults <- uninterpolate_daily_doses(daily_doses_adults_full)
  daily_doses_children <- uninterpolate_daily_doses(daily_doses_children_full)
  
  list(daily_doses_adults_time = daily_doses_adults$time,
       daily_doses_adults_value = daily_doses_adults$value,
       daily_doses_children_time = daily_doses_children$time,
       daily_doses_children_value = daily_doses_children$value)
}
