create_packer <- function(region, assumptions, mixing_matrix,
                          use_ve_D = FALSE, overrides = list()) {

  fitted_pars <- c("alpha_cases",
                   "alpha_deaths",
                   "beta_z_max",
                   "cfr_00_04",
                   "lambda",
                   "phi_05_14",
                   "phi_15_plus",
                   "prop_SW",
                   "R0_hh",
                   "R0_sw_st")
  
  if (region == "sudkivu") {
    fitted_pars <- fitted_pars[fitted_pars != "beta_z_max"]
  }
  
  if (region == "both") {
    groups <- c("equateur", "sudkivu")
    
    fixed <- lapply(groups, 
                    function(r) get_fixed(r, assumptions, mixing_matrix,
                                          use_ve_D, overrides))
    names(fixed) <- groups
    
    ## beta_z_max is "shared" but only applies to equateur
    shared <- c("beta_z_max", "R0_hh")
    
    packer <- monty::monty_packer_grouped(
      groups, scalar = fitted_pars, fixed = fixed, process = process,
      shared = shared)
  } else {
    fixed <- get_fixed(region, assumptions, mixing_matrix, use_ve_D, overrides)
    
    packer <- monty::monty_packer(
      scalar = fitted_pars, fixed = fixed, process = process)
  }
  
  packer
}

get_fixed <- function(region, assumptions, mixing_matrix,
                      use_ve_D = FALSE, overrides = list()) {
  
  fixed <- 
    mpoxseir::parameters_fixed(region = region,
                               mixing_matrix = mixing_matrix,
                               p_HCW = 0,
                               initial_infections = 1,
                               use_ve_D = use_ve_D,
                               overrides = overrides)
  fixed$R0_hh <- NULL
  fixed$R0_sw_st <- NULL
  fixed$beta_z_max <- NULL
  fixed$beta_hcw <- 0
  
  fixed$alpha_cases <- NULL
  fixed$alpha_cases_00_04 <- NULL
  fixed$alpha_cases_05_14 <- NULL
  fixed$alpha_cases_15_plus <- NULL

  fixed$alpha_deaths <- NULL
  fixed$alpha_deaths_00_04 <- NULL
  fixed$alpha_deaths_05_14 <- NULL
  fixed$alpha_deaths_15_plus <- NULL
  
  ## phi_00_04 is keep as default value of 1
  ## we will fit the rest as phi_05_14 or phi_15_plus
  fixed$phi_05_14 <- NULL
  fixed$phi_15_plus <- NULL
  fixed$phi_CSW_12_14 <- NULL
  fixed$phi_CSW_15_17 <- NULL
  fixed$phi_ASW <- NULL
  fixed$phi_PBS <- NULL
  fixed$phi_HCW <- NULL
  
  fixed$mixing_matrix <- mixing_matrix
  
  fixed$N0 <- NULL
  fixed$S0 <- NULL
  fixed$m_sex <- NULL
  fixed$m_gen_pop <- NULL
  
  names(fixed)[names(fixed) == "CFR"] <- "CFR_0"
  
  names(fixed)[names(fixed) == "seed_rate"] <- "seed_profile"
  
  fixed
}

## fitted params are:
## beta_z_max = zoonotic beta for the age-group with the highest risk
##     (used to calculate RR_z)
## R0_hh = R0 for household transmission (used to calculate beta_h)
## R0_sw_st = R0 for sex workers to people who buy sex)

process <- function(x) {
  
  # Converting R0_hh to the beta_hh parameter given mixing matrix
  ## (excluding SW & PBS)
  
  ## Calc. duration of infectiousness by age, weighted by disease severity
  ## We assume here that the R0 calculation that duration_infectious_by_age
  ## is used in is based on duration_infectious_by_age in unvaccinated
  ## individuals (i.e. with the unvaccinated CFR)
  k <- 1 # number of compartments per infectious disease state
  dt <- 1 # timestep
  idx <- mpoxseir::get_compartment_indices()
  
  ## Have to calculate some parameters that are affected by fitting prop_SW
  demographic_params <-
    mpoxseir::parameters_demographic(region = x$region,
                                     mixing_matrix = x$mixing_matrix,
                                     p_SW = x$prop_SW * 0.5,
                                     p_HCW = 0)
  
  N <- demographic_params$province_pop[[x$region]]
  N0 <- round(N * demographic_params$N0 / sum(demographic_params$N0))
  p_unvaccinated <- demographic_params$p_unvaccinated
  S0 <- matrix(0, nrow = x$n_group, ncol = x$n_vax)
  S0[, idx$vax$unvaccinated] <- round(N0 * demographic_params$p_unvaccinated)
  S0[, idx$vax$historic] <- N0 - S0[, idx$vax$unvaccinated]
  
  ## Calculate CFR from fitted cfr_00_04 and relative risk
  RR <- x$CFR_0 / x$CFR_0[idx$group$`0-4`, idx$vax$unvaccinated]
  CFR <- x$cfr_00_04 * RR
  
  duration_infectious_by_age <-
    CFR[, idx$vax$unvaccinated] * ((k * dt) / (1 - exp(-x$gamma_Id * dt))) +
    (1 - CFR[, idx$vax$unvaccinated]) * ((k * dt) / (1 - exp(-x$gamma_Ir * dt)))
  
  
  ## Calculate beta_household given R0, mixing matrix and duration of infectiousness
  m_gen_pop <- demographic_params$m_gen_pop
  beta_h <- x$R0_hh / 
    Re(eigen(m_gen_pop * duration_infectious_by_age)$values[1])
  
  # Converting the relative risk age-spline to the age-specific beta_z
  if (x$region == "equateur") {
    beta_z <- x$RR_z * x$beta_z_max / (sum(x$RR_z * N0) / (sum(N0) / 1e6))
  } else {
    beta_z <- x$RR_z * 0
  }
  
  # Converting the inputted R0_sw_st into the rates required for the mixing matrix
  # we never resolve the product m[i, j] = beta_s * c[i, j] here, so beta_s is
  # set to 1. I.e. we are working with transmission rates, rather than transmission
  # probabilities + contact rates.
  
  N_CSW <- N0["CSW"]
  N_ASW <- N0["ASW"]
  N_PBS <- N0["PBS"]
  
  Delta_PBS <- duration_infectious_by_age[idx$group$PBS]
  
  # transmission rate from CSW / ASW to PBS is the same
  if (N_CSW + N_ASW == 0) {
    m_csw_pbs <- 0
  } else {
    m_csw_pbs <- x$R0_sw_st / Delta_PBS * N_PBS / (N_CSW + N_ASW)
  }
  m_asw_pbs <- m_csw_pbs
  # transmission from PBS -> CSW / ASW is proportional to their group sizes
  m_pbs_csw <- m_csw_pbs * N_CSW / N_PBS
  m_pbs_asw <- m_csw_pbs * N_ASW / N_PBS
  
  m_sex <- demographic_params$m_sex
  
  m_sex[idx$group$CSW, idx$group$PBS] <- m_csw_pbs
  m_sex[idx$group$ASW, idx$group$PBS] <- m_asw_pbs
  m_sex[idx$group$PBS, idx$group$CSW] <- m_pbs_csw
  m_sex[idx$group$PBS, idx$group$ASW] <- m_pbs_asw
  
  beta_s <- 1
  
  seed_rate <- x$lambda * x$seed_profile
  
  list(beta_h = beta_h,
       beta_z = beta_z,
       N0 = N0,
       S0 = S0,
       m_sex = m_sex,
       m_gen_pop = m_gen_pop,
       beta_s = beta_s,
       CFR = CFR,
       alpha_cases_00_04 = x$alpha_cases,
       alpha_cases_05_14 = x$alpha_cases,
       alpha_cases_15_plus = x$alpha_cases,
       alpha_deaths_00_04 = x$alpha_deaths,
       alpha_deaths_05_14 = x$alpha_deaths,
       alpha_deaths_15_plus = x$alpha_deaths,
       seed_rate = seed_rate,
       phi_CSW_12_14 = x$phi_05_14,
       phi_CSW_15_17 = x$phi_15_plus,
       phi_ASW = x$phi_15_plus,
       phi_PBS = x$phi_15_plus,
       phi_HCW = x$phi_15_plus)
  
}
