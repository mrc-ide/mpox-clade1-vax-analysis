save_index <- function(vax = FALSE, rt = FALSE) {
  
  ## baseline index saving
  idx <- c("cases_inc",
           "deaths_inc",
           "cases_inc_00_04",
           "cases_inc_05_14",
           "cases_inc_15_plus", 
           "cases_cumulative_00_04",
           "cases_cumulative_05_14",
           "cases_cumulative_15_plus", 
           "cases_inc_PBS",
           "cases_inc_CSW",
           "cases_inc_ASW",
           "cases_inc_SW",
           "cases_inc_HCW",
           "deaths_inc_00_04",
           "deaths_inc_05_14",
           "deaths_inc_15_plus", 
           "deaths_cumulative_00_04",
           "deaths_cumulative_05_14",
           "deaths_cumulative_15_plus", 
           "deaths_inc_PBS",
           "deaths_inc_SW",
           "deaths_inc_CSW",
           "deaths_inc_ASW",
           "deaths_inc_HCW",
           "observed_cases_inc",
           "observed_cases_inc_00_04",
           "observed_cases_inc_05_14",
           "observed_cases_inc_15_plus",
           "observed_cases_inc_PBS",
           "observed_cases_inc_CSW",
           "observed_cases_inc_ASW",
           "observed_cases_inc_SW",
           "observed_cases_inc_HCW",
           "S_tot",
           "E_tot",
           "I_tot",
           "R_tot",
           "D_tot",
           "N_tot"
  )
  
  names <- idx
  
  if (vax) {
    idx <- c(idx,
             "dose1_inc",
             "dose1_inc_00_04",
             "dose1_inc_05_14",
             "dose1_inc_15_plus", 
             "dose1_inc_PBS",
             "dose1_inc_CSW",
             "dose1_inc_ASW",
             "dose1_inc_SW",
             "dose1_inc_HCW",
             "dose2_inc",
             "dose2_inc_00_04",
             "dose2_inc_05_14",
             "dose2_inc_15_plus", 
             "dose2_inc_PBS",
             "dose2_inc_CSW",
             "dose2_inc_ASW",
             "dose2_inc_SW",
             "dose2_inc_HCW",
             "dose1_cumulative",
             "dose1_cumulative_00_04",
             "dose1_cumulative_05_14",
             "dose1_cumulative_15_plus", 
             "dose1_cumulative_PBS",
             "dose1_cumulative_CSW",
             "dose1_cumulative_ASW",
             "dose1_cumulative_SW",
             "dose1_cumulative_HCW",
             "dose2_cumulative",
             "dose2_cumulative_00_04",
             "dose2_cumulative_05_14",
             "dose2_cumulative_15_plus", 
             "dose2_cumulative_PBS",
             "dose2_cumulative_CSW",
             "dose2_cumulative_ASW",
             "dose2_cumulative_SW",
             "dose2_cumulative_HCW")
  }
  
  names <- idx
  
  if (rt) {
    S_names <- get_S_names()
    
    idx <- c(idx, "S")
    names <- c(names, S_names)
  }
  
  
  list(idx = idx,
       names = names)
}
