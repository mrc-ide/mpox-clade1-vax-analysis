
orderly2::orderly_parameters(short_run = FALSE,
                             fit_by_age = TRUE,
                             deterministic = FALSE,
                             mixing_matrix = "Zimbabwe",
                             fit_KPs = TRUE, 
                             assumptions = "fix_prop_SW",
                             vaccines_onset = "start")

# Load in all dependencies 
orderly2::orderly_dependency(
  name = "pmcmc",
  query = "latest(parameter:region=='sudkivu' && parameter:assumptions == 'fix_prop_SW')",
  files = c("inputs/samples_sudkivu.rds" = "outputs/samples.rds",
            "inputs/fitting_data_sudkivu.rds" = "outputs/fitting_data.rds"))
orderly2::orderly_dependency(
  name = "pmcmc",
  query = "latest(parameter:region=='equateur'  &&  parameter:assumptions == 'fix_prop_SW')",
  files = c("inputs/samples_equateur.rds" = "outputs/samples.rds",
            "inputs/fitting_data_equateur.rds" = "outputs/fitting_data.rds"))
orderly2::orderly_dependency(name = "run_postprocess",
                             "latest(parameter:region =='sudkivu' && 
                                     parameter:deterministic == FALSE && 
                                     parameter:use_both_fit == FALSE && 
                                     parameter:fit_by_age == TRUE && 
                                     parameter:fit_KPs == TRUE &&
                                     parameter:short_run == FALSE && 
                                     parameter:vaccines_onset == 'start' &&
                             parameter:assumptions == 'fix_prop_SW')",
                             files=c("inputs/runs_postprocess_sudkivu.RDS" = "outputs/runs_postprocess.RDS"))
orderly2::orderly_dependency(name = "run_postprocess",
                             "latest(parameter:region =='equateur' && 
                                     parameter:deterministic == FALSE && 
                                     parameter:use_both_fit == FALSE && 
                                     parameter:fit_by_age == TRUE && 
                                     parameter:fit_KPs == TRUE &&
                                     parameter:short_run == FALSE && 
                                     parameter:vaccines_onset == 'start' &&
                             parameter:assumptions == 'fix_prop_SW')",
                             files=c("inputs/runs_postprocess_equateur.RDS" = "outputs/runs_postprocess.RDS"))
orderly2::orderly_dependency(name = "run_postprocess_burundi",
                             "latest(parameter:deterministic == FALSE && 
                             parameter:short_run == FALSE && 
                             parameter:assumptions == 'fix_prop_SW')",
                             files=c("inputs/runs_postprocess_burundi.RDS" = "outputs/runs_postprocess.RDS"))

orderly2::orderly_dependency(name = "pmcmc_burundi",
                             "latest(parameter:deterministic == FALSE && parameter:short_run == FALSE &&
                             parameter:assumptions == 'fix_prop_SW')",
                             files =  c("inputs/samples.rds" = "outputs/samples.rds",
                                        "inputs/fitting_data.rds" = "outputs/fitting_data.rds"))
# Populations of provinces
province_pop<-mpoxseir::parameters_demographic(region="sudkivu")$province_pop

# Create new folder for new plots
dir.create("outputs/rerun")

cols_region <- c(MetBrewer::met.brewer("Cassatt1", 7)[c(5,1)],MetBrewer::met.brewer("Demuth", 3)[3], MetBrewer::met.brewer("Cassatt1", 7)[3]) 
names(cols_region) <- c("equateur", "sudkivu", "both", "burundi")
nms_region <- c(equateur = "Equateur", sudkivu = "Sud Kivu", burundi = "Burundi")

vaccines_onset <- "start"

# Dates 
start_date<-if(vaccines_onset=="end"){
  # as.Date("2024-10-06")
  mpoxseir::mpoxseir_date_as_date(770) # needs updating if fits changed
}else{
  # (as.Date("2024-10-06") - 644 + 371)
  mpoxseir::mpoxseir_date_as_date(371) #+ 14
} 

library(tidyverse)

######## Colour Palettes ########

scenarios_palette_vaccine <- c("No vaccination" = "black",
                               "LC16m8" = "#75884BFF",
                               "LC16m8 & MVA-BN" = "#5B859EFF",
                               "One dose MVA-BN"= "#800000",
                               "Two doses MVA-BN"= "#B38711FF",
                               "One dose fractional MVA-BN" = "#D48F90FF")
#  #AB84A5FF"

# all scenario names A-F
scenario_names_plots <- c("A: LC16m8", 
                          "B: LC16m8 + kids first",
                          "C: LC16m8 & MVA-BN",
                          "D: LC16m8 & MVA-BN + kids first", 
                          "E: LC16m8 & MVA-BN + SWs first",
                          "F: MVA-BN", 
                          "G: MVA-BN + kids first", 
                          "H: MVA-BN + SWs first",
                          "I: MVA-BN + fractional",
                          "J: MVA-BN + fractional + kids first",
                          "K: MVA-BN + fractional + SWs first",
                          "L: Two dose MVA-BN"
                          )

# line type
linetype_priority <- c("All equal"="solid",
                       "<12 years first"="dashed",
                       "SWs first"="dotted")
# shape
shape_priority <- c("All equal"=4,
                    "<12 years first"=16,
                    "SWs first"=17)


########################### EQUATEUR PLOTS #####################################
region <-"equateur"

# load in post process
big_df_equateur <-  readRDS("inputs/runs_postprocess_equateur.RDS") |> 
  mutate(Date = start_date + lubridate::days((TimeStep-1) * 7)) |> 
  mutate(scenario_combo=scenario_name, scenario_name=sub(" \\(.*", "", scenario_new), 
         scenario_dose=total_doses_children+total_doses_adults,
         scenario_labels= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "A: LC16m8",
           scenario_new=="LC16m8 + kids first" ~ "B: LC16m8 + kids first",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "C: LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "D: LC16m8 & MVA-BN + kids first",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "E: LC16m8 & MVA-BN + SWs first",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "F: MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "G: MVA-BN + kids first",
           scenario_new=="MVA-BN  + 1st dose only + SWs first" ~ "H: MVA-BN + SWs first",
           scenario_new=="MVA-BN + fractional" ~ "I: MVA-BN + fractional",
           scenario_new=="MVA-BN + fractional + kids first" ~ "J: MVA-BN + fractional + kids first",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "K: MVA-BN + fractional + SWs first",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "L: Two dose MVA-BN",
           scenario_new=="MVA-BN + kids first" ~ "M: Two dose MVA-BN + kids first",
           scenario_new=="MVA-BN + SWs first" ~ "N: Two dose MVA-BN + SWs first"
         ),
         scenario_vaccine= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "LC16m8",
           scenario_new=="LC16m8 + kids first" ~ "LC16m8",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "LC16m8 & MVA-BN",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN  + 1st dose only + SWs first" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + fractional" ~ "One dose fractional MVA-BN",
           scenario_new=="MVA-BN + fractional + kids first" ~ "One dose fractional MVA-BN",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "One dose fractional MVA-BN",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "Two doses MVA-BN",
           scenario_new=="MVA-BN + kids first" ~ "Two doses MVA-BN",
           scenario_new=="MVA-BN + SWs first" ~ "Two doses MVA-BN"
         ),
         scenario_priority= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "All equal",
           scenario_new=="LC16m8 + kids first" ~ "<12 years first",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "All equal",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "<12 years first",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "SWs first",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "All equal",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN  + 1st dose only + SWs first" ~ "SWs first",
           scenario_new=="MVA-BN + fractional" ~ "All equal",
           scenario_new=="MVA-BN + fractional + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "SWs first",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "All equal",
           scenario_new=="MVA-BN + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + SWs first" ~ "SWs first"
         ),
        
         ) %>%
  mutate(
    scenario_dose_all = scenario_dose,
    scenario_dose = case_when(
    scenario_new=="MVA-BN + fractional" ~ scenario_dose/5,
    scenario_new=="MVA-BN + fractional + kids first" ~ scenario_dose/5,
    scenario_new=="MVA-BN + fractional + SWs first" ~ scenario_dose/5,
    .default=scenario_dose
  )) |>
  mutate(scenario_priority = factor(scenario_priority,
                                    levels=c("All equal","<12 years first",
                                             "SWs first")))|>
  mutate(scenario_doses_per_day_character=factor(format(scenario_doses_per_day_total, big.mark = ",", scientific = FALSE),
                                                 levels=c("   250","   500"," 1,000",
                                                          " 2,500"," 5,000","10,000", "20,000")))


# Plot max dose & middle daily dose total scenario
big_df_plots_eq <- filter(big_df_equateur,
                          scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)], 
                          scenario_doses_per_day_character==" 5,000")


big_df_cases_eq <- big_df_plots_eq |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))


# percentage of Available full doses administered for each scenario
dose_percent_equateur<-big_df_cases_eq|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

sensitivity_data_eq <- big_df_plots_eq |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_eq$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

best_case_scenario_eq <- sensitivity_data_eq|>group_by(scenario_dose)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))

sensitivity_data_eq_plot <- sensitivity_data_eq %>%
  mutate(lower_averted_95_per_dose = case_when(
    lower_averted_95_per_dose < 0 ~ 0,
    .default = lower_averted_95_per_dose
  ))


## Equateur supplementary figures

## line plot
big_df_eq <- big_df_plots_eq |> filter(scenario_dose == 520000, ignore.case = TRUE)|>
  mutate(scenario_labels = factor(scenario_labels,
                                  levels = scenario_names_plots[1:length(scenario_names_plots)]),
         scenario_vaccine = factor(scenario_vaccine,
                                   levels = names(scenarios_palette_vaccine)))

line_plot_equateur <- big_df_eq |> 
  filter(Category =="cases_inc") |>
  ggplot() +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1.3,alpha=0.6) +
  geom_line(aes(x = date, y = mean_value, color = scenario_vaccine, linetype = scenario_priority),linewidth=1.3,alpha=0.6) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Equateur (520,000 doses)",
       subtitle = "Weekly infections",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = "Vaccine", 
                     values=c( scenarios_palette_vaccine,
                               "No vaccination" ="black"))+
  scale_linetype_manual(name="Priority", values=linetype_priority)+   
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=14.5),
        legend.direction = "vertical", 
        legend.box = "vertical")+
  guides(linetype = guide_legend(ncol = 3)) + 
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))

line_plot_equateur_deaths <- big_df_eq |> 
  filter(Category =="deaths_inc") |>
  ggplot() +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1.3,alpha=0.6) +
  geom_line(aes(x = date, y = mean_value, color = scenario_vaccine, linetype = scenario_priority),linewidth=1.3,alpha=0.6) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(subtitle = "Weekly deaths",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = "Vaccine", 
                     values=c("No vaccination" ="black",
                              scenarios_palette_vaccine))+
  scale_linetype_manual(name="Priority", values=linetype_priority)+   
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=14.5),
        legend.direction = "vertical", 
        legend.box = "vertical")+
  guides(linetype = guide_legend(ncol = 3,order = 1))+
  guides(color = guide_legend(ncol = 3,order=2))

eq_lines <- cowplot::plot_grid(
  line_plot_equateur + theme(legend.position = "none"),
  line_plot_equateur_deaths,
  nrow=2,
  align="v",
  rel_heights = c(1,1.3))
ggsave(paste0("outputs/rerun/line_plot_equateur.png"), eq_lines,
       width = 10, height = 10)

sensitivity_data_eq %>% filter(scenario_vaccine=="LC16m8", scenario_priority %in% c("All equal","<12 years first"),scenario_dose==260000) 

sensitivity_data_eq %>% filter(scenario_vaccine %in% c("Two doses MVA-BN","One dose MVA-BN"), scenario_priority=="All equal",scenario_dose==520000) 



## equivalent to figure 3 but for deaths 

big_df_deaths_eq <- big_df_plots_eq |>  filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))


# percentage of Available full doses administered for each scenario
dose_percent_equateur_deaths <-big_df_deaths_eq|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

sensitivity_data_eq_deaths <- big_df_plots_eq |> filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_plots_eq$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

best_case_scenario_eq_deaths <- sensitivity_data_eq_deaths|>group_by(scenario_dose)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))

sensitivity_data_eq_plot_deaths <- sensitivity_data_eq_deaths %>%
  mutate(lower_averted_95_per_dose = case_when(
    lower_averted_95_per_dose < 0 ~ 0,
    .default = lower_averted_95_per_dose
  ))





########################### SUD KIVU PLOTS #####################################
  
# load in post process
big_df_sudkivu <-  readRDS("inputs/runs_postprocess_sudkivu.RDS") |> 
  mutate(Date = start_date + lubridate::days((TimeStep-1) * 7)) |> 
  mutate(scenario_combo=scenario_name, scenario_name=sub(" \\(.*", "", scenario_new), 
         scenario_dose=total_doses_children+total_doses_adults,
         scenario_labels= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "A: LC16m8",
           scenario_new=="LC16m8 + kids first" ~ "B: LC16m8 + kids first",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "C: LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "D: LC16m8 & MVA-BN + kids first",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "E: LC16m8 & MVA-BN + SWs first",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "F: MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "G: MVA-BN + kids first",
           scenario_new=="MVA-BN + 1st dose only + SWs first" ~ "H: MVA-BN + SWs first",
           scenario_new=="MVA-BN + fractional" ~ "I: MVA-BN + fractional",
           scenario_new=="MVA-BN + fractional + kids first" ~ "J: MVA-BN + fractional + kids first",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "K: MVA-BN + fractional + SWs first",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "L: Two dose MVA-BN",
           scenario_new=="MVA-BN + kids first" ~ "M: Two dose MVA-BN + kids first",
           scenario_new=="MVA-BN + SWs first" ~ "N: Two dose MVA-BN + SWs first"
         ),
         scenario_vaccine= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "LC16m8",
           scenario_new=="LC16m8 + kids first" ~ "LC16m8",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "LC16m8 & MVA-BN",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + SWs first" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + fractional" ~ "One dose fractional MVA-BN",
           scenario_new=="MVA-BN + fractional + kids first" ~ "One dose fractional MVA-BN",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "One dose fractional MVA-BN",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "Two doses MVA-BN",
           scenario_new=="MVA-BN + kids first" ~ "Two doses MVA-BN",
           scenario_new=="MVA-BN + SWs first" ~ "Two doses MVA-BN"
         ),
         scenario_priority= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "All equal",
           scenario_new=="LC16m8 + kids first" ~ "<12 years first",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "All equal",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "<12 years first",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "SWs first",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "All equal",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + 1st dose only + SWs first" ~ "SWs first",
           scenario_new=="MVA-BN + fractional" ~ "All equal",
           scenario_new=="MVA-BN + fractional + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "SWs first",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "All equal",
           scenario_new=="MVA-BN + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + SWs first" ~ "SWs first"
         ),
         
  ) |>
  mutate(
    scenario_dose_all = scenario_dose,
    scenario_dose = case_when(
      scenario_new=="MVA-BN + fractional" ~ scenario_dose/5,
      scenario_new=="MVA-BN + fractional + kids first" ~ scenario_dose/5,
      scenario_new=="MVA-BN + fractional + SWs first" ~ scenario_dose/5,
      .default=scenario_dose
    )) |>
  mutate(scenario_priority = factor(scenario_priority,
                                    levels=c("All equal","<12 years first",
                                             "SWs first")))|>
  mutate(scenario_doses_per_day_character=factor(format(scenario_doses_per_day_total, big.mark = ",", scientific = FALSE),
                                                 levels=c("   250","   500"," 1,000",
                                                          " 2,500"," 5,000","10,000", "20,000")))

# Plot max dose scenario and middle daily doses 
big_df_plots_sk <- filter(big_df_sudkivu,
                          scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)])|>
  filter(scenario_dose < 3000000, scenario_doses_per_day_character==" 5,000")


big_df_cases_sk <- big_df_plots_sk |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))

dose_percent_sudkivu<-big_df_cases_sk|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

## SUD KIVU LINE PLOT 
big_df_sk <- big_df_plots_sk |> filter(scenario_dose == 2000000, ignore.case = TRUE)|>
  mutate(scenario_labels = factor(scenario_labels,
                                  levels = scenario_names_plots[1:length(scenario_names_plots)]),
         scenario_vaccine = factor(scenario_vaccine,
                                   levels = names(scenarios_palette_vaccine)))

sensitivity_data_sk <- big_df_plots_sk |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_sk$TimeStep, na.rm = TRUE)) |>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

# sensitivity_data_sk <- big_df_plots_sk |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_sk$TimeStep, na.rm = TRUE))|>
#   mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
#                                levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

# Bar plot of best scenario at each dose 
best_case_scenario_sk <- sensitivity_data_sk|>group_by(scenario_dose)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE), levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))

sensitivity_data_sk_plot <- sensitivity_data_sk %>%
  mutate(lower_averted_95_per_dose = case_when(
    lower_averted_95_per_dose < 0 ~ 0,
    .default = lower_averted_95_per_dose
  ))



## Sud Kivu supplementary plots

## line plot
big_df_sk <- big_df_plots_sk |> filter(scenario_dose == 2000000, ignore.case = TRUE)|>
  mutate(scenario_labels = factor(scenario_labels,
                                  levels = scenario_names_plots[1:length(scenario_names_plots)]),
         scenario_vaccine = factor(scenario_vaccine,
                                   levels = names(scenarios_palette_vaccine)))

line_plot_sk <- big_df_sk |> 
  filter(Category =="cases_inc") |>
  ggplot() +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1.3,alpha=0.6) +
  geom_line(aes(x = date, y = mean_value, color = scenario_vaccine, linetype = scenario_priority),linewidth=1.3,alpha=0.6) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Sud Kivu (2,000,000 doses)",
       subtitle = "Weekly infections",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = "Vaccine", 
                     values=c( scenarios_palette_vaccine,
                               "No vaccination" ="black"))+
  scale_linetype_manual(name="Priority", values=linetype_priority)+   
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=14.5),
        legend.direction = "vertical", 
        legend.box = "vertical")+
  guides(linetype = guide_legend(ncol = 3)) + 
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))

line_plot_sk_deaths <- big_df_sk |> 
  filter(Category =="deaths_inc") |>
  ggplot() +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1.3,alpha=0.6) +
  geom_line(aes(x = date, y = mean_value, color = scenario_vaccine, linetype = scenario_priority),linewidth=1.3,alpha=0.6) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(subtitle = "Weekly deaths",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = "Vaccine", 
                     values=c("No vaccination" ="black",
                              scenarios_palette_vaccine))+
  scale_linetype_manual(name="Priority", values=linetype_priority)+   
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=14.5),
        legend.direction = "vertical", 
        legend.box = "vertical")+
  guides(linetype = guide_legend(ncol = 3,order = 1))+
  guides(color = guide_legend(ncol = 3,order=2))

sk_lines <- cowplot::plot_grid(
  line_plot_sk + theme(legend.position = "none"),
  line_plot_sk_deaths,
  nrow=2,
  align="v",
  rel_heights = c(1,1.3))
ggsave(paste0("outputs/rerun/line_plot_sudkivu.png"), sk_lines,
       width = 10, height = 10)



## equivalent to Fig 4 but for deaths

big_df_deaths_sk <- big_df_plots_sk |>  filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))

# dose_percent_sudkivu<-big_df_cases_sk|>
#   group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
#   summarise(dose_percent=mean_cumulative_doses/scenario_dose)

# ## SUD KIVU LINE PLOT 
# big_df_sk <- big_df_plots_sk |> filter(scenario_dose == 2000000, ignore.case = TRUE)|>
#   mutate(scenario_labels = factor(scenario_labels,
#                                   levels = scenario_names_plots[1:9]),
#          scenario_vaccine = factor(scenario_vaccine,
#                                    levels = names(scenarios_palette_vaccine)))

sensitivity_data_sk_deaths <- big_df_plots_sk |> filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_plots_sk$TimeStep, na.rm = TRUE)) |>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

# sensitivity_data_sk <- big_df_plots_sk |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_sk$TimeStep, na.rm = TRUE))|>
#   mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
#                                levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

# Bar plot of best scenario at each dose 
best_case_scenario_sk_deaths <- sensitivity_data_sk_deaths|>group_by(scenario_dose)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE), levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))

sensitivity_data_sk_plot_deaths <- sensitivity_data_sk_deaths %>%
  mutate(lower_averted_95_per_dose = case_when(
    lower_averted_95_per_dose < 0 ~ 0,
    .default = lower_averted_95_per_dose
  ),
  mean_averted_per_dose = case_when(
    mean_averted_per_dose < 0 ~ 0.000001,
    .default = mean_averted_per_dose
  )
  )





########################### BURUNDI: BUJUMBURA PLOTS #####################################
region <-"bujumbura"

# load in post process
big_df_bujumbura <-  readRDS("inputs/runs_postprocess_burundi.RDS") |> 
  mutate(Date = start_date + lubridate::days((TimeStep-1) * 7)) |> 
  mutate(scenario_combo=scenario_name, scenario_name=sub(" \\(.*", "", scenario_new), 
         scenario_dose=total_doses_children+total_doses_adults,
         scenario_labels= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "A: LC16m8",
           scenario_new=="LC16m8 + kids first" ~ "B: LC16m8 + kids first",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "C: LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "D: LC16m8 & MVA-BN + kids first",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "E: LC16m8 & MVA-BN + SWs first",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "F: MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "G: MVA-BN + kids first",
           scenario_new=="MVA-BN + 1st dose only + SWs first" ~ "H: MVA-BN + SWs first",
           scenario_new=="MVA-BN + fractional" ~ "I: MVA-BN + fractional",
           scenario_new=="MVA-BN + fractional + kids first" ~ "J: MVA-BN + fractional + kids first",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "K: MVA-BN + fractional + SWs first",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "L: Two dose MVA-BN",
           scenario_new=="MVA-BN + kids first" ~ "M: Two dose MVA-BN + kids first",
           scenario_new=="MVA-BN + SWs first" ~ "N: Two dose MVA-BN + SWs first"
         ),
         scenario_vaccine= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "LC16m8",
           scenario_new=="LC16m8 + kids first" ~ "LC16m8",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "LC16m8 & MVA-BN",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "LC16m8 & MVA-BN",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + 1st dose only + SWs first" ~ "One dose MVA-BN",
           scenario_new=="MVA-BN + fractional" ~ "One dose fractional MVA-BN",
           scenario_new=="MVA-BN + fractional + kids first" ~ "One dose fractional MVA-BN",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "One dose fractional MVA-BN",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "Two doses MVA-BN",
           scenario_new=="MVA-BN + kids first" ~ "Two doses MVA-BN",
           scenario_new=="MVA-BN + SWs first" ~ "Two doses MVA-BN"
         ),
         scenario_priority= case_when(
           # LC16m8 to everyone
           scenario_new=="LC16m8" ~ "All equal",
           scenario_new=="LC16m8 + kids first" ~ "<12 years first",
           # LC16m8 to uner 12s & one dose MVA-BN to over 12s
           scenario_new=="LC16m8 + MVA-BN >12y" ~ "All equal",
           scenario_new=="LC16m8 + MVA-BN >12y + kids first" ~ "<12 years first",
           scenario_new=="LC16m8 + MVA-BN >12y + SWs first" ~ "SWs first",
           # One dose MVA-BN to everyone
           scenario_new=="MVA-BN + 1st dose only" ~ "All equal",
           scenario_new=="MVA-BN + 1st dose only + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + 1st dose only + SWs first" ~ "SWs first",
           scenario_new=="MVA-BN + fractional" ~ "All equal",
           scenario_new=="MVA-BN + fractional + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + fractional + SWs first" ~ "SWs first",
           # Two doses of MVA-BN to everyone
           scenario_new=="MVA-BN" ~ "All equal",
           scenario_new=="MVA-BN + kids first" ~ "<12 years first",
           scenario_new=="MVA-BN + SWs first" ~ "SWs first"
         ),
         
  ) |>
  mutate(
    scenario_dose_all = scenario_dose,
    scenario_dose = case_when(
      scenario_new=="MVA-BN + fractional" ~ scenario_dose/5,
      scenario_new=="MVA-BN + fractional + kids first" ~ scenario_dose/5,
      scenario_new=="MVA-BN + fractional + SWs first" ~ scenario_dose/5,
      .default=scenario_dose
    )) |>
  mutate(scenario_priority = factor(scenario_priority,
                                    levels=c("All equal","<12 years first",
                                             "SWs first")))%>%
  mutate(scenario_doses_per_day_character=factor(format(scenario_doses_per_day_total, big.mark = ",", scientific = FALSE),
                                                 levels=c("   25","   50","  100","  250","  500","1,000", "2,000")))


# Plot max dose scenario
big_df_plots_bu <- filter(big_df_bujumbura,
                          scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)])|>
  filter(scenario_dose < (50000 + 1))


big_df_cases_bu <- big_df_plots_bu |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_bujumbura$TimeStep, na.rm = TRUE))


# percentage of Available full doses administered for each scenario
dose_percent_bujumbura<-big_df_cases_bu|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

## Bujumbura LINE PLOT 
big_df_bu <- big_df_plots_bu |> 
  filter(scenario_dose == 50000, ignore.case = TRUE)|>
  mutate(scenario_labels = factor(scenario_labels,
                                  levels = scenario_names_plots[1:length(scenario_names_plots)]),
         scenario_vaccine = factor(scenario_vaccine,
                                   levels = names(scenarios_palette_vaccine)))

sensitivity_data_bu <- big_df_plots_bu |>
  filter(Category %in% c("cases_inc_cum"),
         TimeStep==max(big_df_plots_bu$TimeStep, na.rm = TRUE),
         scenario_doses_per_day_character=="  500")|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))


# Bar plot of best scenario at each dose 
best_case_scenario_bu <- sensitivity_data_bu|>
  group_by(scenario_dose)|>
  top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE), levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))

sensitivity_data_bu_plot <- sensitivity_data_bu %>%
  mutate(lower_averted_95_per_dose = case_when(
    lower_averted_95_per_dose < 0 ~ 0.000001,
    .default = lower_averted_95_per_dose
  ),
  mean_averted_per_dose = case_when(
    mean_averted_per_dose < 0 ~ 0.000001,
    .default = mean_averted_per_dose
  ))



### Bujumbura supplementary figures

line_plot_bujumbura <- big_df_bu |> 
  filter(Category =="cases_inc", 
         scenario_doses_per_day_character=="  500") |>
  ggplot() +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1.3,alpha=0.6) +
  geom_line(aes(x = date, y = mean_value, color = scenario_vaccine, linetype = scenario_priority),linewidth=1.3,alpha=0.6) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Bujumbura (50,000 doses)",
       subtitle = "Weekly infections",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = "Vaccine", 
                     values=c("No vaccination" ="black",
                              scenarios_palette_vaccine))+
  scale_linetype_manual(name="Priority", values=linetype_priority)+   
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=14.5),
        legend.direction = "vertical", 
        legend.box = "vertical")+
  guides(color = guide_legend(ncol = 4,order=2)) +
  #guides(shape = guide_legend(ncol = 3)) +
  guides(linetype = guide_legend(ncol = 3,order=1))
line_plot_bujumbura
ggsave(paste0("outputs/rerun/line_plot_bujumbura.png"), line_plot_bujumbura,
       width = 10, height = 5)


################# Burundi EXTENDED DOSES analysis #######################
# Plot max dose scenario
big_df_plots_bu_all <- filter(big_df_bujumbura,
                              scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)])


big_df_cases_bu_all <- big_df_plots_bu_all |>  filter(Category %in% c("cases_inc_cum"),
                                                      TimeStep==max(big_df_bujumbura$TimeStep, na.rm = TRUE))


# percentage of Available full doses administered for each scenario
dose_percent_bujumbura_all <- big_df_cases_bu_all |>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)


# Bar plot of best scenario at each dose 
sensitivity_data_bu_all <- big_df_plots_bu_all |> filter(Category %in% c("cases_inc_cum"), 
                                                         scenario_doses_per_day_character=="  500",
                                                         TimeStep==max(big_df_plots_bu_all$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

best_case_scenario_bu_all <- sensitivity_data_bu_all |> 
  group_by(dose_character) |> 
  top_n(1, mean_averted) |>
  filter(scenario_new=="MVA-BN + fractional + SWs first")|> #filter because of ties
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE), 
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))



sensitivity_data_bu_plot_all <- sensitivity_data_bu_all %>%
  mutate(lower_averted_95_per_dose = case_when(
    lower_averted_95_per_dose < 0 ~ 0.000001,
    .default = lower_averted_95_per_dose
  ),
  mean_averted_per_dose = case_when(
    mean_averted_per_dose < 0 ~ 0.000001,
    .default = mean_averted_per_dose
  ))




##### Excel tables of averted outcomes ####### -------------------------------------------------
# Plot max dose & middle daily dose total scenario

# equateur data 
big_df_data_eq <- filter(big_df_equateur,
                         scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)],
                         Category %in% c("cases_inc_cum"), 
                         TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))


big_df_data_eq_deaths <- filter(big_df_equateur,
                                scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)],
                                Category %in% c("deaths_inc_cum"), 
                                TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))
# sud kicu data
big_df_data_sk <- filter(big_df_sudkivu, scenario_dose < 3000000,
                         scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)],
                         Category %in% c("cases_inc_cum"), 
                         TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))


big_df_data_sk_deaths <- filter(big_df_sudkivu, scenario_dose < 3000000,
                                scenario_labels%in%scenario_names_plots[1:length(scenario_names_plots)],
                                Category %in% c("deaths_inc_cum"), 
                                TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))

# Equateur data 
data_averted_cases_eq<-big_df_data_eq[,c("scenario_vaccine","scenario_priority",
                                         "scenario_doses_per_day_character",
                                         "scenario_dose", "scenario_dose_all",
                                         "doses_per_day_total",
                                         "mean_cumulative_doses",
                                         "lower_95_cumulative_doses",
                                         "upper_95_cumulative_doses",
                                         "mean_averted",
                                         "lower_95_averted","upper_95_averted",
                                         "mean_percent",
                                         "lower_95_percent","upper_95_percent")]

data_averted_deaths_eq<-big_df_data_eq_deaths[,c("scenario_vaccine","scenario_priority",
                                                 "scenario_doses_per_day_character",
                                                 "scenario_dose", "scenario_dose_all",
                                                 "doses_per_day_total",
                                                 "mean_cumulative_doses",
                                                 "lower_95_cumulative_doses",
                                                 "upper_95_cumulative_doses",
                                                 "mean_averted",
                                                 "lower_95_averted","upper_95_averted",
                                                 "mean_percent",
                                                 "lower_95_percent","upper_95_percent")]


colnames(data_averted_cases_eq)<-c( "Vaccine","Priority",
                                    "Strategy doses per day",
                                    "Full doses available", "Doses available",
                                    "Doses per day administered",
                                    "Doses administered mean",
                                    "Doses administered lower 95% CrI",
                                    "Doses administered upper 95% CrI",
                                    "Infections averted mean",
                                    "Infections averted lower 95% CrI","Infections averted upper 95% CrI",
                                    "Percentage of infections averted mean",
                                    "Percentage of infections averted lower 95% CrI",
                                    "Percentage of infections averted upper 95% CrI"
)

colnames(data_averted_deaths_eq)<-c("Vaccine","Priority",
                                    "Strategy doses per day",
                                    "Full doses available", "Doses available",
                                    "Doses per day administered",
                                    "Doses administered mean",
                                    "Doses administered lower 95% CrI",
                                    "Doses administered upper 95% CrI" ,
                                    "Deaths averted mean",
                                    "Deaths averted lower 95% CrI",
                                    "Deaths averted upper 95% CrI",
                                    "Percentage of deaths averted mean",
                                    "Percentage of deaths averted lower 95% CrI",
                                    "Percentage of deaths averted upper 95% CrI")

data_averted_eq<-full_join(data_averted_cases_eq,data_averted_deaths_eq, 
                           by=c("Vaccine","Priority", "Doses available",
                                "Doses per day administered", "Strategy doses per day",
                                "Doses administered mean", "Doses administered lower 95% CrI","Doses administered upper 95% CrI"))|>
  mutate(Province="Equateur")

# Sud Kivu data 
data_averted_cases_sk<-big_df_data_sk[,c("scenario_vaccine","scenario_priority",
                                         "scenario_doses_per_day_character",
                                         "scenario_dose", "scenario_dose_all",
                                         "doses_per_day_total",
                                         "mean_cumulative_doses",
                                         "lower_95_cumulative_doses",
                                         "upper_95_cumulative_doses",
                                         "mean_averted",
                                         "lower_95_averted","upper_95_averted",
                                         "mean_percent",
                                         "lower_95_percent","upper_95_percent")]

data_averted_deaths_sk<-big_df_data_sk_deaths[,c("scenario_vaccine",
                                                 "scenario_priority",
                                                 "scenario_doses_per_day_character",
                                                 "scenario_dose", "scenario_dose_all",
                                                 "doses_per_day_total",
                                                 "mean_cumulative_doses",
                                                 "lower_95_cumulative_doses",
                                                 "upper_95_cumulative_doses",
                                                 "mean_averted",
                                                 "lower_95_averted","upper_95_averted",
                                                 "mean_percent",
                                                 "lower_95_percent","upper_95_percent")]


colnames(data_averted_cases_sk)<-c("Vaccine","Priority",
                                   "Strategy doses per day",
                                   "Full doses available", "Doses available",
                                   "Doses per day administered",
                                   "Doses administered mean",
                                   "Doses administered lower 95% CrI",
                                   "Doses administered upper 95% CrI",
                                   "Infections averted mean",
                                   "Infections averted lower 95% CrI",
                                   "Infections averted upper 95% CrI",
                                   "Percentage of infections averted mean",
                                   "Percentage of infections averted lower 95% CrI",
                                   "Percentage of infections averted upper 95% CrI")

colnames(data_averted_deaths_sk)<-c("Vaccine","Priority",
                                    "Strategy doses per day",
                                    "Full doses available", "Doses available",
                                    "Doses per day administered",
                                    "Doses administered mean",
                                    "Doses administered lower 95% CrI",
                                    "Doses administered upper 95% CrI",
                                    "Deaths averted mean",
                                    "Deaths averted lower 95% CrI",
                                    "Deaths averted upper 95% CrI",
                                    "Percentage of deaths averted mean",
                                    "Percentage of deaths averted lower 95% CrI",
                                    "Percentage of deaths averted upper 95% CrI")

data_averted_sk<-full_join(data_averted_cases_sk, data_averted_deaths_sk, 
                           by=c("Vaccine","Priority", "Doses administered mean","Doses administered lower 95% CrI","Doses administered upper 95% CrI",
                                "Doses per day administered", "Strategy doses per day","Doses available"))|>
  mutate(Province="Sud Kivu", `Doses administered mean`=floor(`Doses administered mean`))


### DRC tables 
data_averted_drc_both <- full_join(data_averted_eq,data_averted_sk)
data_averted_drc <- data_averted_drc_both[,c(dim(data_averted_drc_both)[2],1:(dim(data_averted_drc_both)[2]-1))]
data_averted_drc <- data_averted_drc[order(data_averted_drc$`Doses available`),]
data_averted_drc <- data_averted_drc[order(data_averted_drc$Priority),]
data_averted_drc <- data_averted_drc[order(data_averted_drc$Vaccine),]
data_averted_drc <- data_averted_drc[order(data_averted_drc$Province),]


# Bujumbura data 
data_averted_cases_bu<-big_df_cases_bu_all[,c("scenario_vaccine","scenario_priority",
                                              "scenario_doses_per_day_character",
                                              "scenario_dose", "scenario_dose_all",
                                              "doses_per_day_total",
                                              "mean_cumulative_doses",
                                              "lower_95_cumulative_doses",
                                              "upper_95_cumulative_doses",
                                              "mean_averted",
                                              "lower_95_averted","upper_95_averted",
                                              "mean_percent",
                                              "lower_95_percent","upper_95_percent")]



colnames(data_averted_cases_bu)<-c("Vaccine","Priority",
                                   "Strategy doses per day",
                                   "Full doses available", "Doses available",
                                   "Doses per day administered",
                                   "Doses administered mean",
                                   "Doses administered lower 95% CrI",
                                   "Doses administered upper 95% CrI",
                                   "Infections averted mean",
                                   "Infections averted lower 95% CrI","Infections averted upper 95% CrI",
                                   "Percentage of infections averted mean",
                                   "Percentage of infections averted lower 95% CrI",
                                   "Percentage of infections averted upper 95% CrI"
)
### Burundi tables 
data_averted_cases_bu<-data_averted_cases_bu|>
  mutate(Province="Bujumbura", `Doses administered mean`=floor(`Doses administered mean`))

data_averted_bu <- data_averted_cases_bu[order(data_averted_cases_bu$`Doses available`),c(dim(data_averted_cases_bu)[2],1:(dim(data_averted_cases_bu)[2]-1))]
data_averted_bu <- data_averted_bu[order(data_averted_bu$Priority),]
data_averted_bu <- data_averted_bu[order(data_averted_bu$Vaccine),]
data_averted_bu <- data_averted_bu[order(data_averted_bu$Province),]



# Save the data frame as an Excel file
library(writexl)

write_xlsx(data_averted_drc, "outputs/rerun/DRC_results.xlsx")

write_xlsx(data_averted_bu, "outputs/rerun/Burundi_results.xlsx")



################## DAILY DOSE ANALYSIS ################################
library(RColorBrewer)
#### Data for daily dose analysis 
## EQUATEUR - cases
big_df_plots_eq_DD <- filter(big_df_equateur,
                             scenario_labels%in%scenario_names_plots[1:9])


big_df_cases_eq_DD <- big_df_plots_eq_DD |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))


dose_percent_equateur_DD <- big_df_cases_eq_DD|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

sensitivity_data_eq_DD <- big_df_plots_eq_DD |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_eq$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

## Daily dose heatmap data - Equateur
eq_heatmap_data <- sensitivity_data_eq_DD|>group_by(scenario_dose, scenario_doses_per_day_character,scenario_priority)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))),
         dose_thousands=factor((scenario_dose/1000),
                               levels=unique((scenario_dose/1000))))


## EQUATEUR - deaths
big_df_deaths_eq_DD <- big_df_plots_eq_DD |>  filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))

dose_percent_equateur_deaths_DD <-big_df_deaths_eq_DD |>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

sensitivity_data_eq_deaths_DD <- big_df_plots_eq_DD |> filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_plots_eq$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

## Daily dose heatmap data - Equateur deaths
eq_heatmap_data_deaths <- sensitivity_data_eq_deaths_DD|>group_by(scenario_dose, scenario_doses_per_day_character,scenario_priority)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))),
         dose_thousands=factor((scenario_dose/1000),
                               levels=unique((scenario_dose/1000))))

## SUD KIVU - cases
big_df_plots_sk_DD <- filter(big_df_sudkivu,
                             scenario_labels%in%scenario_names_plots[1:9])|>
  filter(scenario_dose < 3000000)


big_df_cases_sk_DD <- big_df_plots_sk_DD |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))

dose_percent_sudkivu_DD<-big_df_cases_sk_DD|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)


sensitivity_data_sk_DD <- big_df_plots_sk_DD |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_sk$TimeStep, na.rm = TRUE)) |>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))


## Daily dose heatmap data - Sud Kivu
sk_heatmap_data <- sensitivity_data_sk_DD|>group_by(scenario_dose, scenario_doses_per_day_character,scenario_priority)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))),
         dose_thousands=factor((scenario_dose/1000),
                               levels=unique((scenario_dose/1000))))


## SUD KIVU - deaths
big_df_deaths_sk_DD <- big_df_plots_sk_DD |>  filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))

dose_percent_sudkivu_deaths_DD <-big_df_deaths_sk_DD|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

sensitivity_data_sk_deaths_DD <- big_df_plots_sk_DD |> filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_plots_sk$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))

best_case_scenario_sk_deaths_DD <- sensitivity_data_sk_deaths_DD|>group_by(scenario_dose, scenario_doses_per_day_character)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))))

## Daily dose heatmap data - Sud Kivu deaths
sk_heatmap_data_deaths <- sensitivity_data_sk_deaths_DD|>group_by(scenario_dose, scenario_doses_per_day_character,scenario_priority)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))),
         dose_thousands=factor((scenario_dose/1000),
                               levels=unique((scenario_dose/1000))))|>
  mutate(  mean_percent = case_when(
    mean_percent < 0 ~ 0,
    .default = mean_percent
  ))



## BUJUMBURA - cases
big_df_plots_bu_DD <- filter(big_df_bujumbura,
                             scenario_labels%in%scenario_names_plots[1:9])


big_df_cases_bu_DD <- big_df_plots_bu_DD |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_bujumbura$TimeStep, na.rm = TRUE))


dose_percent_bujumbura_DD <-big_df_cases_bu_DD|>
  group_by(Category,scenario_labels, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

sensitivity_data_bu_DD <- big_df_plots_bu_DD |> filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_plots_bu$TimeStep, na.rm = TRUE))|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=format(unique(sort(scenario_dose)), big.mark = ",", scientific = FALSE)))


# Daily dose heatmap data - bujumbura
bu_heatmap_data <- sensitivity_data_bu_DD|>group_by(scenario_dose, scenario_doses_per_day_character,scenario_priority)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))),
         dose_thousands=factor((scenario_dose/1000),
                               levels=unique((scenario_dose/1000))))%>%
  mutate(mean_percent = case_when(
    mean_percent < 0 ~ 0,
    .default = mean_percent
  ))


bu_heatmap_data_extended <- sensitivity_data_bu_DD|>group_by(scenario_dose, scenario_doses_per_day_character,scenario_priority)|>top_n(1, mean_averted)|>
  mutate(dose_character=factor(format(scenario_dose, big.mark = ",", scientific = FALSE),
                               levels=unique((format(scenario_dose, big.mark = ",", scientific = FALSE)))),
         dose_thousands=factor((scenario_dose/1000),
                               levels=unique((scenario_dose/1000))))%>%
  mutate(mean_percent = case_when(
    mean_percent < 0 ~ 0,
    .default = mean_percent
  ))


############################## MAIN FIGURES ####################################
## EQUATEUR - FIGURE 3
## % infections averted
fig3a <- ggplot(data=sensitivity_data_eq, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_percent, colour=scenario_vaccine, shape=scenario_priority), 
             size=2.75, stroke=1.2, alpha=0.75, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_95_percent, 
                    ymax=upper_95_percent, colour=scenario_vaccine),width = 0,
                alpha=0.65, linewidth = 0.9, 
                position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_eq$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Percentage of infections averted", title="Equateur")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete()+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))
fig3a

ggsave(paste0("outputs/rerun/fig3a.png"), fig3a,
       width = 10, height = 6)


# best case scenario for each dose number 
fig3b <- best_case_scenario_eq |> 
  ggplot() +
  geom_bar(stat = "identity", aes(x =  dose_character, y = mean_averted, 
                                  fill =  scenario_vaccine)) + 
 # scale_linetype_manual(name="Priority", values=linetype_priority)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  geom_point(aes(x =  dose_character, y = mean_averted, 
                 shape=scenario_priority), size=2.75, stroke=1.2) +
  geom_errorbar(aes(x =  dose_character, y = mean_averted, #linetype=scenario_priority,
                    ymin = lower_95_averted, ymax = upper_95_averted), width = 0.5, linewidth = 0.9, linewidth = 0.9) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  labs(title = NULL,
       subtitle="Maximum infections averted",
       x = "Available full doses",
       y = NULL)  +
  theme_minimal()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5), 
        legend.position = "bottom")+
  scale_y_continuous(labels = scales::comma)

# daily dose heat map  
equateur_heatmap_maxdoses <- ggplot(eq_heatmap_data, aes(x=dose_thousands, y=scenario_doses_per_day_character, fill= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  labs(x="Available full doses (thousands)", y="Doses per day", title=NULL,
       subtitle="Percentage of infections averted", fill=NULL) +
  scale_fill_viridis_c(name="", option="plasma", limits = c(0,1),
                       labels = scales::percent, begin = 0.1, end = 0.9,
                       breaks = seq(from=0, to=1, by = 0.1))+
  guides(
    fill = guide_colorbar(
      barwidth = unit(15, "cm"))) +
  theme(legend.position = "bottom",legend.text = element_text(size=13),
            axis.text = element_text(size=14),
            strip.text = element_text(size=13),
            plot.subtitle = element_text(size=15),
            plot.title = element_text(size=16),
            axis.title = element_text(size=14.5),
            legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(eq_heatmap_data$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")


# FIGURE 3 
# combine all panels
fig3 <- cowplot::plot_grid(
  fig3a,
  fig3b + theme(legend.position = "none"),
  equateur_heatmap_maxdoses,
  nrow=3,
  rel_heights = c(1.5,1,1.25),align="v")
ggsave(paste0("outputs/rerun/Figure3.png"), fig3,
       width = 15, height = 18, scale=1.5, units="cm")


equateur_heatmap_maxdoses_vaccine <- ggplot(eq_heatmap_data, aes(x=dose_thousands, fill=scenario_vaccine,
                                                                 y=scenario_doses_per_day_character, alpha= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  scale_alpha_continuous(name="Percent averted")+
  labs(x="Available full doses (thousands)", y="Doses per day", title="Equateur",
       subtitle="Percentage of infections averted", fill=NULL) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  theme(legend.position = "bottom",legend.text = element_text(size=13),
            axis.text = element_text(size=14),
            strip.text = element_text(size=13),
            plot.subtitle = element_text(size=15),
            plot.title = element_text(size=16),
            axis.title = element_text(size=14.5),
            legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(eq_heatmap_data$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")


### EQUATEUR DEATHS 

## % deaths averted
fig3a_deaths <- ggplot(data=sensitivity_data_eq_deaths, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_percent, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_95_percent, 
                    ymax=upper_95_percent, colour=scenario_vaccine),width = 0,
                alpha=0.65, linewidth = 0.9, 
                position = position_dodge(width = 0.975))+ 
  theme_minimal()+

  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_eq_deaths$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Percentage of deaths averted", title="Equateur")+
  scale_y_continuous(labels = scales::percent)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))


# best case scenario per dose number 
fig3b_deaths <- best_case_scenario_eq_deaths |> 
  ggplot() +
  geom_bar(stat = "identity", aes(x =  dose_character, y = mean_averted, 
                                  fill =  scenario_vaccine)) + 
  scale_linetype_manual(name="Priority", values=linetype_priority)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  geom_point(aes(x =  dose_character, y = mean_averted, 
                 shape=scenario_priority), size=2.75, stroke=1.2) +
  geom_errorbar(aes(x =  dose_character, y = mean_averted, linetype=scenario_priority,
                    ymin = lower_95_averted, ymax = upper_95_averted), width = 0.5, linewidth = 0.9, linewidth = 0.9) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  labs(title = NULL,
       subtitle="Maximum deaths averted",
       x = "Available full doses",
       y = NULL)  +
  theme_minimal()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5), 
        legend.position = "bottom")+
  scale_y_continuous(labels = scales::comma)


# Daily dose heatmap deaths
equateur_heatmap_deaths <- ggplot(eq_heatmap_data_deaths, aes(x=dose_thousands, y=scenario_doses_per_day_character, fill= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  labs(x="Available full doses (thousands)", y="Doses per day",
       subtitle="Percentage of deaths averted", fill=NULL) +#
  scale_fill_viridis_c(name="", option="plasma", limits = c(0,1),
                       labels = scales::percent, begin = 0.1, end = 0.9,
                       breaks = seq(from=0, to=1, by = 0.1)) +
  guides( 
    fill = guide_colorbar(
      barwidth = unit(15, "cm"))) +
  theme(legend.position = "bottom",legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        strip.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_hline(yintercept =  which(levels(eq_heatmap_data_deaths$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")



# combine 
fig3_deaths <- cowplot::plot_grid(
  fig3a_deaths ,
  fig3b_deaths + theme(legend.position = "none"),
  equateur_heatmap_deaths,
  nrow=3,
  rel_heights = c(1.5,1,1.25),align="v")
ggsave(paste0("outputs/rerun/Figure3_Deaths.png"), fig3_deaths,
       width = 15, height = 18, scale=1.5, units="cm")




equateur_heatmap_maxdoses_vaccine_deaths <- ggplot(eq_heatmap_data_deaths, aes(x=dose_thousands, fill=scenario_vaccine,
                                                                        y=scenario_doses_per_day_character, alpha= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  scale_alpha_continuous(name="Percent averted")+
  labs(x="Available full doses (thousands)", y="Doses per day", title="Equateur",
       subtitle="Percentage of deaths averted", fill=NULL) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  theme(legend.position = "bottom",legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        strip.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(eq_heatmap_data_deaths$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")


### SUD KIVU
## Figure 4

## line plot of percentage of vaccines averted
fig4a<-ggplot(data=sensitivity_data_sk, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_percent, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_sk$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Percentage of infections averted", title="Sud Kivu")+
  scale_y_continuous(labels = scales::percent)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))

ggsave(paste0("outputs/rerun/fig4a.png"), fig4a,
       width = 10, height = 6)

# best case scenario for dose numbers 
fig4b <- best_case_scenario_sk |> 
  ggplot() +
  geom_bar(stat = "identity", aes(x =  dose_character, y = mean_averted, 
                                  fill =  scenario_vaccine)) + 
#  scale_linetype_manual(name="Priority", values=linetype_priority)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  geom_point(aes(x =  dose_character, y = mean_averted, 
                 shape=scenario_priority), size=2.75, stroke=1.2) +
  geom_errorbar(aes(x =  dose_character, y = mean_averted, #linetype=scenario_priority,
                    ymin = lower_95_averted, ymax = upper_95_averted), width = 0.5, linewidth = 0.9) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  labs( title = NULL,
        subtitle="Maximum infections averted",
        x = "Available full doses",
        y = NULL)  +
  theme_minimal()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5), 
        legend.position = "bottom")+
  scale_y_continuous(labels = scales::comma)


# Daily dose heatmap 
sudkivu_heatmap_maxdoses <-ggplot(sk_heatmap_data, aes(x=dose_thousands, y=scenario_doses_per_day_character, fill= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  labs(x="Available full doses (thousands)", y="Doses per day", subtitle="Percentage of infections averted",
       title=NULL, fill=NULL) +
  scale_fill_viridis_c(name="", option="plasma", limits = c(0,1),
                       labels = scales::percent, begin = 0.1, end = 0.9,
                       breaks = seq(from=0, to=1, by = 0.1))+
  guides(
    fill = guide_colorbar(
      barwidth = unit(15, "cm"))) +
  theme(legend.position = "bottom",legend.text = element_text(size=13),
           axis.text = element_text(size=14),
           strip.text = element_text(size=13),
           plot.subtitle = element_text(size=15),
           plot.title = element_text(size=16),
           axis.title = element_text(size=14.5),
           legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(sk_heatmap_data$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")


# combine 
fig4 <- cowplot::plot_grid(
  fig4a,
  fig4b + theme(legend.position = "none"),
  sudkivu_heatmap_maxdoses,
  nrow=3,
  rel_heights = c(1.5,1,1.25),align="v")
ggsave(paste0("outputs/rerun/Figure4.png"), fig4,
       width = 15, height = 18, scale=1.5, units="cm")


sudkivu_heatmap_maxdoses_vaccine <- ggplot(sk_heatmap_data, aes(x=dose_thousands, fill=scenario_vaccine,
                                                                               y=scenario_doses_per_day_character, alpha= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  scale_alpha_continuous(name="Percent averted")+
  labs(x="Available full doses (thousands)", y="Doses per day", title="Sud Kivu",
       subtitle="Percentage of infections averted", fill=NULL) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  theme(legend.position = "bottom",legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        strip.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(sk_heatmap_data$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")

## SUD KIVU DEATHS 
## Figure 4

## line plot of percentage of vaccines averted
fig4a_deaths <-ggplot(data=sensitivity_data_sk_deaths, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_percent, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_sk_deaths$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Percentage of deaths averted", title="Sud Kivu")+
  scale_y_continuous(labels = scales::percent)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))


# best case scenario for dose number deaths
fig4b_deaths <- best_case_scenario_sk_deaths |> 
  ggplot() +
  geom_bar(stat = "identity", aes(x =  dose_character, y = mean_averted, 
                                  fill =  scenario_vaccine)) + 
 # scale_linetype_manual(name="Priority", values=linetype_priority)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  geom_point(aes(x =  dose_character, y = mean_averted, 
                 shape=scenario_priority), size=2.75, stroke=1.2) +
  geom_errorbar(aes(x =  dose_character, y = mean_averted,# linetype=scenario_priority,
                    ymin = lower_95_averted, ymax = upper_95_averted), width = 0.5, linewidth = 0.9) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  labs( title = NULL,
        subtitle="Maximum deaths averted",
        x = "Available full doses",
        y = NULL)  +
  theme_minimal()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom")+
  scale_y_continuous(labels = scales::comma)


# Daily dose heatmap deaths  
sudkivu_heatmap_deaths <- ggplot(sk_heatmap_data_deaths, aes(x=dose_thousands, y=scenario_doses_per_day_character, fill= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  labs(x="Available full doses (thousands)", y="Doses per day", subtitle="Percentage of deaths averted",
       title=NULL, fill=NULL) +
  scale_fill_gradientn(colors = c(brewer.pal(n = 5, name = "Purples")),
                       limits = c(0,1),
                       labels = scales::percent,
                       name = NULL # Fixed range for colour scale
  ) + theme(legend.position = "bottom",legend.text = element_text(size=13),
            axis.text = element_text(size=14),
            strip.text = element_text(size=13),
            plot.subtitle = element_text(size=15),
            plot.title = element_text(size=16),
            axis.title = element_text(size=14.5),
            legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_hline(yintercept =  which(levels(sk_heatmap_data_deaths$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")



# combine 
fig4_deaths <- cowplot::plot_grid(
  fig4a_deaths ,
  fig4b_deaths + theme(legend.position = "none"),
  sudkivu_heatmap_deaths,
  nrow=3,
  rel_heights = c(1.5,1,1.25),align="v")
ggsave(paste0("outputs/rerun/Figure4_deaths.png"), fig4_deaths,
       width = 15, height = 18, scale=1.5, units="cm")



sudkivu_heatmap_maxdoses_vaccine_deaths <- ggplot(sk_heatmap_data_deaths, aes(x=dose_thousands, fill=scenario_vaccine,
                                                                y=scenario_doses_per_day_character, alpha= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  scale_alpha_continuous(name="Percent averted")+
  labs(x="Available full doses (thousands)", y="Doses per day", title="Sud Kivu",
       subtitle="Percentage of deaths averted", fill=NULL) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  theme(legend.position = "bottom",legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        strip.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(sk_heatmap_data_deaths$scenario_doses_per_day_character) == " 5,000"), colour="black", linetype="dashed")


#### BURUNDI 


## Figure 5
## line plot of percentage of vaccines averted
fig5a<-ggplot(data=sensitivity_data_bu, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_percent, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_bu$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Percentage of infections averted", title="Bujumbura")+
  scale_y_continuous(labels = scales::percent)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))

ggsave(paste0("outputs/rerun/fig5a.png"), fig5a,
            width = 10, height = 6)

# best case scenario for dose numbers 
fig5b <- best_case_scenario_bu |> 
  ggplot() +
  geom_bar(stat = "identity", aes(x =  dose_character, y = mean_averted, 
                                  fill =  scenario_vaccine)) + 
 # scale_linetype_manual(name="Priority", values=linetype_priority)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  geom_point(aes(x =  dose_character, y = mean_averted, 
                 shape=scenario_priority), size=2.75, stroke=1.2) +
  geom_errorbar(aes(x =  dose_character, y = mean_averted, #linetype=scenario_priority,
                    ymin = lower_95_averted, ymax = upper_95_averted), width = 0.5, linewidth = 0.9) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  labs( title = NULL,
        subtitle="Maximum infections averted",
        x = "Available full doses",
        y = NULL)  +
  theme_minimal()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5), 
        legend.position = "bottom")+
  scale_y_continuous(labels = scales::comma)


# Daily dose heatmap
burundi_heatmap_maxdoses <- ggplot(bu_heatmap_data|>filter(scenario_dose<=50000), aes(x=dose_thousands, y=scenario_doses_per_day_character, fill= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  labs(x="Available full doses (thousands)", y="Doses per day", subtitle="Percentage of infections averted",
       title=NULL, fill=NULL) +
  # scale_fill_gradientn(colors = c(brewer.pal(n = 5, name = "Purples")),
  #                      limits = c(0,1),
  #                      labels = scales::percent,
  #                      name = NULL # Fixed range for colour scale
  # ) +
  scale_fill_viridis_c(name="", option="plasma", limits = c(0,1),
                       labels = scales::percent, begin = 0.1, end = 0.9,
                       breaks = seq(from=0, to=1, by = 0.1))+
  guides(
    fill = guide_colorbar(
      barwidth = unit(15, "cm"))) +
  theme(legend.position = "bottom",legend.text = element_text(size=13),
            axis.text = element_text(size=14),
            plot.subtitle = element_text(size=15),
            plot.title = element_text(size=16),
            strip.text = element_text(size=13),
            axis.title = element_text(size=14.5),
            legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(bu_heatmap_data$scenario_doses_per_day_character) == "  500"), colour="black", linetype="dashed")



# combine 
fig5 <- cowplot::plot_grid(
  fig5a,
  fig5b + theme(legend.position = "none"),
  burundi_heatmap_maxdoses,
  nrow=3,
  rel_heights = c(1.5,1,1.25),align="v")
ggsave(paste0("outputs/rerun/Figure5.png"), fig5,
       width = 15, height = 18, scale=1.5, units="cm")


bu_heatmap_maxdoses_vaccine <- ggplot(bu_heatmap_data, aes(x=dose_thousands, fill=scenario_vaccine,
                                                                y=scenario_doses_per_day_character, alpha= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  scale_alpha_continuous(name="Percent averted")+
  labs(x="Available full doses (thousands)", y="Doses per day", title="Bujumbura",
       subtitle="Percentage of infections averted", fill=NULL) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  theme(legend.position = "bottom",legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        strip.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(bu_heatmap_data$scenario_doses_per_day_character) == "  500"), colour="black", linetype="dashed")


### BURUNDI EXTENDED 
## Figure 5

## line plot of percentage of vaccines averted
fig5a_ext <-ggplot(data=sensitivity_data_bu_all, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_percent, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_bu_all$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Percentage of infections averted", title="Bujumbura")+
  scale_y_continuous(labels = scales::percent)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))

ggsave(paste0("outputs/rerun/fig5a_ext.png"), fig5a_ext,
       width = 12, height = 6)


# best case scenario for dose numbers
fig5b_ext <- best_case_scenario_bu_all |> 
  ggplot() +
  geom_bar(stat = "identity", aes(x =  dose_character, y = mean_averted, 
                                  fill =  scenario_vaccine)) + 
 # scale_linetype_manual(name="Priority", values=linetype_priority)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  geom_point(aes(x =  dose_character, y = mean_averted, 
                 shape=scenario_priority), size=2.75, stroke=1.2) +
  geom_errorbar(aes(x =  dose_character, y = mean_averted,# linetype=scenario_priority,
                    ymin = lower_95_averted, ymax = upper_95_averted), width = 0.5, linewidth = 0.9) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  labs( title = NULL,
        subtitle="Maximum infections averted",
        x = "Available full doses",
        y = NULL)  +
  theme_minimal()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical")+
  scale_y_continuous(labels = scales::comma)


# daily dose heatmap 
burundi_heatmap_maxdoses_extended <- ggplot(bu_heatmap_data_extended, aes(x=dose_thousands, y=scenario_doses_per_day_character, fill= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+ 
  #guides(fill=guide_legend(keywidth=2)) +
  labs(x="Available full doses (thousands)", y="Doses per day", subtitle="Percentage of infections averted",
       title=NULL, fill=NULL) +
  scale_fill_viridis_c(name="", option="plasma", limits = c(0,1),
                       labels = scales::percent, begin = 0.1, end = 0.9,
                       breaks = seq(from=0, to=1, by = 0.1))+
  guides(
    fill = guide_colorbar(
      barwidth = unit(15, "cm"))) +
  theme(legend.position = "bottom",legend.text = element_text(size=13),
            axis.text = element_text(size=14),
            strip.text = element_text(size=13),
            plot.subtitle = element_text(size=15),
            plot.title = element_text(size=16),
            axis.title = element_text(size=14.5),
            legend.key.width = unit(1, "cm")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_hline(yintercept =  which(levels(bu_heatmap_data_extended$scenario_doses_per_day_character) == "  500"), colour="black", linetype="dashed")



# combine 
fig5_ext <- cowplot::plot_grid(
  fig5a_ext ,
  fig5b_ext + theme(legend.position = "none"),
  burundi_heatmap_maxdoses_extended,
  nrow=3,
  rel_heights = c(1.5,1,1.25),align="v")
ggsave(paste0("outputs/rerun/Figure5_extended.png"), fig5_ext,
       width = 17, height = 18, scale=1.5, units="cm")


bu_heatmap_maxdoses_vaccine_ext <- ggplot(bu_heatmap_data_extended, aes(x=dose_thousands, fill=scenario_vaccine,
                                                           y=scenario_doses_per_day_character, alpha= mean_percent)) + 
  geom_tile()+
  theme_minimal()+
  facet_wrap(~scenario_priority)+
  scale_alpha_continuous(name="Percent averted")+
  labs(x="Available full doses (thousands)", y="Doses per day", title="Bujumbura",
       subtitle="Percentage of infections averted", fill=NULL) +
  scale_fill_manual(name="Vaccine", values=scenarios_palette_vaccine)+
  theme(legend.position = "bottom",legend.text = element_text(size=13),
        axis.text = element_text(size=14),
        strip.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.key.width = unit(1, "cm")) + scale_x_discrete(guide = guide_axis(angle = 90))+
  geom_hline(yintercept =  which(levels(bu_heatmap_data_extended$scenario_doses_per_day_character) == "  500"), colour="black", linetype="dashed")



### APPENDIX PLOT - AVERTED PER DOSE  ###
##### INFECTIONS AVERTED PER DOSE PLOTS 

# EQUATEUR - infections averted per dose 
equateur_averted_per_dose_a <-
  ggplot(data=sensitivity_data_eq_plot, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_averted_per_dose, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_averted_95_per_dose, ymax=upper_averted_95_per_dose, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_eq_plot$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x=NULL, y = NULL, subtitle="Equateur",title="Infections averted per dose (log scale)")+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))+
  scale_y_log10(labels = scales::comma)



# SUD KIVU - infections averted per dose 
sudkivu_averted_per_dose_b<- ggplot(data=sensitivity_data_sk_plot, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_averted_per_dose, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_averted_95_per_dose, ymax=upper_averted_95_per_dose, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_sk_plot$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x=NULL, y=NULL,
       subtitle="Sud Kivu",title=NULL)+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))+
  scale_y_log10(labels = scales::comma)



# BURUNDI - infections averted per dose 
burundi_averted_per_dose_d <- ggplot(data=sensitivity_data_bu_plot_all, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_averted_per_dose, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_averted_95_per_dose, ymax=upper_averted_95_per_dose, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_bu_plot_all$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL,subtitle="Bujumbura")+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))+
  scale_y_log10(labels = scales::comma)



# BURUNDI - averted per dose (extended)
burundi_averted_per_dose_c<- ggplot(data=sensitivity_data_bu_plot, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_averted_per_dose, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_averted_95_per_dose, ymax=upper_averted_95_per_dose, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_bu_plot$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL, subtitle = "Bujumbura")+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4,order=2)) +
  guides(shape = guide_legend(ncol = 3,order=1))+
  scale_y_log10(labels = scales::comma)

# combine cases averted plots per dose
fig_averted_per_dose <- cowplot::plot_grid(
  equateur_averted_per_dose_a + theme(legend.position = "none"),
  sudkivu_averted_per_dose_b + theme(legend.position = "none") ,
  burundi_averted_per_dose_d + theme(legend.position = "bottom"),
  nrow=3,
  rel_heights = c(1.1,1,1.8),align="v")
ggsave(paste0("outputs/rerun/Figure_cases_averted_per_dose.png"), fig_averted_per_dose,
       width = 21, height = 20, scale=1, units="cm")



##### DEATHS AVERTED PER DOSE PLOT 

equateur_deaths_averted_per_dose_a<-
  ggplot(data=sensitivity_data_eq_plot_deaths, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_averted_per_dose, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_averted_95_per_dose, ymax=upper_averted_95_per_dose, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_eq_plot_deaths$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x=NULL, y=NULL,
       subtitle="Equateur",title="Deaths averted per dose (log scale)")+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))+
  scale_y_log10(labels = scales::comma)


# SUD KIVU - deaths averted per dose 
sudkivu_deaths_averted_per_dose_b <- ggplot(data=sensitivity_data_sk_plot_deaths, aes(group=scenario_labels))+
  geom_point(aes(x=dose_character, y=mean_averted_per_dose, colour=scenario_vaccine, shape=scenario_priority), size=2.75, stroke=1.2, alpha=0.65, position = position_dodge(width = 0.975))+
  scale_color_manual(name = "Vaccine", values= scenarios_palette_vaccine)+
  geom_errorbar(aes(x=dose_character, ymin=lower_averted_95_per_dose, ymax=upper_averted_95_per_dose, colour=scenario_vaccine),width = 0, alpha=0.65, linewidth = 0.9, position = position_dodge(width = 0.975))+ 
  theme_minimal()+
  geom_vline(xintercept = seq(1.5, length(unique(sensitivity_data_sk_plot_deaths$dose_character)) - 0.5, 1),
             linetype = "dashed", color = "grey92")+
  labs(x="Available full doses", y=NULL, title=NULL, subtitle="Sud Kivu")+
  scale_shape_manual(name="Priority", values=shape_priority)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=15),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14.5),
        legend.position = "bottom", 
        legend.direction = "vertical",
        legend.box = "vertical") +
  guides(color = guide_legend(ncol = 4)) +
  guides(shape = guide_legend(ncol = 3))+
  scale_y_log10(labels = scales::comma,n.breaks=4)


# combine deaths averted plots per dose

fig_deaths_averted_per_dose <- cowplot::plot_grid(
  equateur_deaths_averted_per_dose_a + theme(legend.position = "none"),
  sudkivu_deaths_averted_per_dose_b + theme(legend.position = "bottom"),
  nrow=2,
  rel_heights = c(1.1,1.8),align="v")
ggsave(paste0("outputs/rerun/Figure_deaths_averted_per_dose.png"), fig_deaths_averted_per_dose,
       width = 21, height = 16, scale=1, units="cm")


## SAVE HEATMAPS

heatmaps_vaccines <- cowplot::plot_grid(
  equateur_heatmap_maxdoses_vaccine,
  sudkivu_heatmap_maxdoses_vaccine,
  bu_heatmap_maxdoses_vaccine,
  nrow=3,
  rel_heights = c(1,1,1),align="v")
ggsave(paste0("outputs/rerun/heatmaps_infections_vaccines.png"), heatmaps_vaccines,
       width = 20, height = 20, scale=1.5, units="cm")

heatmaps <- cowplot::plot_grid(
  equateur_heatmap_maxdoses+labs(title="Equateur"),
  sudkivu_heatmap_maxdoses+labs(title="Sud Kivu"),
  burundi_heatmap_maxdoses+labs(title="Bujumbura"),
  nrow=3,
  rel_heights = c(1,1,1),align="v")
ggsave(paste0("outputs/rerun/heatmaps_infections.png"), heatmaps,
       width = 20, height = 20, scale=1.5, units="cm")



## Sup plot for fit assumption ------------------
fig_averted <- cowplot::plot_grid(fig3a + theme(legend.position = "none")+
                                    labs(title="Percentage of infections averted",
                                         subtitle="Equateur"), 
                     fig4a + theme(legend.position = "none")+
                       labs(title=NULL,
                            subtitle="Sud Kivu"),
                     fig5a+
                       labs(title=NULL,
                            subtitle="Bujumbura"),
                     nrow=3, labels=c("D", "E", "F"),
                     rel_heights = c(1.1,1.1,2))


fig_averted
save(fig_averted, file="outputs/fig_averted.RData")
ggsave("outputs/Figure_averted.png", fig_averted,
       width = 20, height = 15, scale = 1.5,
       units = "cm")

### NEED TO FUN PAPER.FIGURES.R FIRST
# combine figures 
load("Z:/Alba/mpox-drc-outputs/src/paper_figures/outputs/fig_fits.RData")

fig_assumption <- cowplot::plot_grid(fig_fits,
                                  fig_averted,
                                  nrow=2, labels=c("", ""),
                                  rel_heights = c(1,3))


fig_assumption
save(fig_assumption, file="outputs/fig_assumption.RData")
ggsave(paste("outputs/Figure_",assumption,"_assumption.png"), fig_assumption,
       width = 20, height = 20, scale = 1.5,
       units = "cm")

