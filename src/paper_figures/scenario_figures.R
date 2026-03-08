
orderly2::orderly_parameters(short_run = FALSE,
                             fit_by_age = TRUE,
                             deterministic = FALSE,
                             mixing_matrix = "Zimbabwe",
                             fit_KPs = TRUE, 
                             assumptions = "standard",
                             vaccines_onset = "start")

# Load in all dependencies 
orderly2::orderly_dependency(
  name = "pmcmc",
  query = "latest(parameter:region=='sudkivu' && parameter:assumptions == 'standard')",
  files = c("inputs/samples_sudkivu.rds" = "outputs/samples.rds",
            "inputs/fitting_data_sudkivu.rds" = "outputs/fitting_data.rds",
            "outputs/traceplots_sudkivu.png" = "outputs/traceplots.png"))
orderly2::orderly_dependency(
  name = "pmcmc",
  query = "latest(parameter:region=='equateur'  &&  parameter:assumptions == 'standard')",
  files = c("inputs/samples_equateur.rds" = "outputs/samples.rds",
            "inputs/fitting_data_equateur.rds" = "outputs/fitting_data.rds",
            "outputs/traceplots_equateur.png" = "outputs/traceplots.png"))
orderly2::orderly_dependency(name = "run_postprocess",
                             "latest(parameter:region =='sudkivu' && 
                                     parameter:deterministic == FALSE && 
                                     parameter:use_both_fit == FALSE && 
                                     parameter:fit_by_age == TRUE && 
                                     parameter:fit_KPs == TRUE &&
                                     parameter:short_run == FALSE && 
                                     parameter:vaccines_onset == 'start')",
                             files=c("inputs/runs_postprocess_sudkivu.RDS" = "outputs/runs_postprocess.RDS"))
orderly2::orderly_dependency(name = "run_postprocess",
                             "latest(parameter:region =='equateur' && 
                                     parameter:deterministic == FALSE && 
                                     parameter:use_both_fit == FALSE && 
                                     parameter:fit_by_age == TRUE && 
                                     parameter:fit_KPs == TRUE &&
                                     parameter:short_run == FALSE && 
                                     parameter:vaccines_onset == 'start')",
                             files=c("inputs/runs_postprocess_equateur.RDS" = "outputs/runs_postprocess.RDS"))
orderly2::orderly_dependency(name = "run_postprocess_burundi",
                             "latest(parameter:deterministic == this:deterministic && 
                             parameter:short_run == this:short_run)",
                             files=c("inputs/runs_postprocess_burundi.RDS" = "outputs/runs_postprocess.RDS"))
orderly2::orderly_dependency(name = "pmcmc_burundi",
                             "latest(parameter:deterministic == this:deterministic && parameter:short_run == this:short_run)",
                             files =  c("inputs/samples.rds" = "outputs/samples.rds",
                                        "inputs/fitting_data.rds" = "outputs/fitting_data.rds"))
# Populations of provinces
province_pop<-mpoxseir::parameters_demographic(region="sudkivu")$province_pop

# Scenario colour palette
scenarios_palette_slides <- c("No vaccination"="black",
                              "A: LC16m8"="#5B859EFF", 
                              "B: A + MVA-BN >12y"= "#75884BFF",
                              "C: MVA-BN"= "#D48F90FF" ,
                              "D: C + 1st dose only"="#AB84A5FF",
                              "E: D + <12y first"="#B38711FF",
                              "F: D + key pops first"="#800000")#"#429264")#"DF8D71FF")
# Region colour pallette
cols_region <- c(MetBrewer::met.brewer("Cassatt1", 7)[c(5,1)],MetBrewer::met.brewer("Demuth", 3)[3], MetBrewer::met.brewer("Cassatt1", 7)[3]) 
names(cols_region) <- c("equateur", "sudkivu", "both", "burundi")
nms_region <- c(equateur = "Equateur", sudkivu = "Sud Kivu", burundi = "Burundi")

# Dates 
start_date<-if(vaccines_onset=="end"){
  mpoxseir::mpoxseir_date_as_date(770) # needs updating if fits changed
}else{
  mpoxseir::mpoxseir_date_as_date(371) #+ 14
} 

library(tidyverse)

# Figure 3: Equateur ---------------------------------------------------------------------
region <-"equateur"

# load in post process
big_df_equateur <- readRDS("inputs/runs_postprocess_equateur.RDS")|> 
  mutate(Date = start_date + lubridate::days(TimeStep * 7)) |> 
  mutate(scenario_combo=scenario_name, 
         scenario_name=sub(" \\(.*", "", scenario_new), 
         scenario_dose=total_doses_children+total_doses_adults)|>
  filter(scenario_num>49)

# load in fitting data
fitting_data_equateur <-  readRDS("inputs/fitting_data_equateur.RDS")

saveRDS(big_df_equateur, 
        file=paste0("outputs/big_df_", region,"_",vaccines_onset, ".rds" ))


# all scenario names A-F
scenario_names <- names(scenarios_palette_slides)

# Plot max dose scenario
big_df_plots_eq <- filter(big_df_equateur, scenario_num%in%c((max(big_df_equateur$scenario_num)-5):max(big_df_equateur$scenario_num)), ignore.case = TRUE)

big_df_eq <- big_df_plots_eq |>
  mutate(scenario_name = factor(scenario_name,
                                levels = c(scenario_names[1], sort(scenario_names[-1]))))


line_plot_equateur <- big_df_eq |> 
  filter(Category =="cases_inc") |>
  filter(scenario_new %in% sort(unique(big_df_eq$scenario_new))[1:(6)]) |>
  mutate(Vaccination=factor(scenario_name, levels=sort(unique(big_df_eq$scenario_name))[1:(6)])) |>
  ggplot() +
  geom_bar(data = fitting_data_equateur$data |>
             filter(region==region) |>
             mutate(date = mpoxseir::mpoxseir_date_as_date(date)),
           aes(x = date, y = cases), colour = "grey50", fill = cols_region["equateur"],
           stat = "identity") +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1,alpha=0.7) +
  geom_line(aes(x = date, y = mean_value, color = Vaccination),linewidth=1,alpha=0.7) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Weekly infections",
       subtitle ="Maximum doses",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = NULL, 
                     values=c("No vaccination" ="black",
                              scenarios_palette_slides[c(1:7)]))+
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=13),
        plot.title = element_text(size=14),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13.5))
line_plot_equateur


###  Bar charts of cases averted

big_df_cases_eq <- big_df_equateur |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))

# percentage of available doses administered for each scenario
dose_percent_equateur<-big_df_cases_eq|>
  group_by(Category, scenario_name, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

cases_averted_per_dose_equatuer<-big_df_cases_eq |> filter( scenario_num%in%c(62:67))|>
  ggplot(aes(x =  substring(scenario_name,1,1), y = mean_averted_per_dose,
             fill =  substring(scenario_name,1,1))) +
  geom_bar(stat = "identity") +  
  geom_errorbar(aes( x = substring(scenario_name,1,1),ymin = lower_averted_95_per_dose, ymax = upper_averted_95_per_dose), width = 0.2) +
  scale_fill_manual(name=NULL, values=c("A"=scenarios_palette_slides[[2]], 
                                        "B"=scenarios_palette_slides[[3]], 
                                        "C"=scenarios_palette_slides[[4]],
                                        "D"=scenarios_palette_slides[[5]],
                                        "E"=scenarios_palette_slides[[6]], 
                                        "F"=scenarios_palette_slides[[7]]))+
  labs( title = "Infections averted per dose",
        x = NULL,
        y = NULL)  +
  theme(axis.text = element_text(size=13),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme_minimal(base_size = 13) + guides(fill = "none") +
  scale_y_continuous(labels = scales::comma)



cases_averted_equateur <- big_df_cases_eq |> filter( scenario_num%in%c(62:67))|>
  ggplot(aes(x =  substring(scenario_name,1,1), y = mean_averted,
             fill =  substring(scenario_name,1,1))) +
  geom_bar(stat = "identity") +  
  geom_errorbar(aes( x = substring(scenario_name,1,1),ymin = lower_95_averted, ymax = upper_95_averted), width = 0.2) +
  scale_fill_manual(name=NULL, values=c("A"=scenarios_palette_slides[[2]], 
                                        "B"=scenarios_palette_slides[[3]], 
                                        "C"=scenarios_palette_slides[[4]],
                                        "D"=scenarios_palette_slides[[5]],
                                        "E"=scenarios_palette_slides[[6]], 
                                        "F"=scenarios_palette_slides[[7]]))+
  labs( title = "Infections averted",
        x = NULL,
        y = NULL)  +
  theme(axis.text = element_text(size=13),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme_minimal(base_size = 13) + guides(fill = "none") +
  scale_y_continuous(labels = scales::comma)



## CASES AVERTED SESNSITIVITY PLOTS

cases_averted_sensitivity_eq<-ggplot(data=big_df_cases_eq)+
  geom_point(aes(x=mean_cumulative_doses/province_pop$equateur, y=mean_percent, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  geom_line(aes(x=mean_cumulative_doses/province_pop$equateur, y=mean_percent, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  geom_errorbar(aes(x=mean_cumulative_doses/province_pop$equateur, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_name),width = 0)+ #,position = position_dodge(width = 0.5
  theme_minimal()+labs(x="Population vaccine coverage", title="Equateur", subtitle="Percentage of infections averted", y=NULL)+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent, breaks=waiver(), n.breaks=6)+
  theme(legend.position = "none",
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))


cases_averted_per_dose_sensitivity_eq<-ggplot(data=big_df_cases_eq)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", y=NULL, subtitle="Number of infections averted per dose")+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=6)+
  theme(legend.position = "none",
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))
        


# FIGURE PLOT

fig3<-cowplot::plot_grid(cases_averted_sensitivity_eq,cases_averted_per_dose_sensitivity_eq,line_plot_equateur, align = "v",
                         rel_heights=c(1.25,1,1.5),labels=c("A","B","C"), ncol=1)
fig3

ggsave(paste0("outputs/Figure_3.png"), fig3,
       width = 12.5, height = 15, scale=1.5, units="cm")


# Figure 4: Sud Kivu ---------------------------------------------------------------------
region <-"sudkivu"

# load in post process
big_df_sudkivu <- readRDS("inputs/runs_postprocess_sudkivu.RDS")|> 
  mutate(Date = start_date + lubridate::days((TimeStep-1) * 7)) |> 
  mutate(scenario_combo=scenario_name, scenario_name=sub(" \\(.*", "", scenario_new), 
         scenario_dose=total_doses_children+total_doses_adults)|>
  filter(scenario_num>49)

# load in fitting data 
fitting_data_sudkivu <-  readRDS("inputs/fitting_data_sudkivu.RDS")


saveRDS(big_df_sudkivu, 
        file=paste0("outputs/big_df_", region,"_",vaccines_onset, ".rds" ))


# all scenario names A-F
scenario_names <- names(scenarios_palette_slides)

# Plot max dose scenario
big_df_plots_sk <- filter(big_df_sudkivu, scenario_num%in%c((max(big_df_sudkivu$scenario_num)-5):max(big_df_sudkivu$scenario_num)), ignore.case = TRUE)

big_df_sk <- big_df_plots_sk |>
  mutate(scenario_name = factor(scenario_name,
                                levels = c(scenario_names[1], sort(scenario_names[-1]))))


line_plot_sudkivu <- big_df_sk |> 
  filter(Category =="cases_inc") |>
  filter(scenario_new %in% sort(unique(big_df_sk$scenario_new))[1:(6)]) |>
  mutate(Vaccination=factor(scenario_name, levels=sort(unique(big_df_sk$scenario_name))[1:(6)])) |>
  ggplot() +
  geom_bar(data = fitting_data_sudkivu$data |>
             filter(region==region) |>
             mutate(date = mpoxseir::mpoxseir_date_as_date(date)),
           aes(x = date, y = cases), colour = "grey50", fill = cols_region["sudkivu"],
           stat = "identity") +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1,alpha=0.7) +
  geom_line(aes(x = date, y = mean_value, color = Vaccination),linewidth=1,alpha=0.7) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Weekly infections",
       subtitle = "Maximum doses",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = NULL, 
                     values=c("No vaccination" ="black",
                              scenarios_palette_slides[c(1:7)]))+
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=13),
        plot.title = element_text(size=14),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13.5))
line_plot_sudkivu


###  Bar charts of cases averted

big_df_cases_sk <- big_df_sudkivu |>  filter(Category %in% c("cases_inc_cum"), TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))


# percentage of available doses administered for each scenario
dose_percent_sudkivu<-big_df_cases_sk|>
  group_by(Category,scenario_name, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

cases_averted_per_dose_sudkivu<-big_df_cases_sk |> filter( scenario_num%in%c(62:67))|>
  ggplot(aes(x =  substring(scenario_name,1,1), y = mean_averted_per_dose,
             fill =  substring(scenario_name,1,1))) +
  geom_bar(stat = "identity") +  
  geom_errorbar(aes( x = substring(scenario_name,1,1),ymin = lower_averted_95_per_dose, ymax = upper_averted_95_per_dose), width = 0.2) +
  scale_fill_manual(name=NULL, values=c("A"=scenarios_palette_slides[[2]], 
                                        "B"=scenarios_palette_slides[[3]], 
                                        "C"=scenarios_palette_slides[[4]],
                                        "D"=scenarios_palette_slides[[5]],
                                        "E"=scenarios_palette_slides[[6]], 
                                        "F"=scenarios_palette_slides[[7]]))+
  labs( title = "Infections averted per dose",
        x = NULL,
        y = NULL)  +
  theme(axis.text = element_text(size=13),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme_minimal(base_size = 13) + guides(fill = "none") +
  scale_y_continuous(labels = scales::comma)


cases_averted_sudkivu <- big_df_cases_sk |> filter( scenario_num%in%c(62:67))|>
  ggplot(aes(x =  substring(scenario_name,1,1), y = mean_averted,
             fill =  substring(scenario_name,1,1))) +
  geom_bar(stat = "identity") +  
  geom_errorbar(aes( x = substring(scenario_name,1,1),ymin = lower_95_averted, ymax = upper_95_averted), width = 0.2) +
  scale_fill_manual(name=NULL, values=c("A"=scenarios_palette_slides[[2]], 
                                        "B"=scenarios_palette_slides[[3]], 
                                        "C"=scenarios_palette_slides[[4]],
                                        "D"=scenarios_palette_slides[[5]],
                                        "E"=scenarios_palette_slides[[6]], 
                                        "F"=scenarios_palette_slides[[7]]))+
  labs( title = "Infections averted",
        x = NULL,
        y = NULL)  +
  theme(axis.text = element_text(size=13),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme_minimal(base_size = 13) + guides(fill = "none") +
  scale_y_continuous(labels = scales::comma)

## CASES AVERTED SESNSITIVITY PLOTS

cases_averted_sensitivity_sk<-ggplot(data=big_df_cases_sk)+
  geom_point(aes(x=mean_cumulative_doses/province_pop$sudkivu, y=mean_percent, colour=scenario_name), alpha=0.7)+#,position = position_dodge(width = 0.5)
  geom_line(aes(x=mean_cumulative_doses/province_pop$sudkivu, y=mean_percent, colour=scenario_name), alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  geom_errorbar(aes(x=mean_cumulative_doses/province_pop$sudkivu, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_name),width = 0)+ #,position = position_dodge(width = 0.5
  theme_minimal()+labs(x="Population vaccine coverage", y=NULL,subtitle="Percentage of infections averted", title="Sud Kivu")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent, breaks=waiver(), n.breaks=6)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))
cases_averted_sensitivity_sk

cases_averted_per_dose_sensitivity_sk<-ggplot(data=big_df_cases_sk)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name), alpha=0.7)+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", subtitle="Number of infections averted per dose", y=NULL)+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=4)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))
cases_averted_per_dose_sensitivity_sk


cases_averted_per_dose_sensitivity_sk_log<-ggplot(data=big_df_cases_sk)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name), alpha=0.7)+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", subtitle="Number of infections averted per dose (log scale)", y=NULL)+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=4)+
  scale_y_log10()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))
cases_averted_per_dose_sensitivity_sk_log

# FIGURE PLOT

fig4<-cowplot::plot_grid(cases_averted_sensitivity_sk,cases_averted_per_dose_sensitivity_sk_log,line_plot_sudkivu, align = "v",
                         rel_heights=c(1,1,1.5),labels=c("A","B","C"), ncol=1)

fig4

ggsave(paste0("outputs/Figure_4.png"), fig4,
       width = 12.5, height = 15, scale=1.5, units="cm")



## Figure 5: Burundi plots -------------------------------------------------------------------

region <- "bujumbura"
## Fitting outputs
fitting_data_bujumbura <- readRDS("inputs/fitting_data.rds")
samples <- readRDS("inputs/samples.rds")

# Data frame of all scenarios.
start_date<-if(vaccines_onset=="start"){
  mpoxseir::mpoxseir_date_as_date(595) 
}else{
  mpoxseir::mpoxseir_date_as_date(735)
  
}
big_df_all <- readRDS("inputs/runs_postprocess_burundi.RDS") |> 
  mutate(scenario_combo=scenario_name, scenario_name=sub(" \\(.*", "", scenario_new), scenario_dose=total_doses_children+total_doses_adults)

saveRDS(big_df_all, 
        file=paste0("outputs/big_df_", region,"_start", ".rds" ))


big_df_burundi_plots <- filter(big_df_all,
                               scenario_num%in%c(44:49,NA))

big_df_burundi <- big_df_burundi_plots |>
  mutate(scenario_name = factor(scenario_name,
                                levels = c(scenario_names[1], 
                                           sort(scenario_names[-1]))))


line_plots_burundi <- big_df_burundi |> 
  filter(Category =="cases_inc") |>
  mutate(Vaccination=scenario_name) |>
  ggplot() +
  geom_bar(data = fitting_data_bujumbura$data |>
             mutate(date = mpoxseir::mpoxseir_date_as_date(date)),
           aes(x = date, y = cases), colour = "grey50", fill=cols_region["burundi"],
           stat = "identity") +
  geom_line(aes(x = date, y = mean_value, color = Vaccination),linewidth=1,alpha=0.7) +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1,alpha=0.7) +
  labs(title =  "Weekly infections",
       subtitle = "Maximum doses",
       x = NULL,
       y = NULL) + 
  scale_x_date(date_breaks="4 months", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=13),
        plot.title = element_text(size=14),
        axis.title = element_text(size=13.5))+
  guides(colour=guide_legend(ncol=3))
line_plots_burundi


###  Bar charts of cases averted

big_df_cases_bu <- big_df_all |>  filter(Category %in% c("cases_inc_cum"), 
                                         TimeStep==max(big_df_all$TimeStep, na.rm = TRUE),
                                         scenario_num<50)

# percentage of available doses administered for each scenario
dose_percent_burundi<-big_df_cases_bu|>
  group_by(Category,scenario_name, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

cases_averted_per_dose<-big_df_cases_bu |> filter(scenario_num%in%c(44:49))|>
  ggplot(aes(x =  substring(scenario_name,1,1), y = mean_averted_per_dose,
             fill =  substring(scenario_name,1,1))) +
  geom_bar(stat = "identity") +  
  geom_errorbar(aes( x = substring(scenario_name,1,1),ymin = lower_95_averted_per_dose, ymax = upper_95_averted_per_dose), width = 0.2) +
  scale_fill_manual(name=NULL, values=c("A"=scenarios_palette_slides[[2]], 
                                        "B"=scenarios_palette_slides[[3]], 
                                        "C"=scenarios_palette_slides[[4]],
                                        "D"=scenarios_palette_slides[[5]],
                                        "E"=scenarios_palette_slides[[6]], 
                                        "F"=scenarios_palette_slides[[7]]))+
  labs( title = "Infections averted per dose",
        x = NULL,
        y = NULL)  +
  theme(axis.text = element_text(size=13),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme_minimal(base_size = 13) + guides(fill = "none") +
  scale_y_continuous(labels = scales::comma)


cases_averted <- big_df_cases_bu |> filter(scenario_num%in%c(44:49))|>
  ggplot(aes(x =  substring(scenario_name,1,1), y = mean_averted,
             fill =  substring(scenario_name,1,1))) +
  geom_bar(stat = "identity") +  
  geom_errorbar(aes( x = substring(scenario_name,1,1),ymin = lower_95_averted, ymax = upper_95_averted), width = 0.2) +
  scale_fill_manual(name=NULL, values=c("A"=scenarios_palette_slides[[2]], 
                                        "B"=scenarios_palette_slides[[3]], 
                                        "C"=scenarios_palette_slides[[4]],
                                        "D"=scenarios_palette_slides[[5]],
                                        "E"=scenarios_palette_slides[[6]], 
                                        "F"=scenarios_palette_slides[[7]]))+
  labs( title = "Infections averted",
        x = NULL,
        y = NULL)  +
  theme(axis.text = element_text(size=13),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme_minimal(base_size = 13) + guides(fill = "none") +
  scale_y_continuous(labels = scales::comma)

## CASES AVERTED SESNSITIVITY PLOTS

cases_averted_sensitivity<-ggplot(data=big_df_cases_bu)+
  geom_point(aes(x=mean_cumulative_doses/province_pop$bujumbura, y=mean_percent, colour=scenario_name))+#,position = position_dodge(width = 0.5)
  geom_line(aes(x=mean_cumulative_doses/province_pop$bujumbura, y=mean_percent, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  geom_errorbar(aes(x=mean_cumulative_doses/province_pop$bujumbura, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_name),width = 0)+ #,position = position_dodge(width = 0.5
  theme_minimal()+labs(x="Population vaccine coverage", y=NULL, title="Bujumbura" ,subtitle ="Percentage of infections averted")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent, breaks=waiver(), n.breaks=6)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none", plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"))

cases_averted_per_dose_sensitivity<-ggplot(data=big_df_cases_bu)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name))+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", y=NULL, subtitle="Number of infections averted per dose")+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=6)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))

cases_averted_per_dose_sensitivity_log<-ggplot(data=big_df_cases_bu)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name))+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", y=NULL, subtitle="Number of infections averted per dose (log scale)")+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=6)+
  scale_y_log10()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))

## FIGURE PLOT 

fig5<-cowplot::plot_grid(cases_averted_sensitivity,cases_averted_per_dose_sensitivity_log,line_plots_burundi, align = "v",
                         rel_heights=c(1.25,1,1.5),labels=c("A","B","C"), ncol=1)

fig5

ggsave(paste0("outputs/Figure_5.png"), fig5,
       width = 12.5, height = 15, scale=1.5, units="cm")


## BUJUMBURA APPENDIX PLOT - EXTENDED DOSE ANALYSIS 

big_df_cases_bu_app <- big_df_all |>  filter(Category %in% c("cases_inc_cum"), 
                                         TimeStep==max(big_df_all$TimeStep, na.rm = TRUE))

# percentage of available doses administered for each scenario
dose_percent_burundi_app<-big_df_cases_bu_app|>
  group_by(Category,scenario_name, scenario_dose,mean_cumulative_doses)|>
  summarise(dose_percent=mean_cumulative_doses/scenario_dose)

## CASES AVERTED SESNSITIVITY PLOTS
cases_averted_sensitivity_app<-ggplot(data=big_df_cases_bu_app)+
  geom_point(aes(x=mean_cumulative_doses/province_pop$bujumbura, y=mean_percent, colour=scenario_name))+#,position = position_dodge(width = 0.5)
  geom_line(aes(x=mean_cumulative_doses/province_pop$bujumbura, y=mean_percent, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  geom_errorbar(aes(x=mean_cumulative_doses/province_pop$bujumbura, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_name),width = 0)+ #,position = position_dodge(width = 0.5
  theme_minimal()+labs(x="Population vaccine coverage", y=NULL, title="Bujumbura" ,subtitle ="Percentage of infections averted")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent, breaks=waiver(), n.breaks=6)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none", plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"))
cases_averted_sensitivity_app

# CASES AVERTED PER DOSE SENSITIVITY
cases_averted_per_dose_sensitivity_app<-ggplot(data=big_df_cases_bu_app)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name))+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", y=NULL, subtitle="Number of infections averted per dose")+
  scale_x_continuous(labels = scales::comma, breaks =waiver(), n.breaks=6 )+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))
cases_averted_per_dose_sensitivity_app

cases_averted_per_dose_sensitivity_app_log<-ggplot(data=big_df_cases_bu_app)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name))+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", y=NULL, subtitle="Number of infections averted per dose (log scale)")+
  scale_x_continuous(labels = scales::comma, breaks =waiver(), n.breaks=6 )+
  scale_y_log10()+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))
cases_averted_per_dose_sensitivity_app_log

# LINE PLOT 
big_df_burundi_plots_app <- filter(big_df_all, scenario_num%in%c(62:67, NA))

big_df_burundi_app <- big_df_burundi_plots_app |>
  mutate(scenario_name = factor(scenario_name,
                                levels = c(scenario_names[1], sort(scenario_names[-1]))))

line_plots_burundi_app <- big_df_burundi_app |> 
  filter(Category =="cases_inc") |>
  mutate(Vaccination=scenario_name) |>
  ggplot() +
  geom_bar(data = fitting_data_bujumbura$data |>
             mutate(date = mpoxseir::mpoxseir_date_as_date(date)),
           aes(x = date, y = cases), colour = "grey50", fill=cols_region["burundi"],
           stat = "identity") +
  geom_line(aes(x = date, y = mean_value, color = Vaccination),linewidth=1,alpha=0.7) +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1,alpha=0.7) +
  labs(title =  "Weekly infections",
       subtitle = "Maximum doses",
       x = NULL,
       y = NULL) + 
  scale_x_date(date_breaks="4 months", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=13),
        plot.title = element_text(size=14),
        axis.title = element_text(size=13.5))+
  guides(colour=guide_legend(ncol=3))
line_plots_burundi_app

fig5_app<-cowplot::plot_grid(cases_averted_sensitivity_app,cases_averted_per_dose_sensitivity_app_log,line_plots_burundi_app, align = "v",
                         rel_heights=c(1.25,1,1.5),labels=c("A","B","C"), ncol=1)

fig5_app

ggsave(paste0("outputs/Figure_5_appendix_version.png"), fig5_app,
       width = 12.5, height = 15, scale=1.5, units="cm")

## DOSE ROLLOUT PLOTS ----------------------------------------------------------------------------------------------

# Define the plotting function
plot_cumulative_doses_by_dose <- function(data){
  # Define categories for grouping
  age_categories <- c("dose1_cumulative_00_04", "dose1_cumulative_05_14", "dose1_cumulative_15_plus",
                      "dose2_cumulative_00_04", "dose2_cumulative_05_14", "dose2_cumulative_15_plus")
  # Filter the data
  filtered_data <- data %>%
    filter(Category %in% age_categories)
  
  
  # Add readable labels
  labeled_data <- filtered_data %>%
    mutate( group = case_when(
      grepl("dose1", Category) & grepl("00_04", Category) ~ "Dose 1: 0-4",
      grepl("dose1", Category) & grepl("05_14", Category) ~ "Dose 1: 5-14",
      grepl("dose1", Category) & grepl("15_plus", Category) ~ "Dose 1: 15+",
      grepl("dose2", Category) & grepl("00_04", Category) ~ "Dose 2: 0-4",
      grepl("dose2", Category) & grepl("05_14", Category) ~ "Dose 2: 5-14",
      grepl("dose2", Category) & grepl("15_plus", Category) ~ "Dose 2: 15+"
    ))
  
  # Summarize total doses if "both" doses are selected
  summarized_data <- labeled_data %>%
    group_by(TimeStep, scenario_new, group) %>%
    summarise(total_cumulative_doses = sum(mean_value), .groups = "drop")
  
  summarized_data<-summarized_data%>%mutate(group=factor(group, levels=c("Dose 1: 0-4","Dose 1: 5-14","Dose 1: 15+",
                                                                         "Dose 2: 0-4","Dose 2: 5-14","Dose 2: 15+")))
  # Plot the data
  ggplot(summarized_data, aes(x = TimeStep, y = total_cumulative_doses, fill = group)) +
    geom_area(alpha=0.7) +
    facet_wrap(~ scenario_new) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste("Cumulative Doses Over Time"),
      x = "Time Step",
      y = "Cumulative Doses",
      fill = NULL,
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_y_continuous(labels = scales::comma)
  
}

plot_cumulative_doses_by_dose(data=big_df_sk)
plot_cumulative_doses_by_dose(data=big_df_eq)
plot_cumulative_doses_by_dose(data=big_df_burundi)

##### SUPPLEMENT PLOTS #####

## scenarios for deaths 

# Figure 3 deaths: Equateur ---------------------------------------------------------------------
region <-"equateur"

# Plot max dose scenario
big_df_plots_eq <- filter(big_df_equateur, scenario_num%in%c((max(big_df_equateur$scenario_num)-5):max(big_df_equateur$scenario_num)), ignore.case = TRUE)

big_df_eq <- big_df_plots_eq |>
  mutate(scenario_name = factor(scenario_name,
                                levels = c(scenario_names[1], sort(scenario_names[-1]))))


line_plot_equateur_deaths <- big_df_eq |> 
  filter(Category =="deaths_inc") |>
  filter(scenario_new %in% sort(unique(big_df_eq$scenario_new))[1:(6)]) |>
  mutate(Vaccination=factor(scenario_name, levels=sort(unique(big_df_eq$scenario_name))[1:(6)])) |>
  ggplot() +
  geom_bar(data = fitting_data_equateur$data |>
             filter(region==region) |>
             mutate(date = mpoxseir::mpoxseir_date_as_date(date)),
           aes(x = date, y = deaths), colour = "grey50", fill = cols_region["equateur"],
           stat = "identity") +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1,alpha=0.7) +
  geom_line(aes(x = date, y = mean_value, color = Vaccination),linewidth=1,alpha=0.7) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Weekly deaths",
       subtitle = "Maximum doses",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = NULL, 
                     values=c("No vaccination" ="black",
                              scenarios_palette_slides[c(1:7)]))+
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=13),
        plot.title = element_text(size=14),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13.5))
line_plot_equateur_deaths


###  BAR CHARTS OF CASES AVERTED

big_df_cases_eq_deaths <- big_df_equateur |>  filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_equateur$TimeStep, na.rm = TRUE))


## CASES AVERTED SESNSITIVITY PLOTS

cases_averted_sensitivity_eq_deaths<-ggplot(data=big_df_cases_eq_deaths)+
  geom_point(aes(x=mean_cumulative_doses/province_pop$equateur, y=mean_percent, colour=scenario_name), alpha=0.7)+#,position = position_dodge(width = 0.5)
  geom_line(aes(x=mean_cumulative_doses/province_pop$equateur, y=mean_percent, colour=scenario_name), alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  geom_errorbar(aes(x=mean_cumulative_doses/province_pop$equateur, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_name),width = 0)+ #,position = position_dodge(width = 0.5
  theme_minimal()+labs(x="Population vaccine covaerage", y=NULL,subtitle="Percentage of deaths averted", title="Equateur")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent, breaks=waiver(), n.breaks=6)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))



cases_averted_per_dose_sensitivity_eq_deaths<-ggplot(data=big_df_cases_eq_deaths)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name), alpha=0.7)+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", subtitle="Number of deaths averted per dose", y=NULL)+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=6)+
  scale_y_continuous(labels = scales::comma)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))


# FIGURE PLOT

fig3_deaths<-cowplot::plot_grid(cases_averted_sensitivity_eq_deaths,cases_averted_per_dose_sensitivity_eq_deaths,line_plot_equateur_deaths, align = "v",
                                rel_heights=c(1.25,1,1.5),labels=c("A","B","C"), ncol=1)

fig3_deaths

ggsave(paste0("outputs/Figure_equateur_deaths.png"), fig3_deaths,
       width = 12.5, height = 15, scale=1.5, units="cm")


# Figure 4 deaths: Sud Kivu ---------------------------------------------------------------------
region <-"sudkivu"

# Plot max dose scenario
big_df_plots_sk <- filter(big_df_sudkivu, scenario_num%in%c((max(big_df_sudkivu$scenario_num)-5):max(big_df_sudkivu$scenario_num)), ignore.case = TRUE)

big_df_sk <- big_df_plots_sk |>
  mutate(scenario_name = factor(scenario_name,
                                levels = c(scenario_names[1], sort(scenario_names[-1]))))


line_plot_sudkivu_deaths <- big_df_sk |> 
  filter(Category =="deaths_inc") |>
  filter(scenario_new %in% sort(unique(big_df_sk$scenario_new))[1:(6)]) |>
  mutate(Vaccination=factor(scenario_name, levels=sort(unique(big_df_sk$scenario_name))[1:(6)])) |>
  ggplot() +
  geom_bar(data = fitting_data_sudkivu$data |>
             filter(region==region) |>
             mutate(date = mpoxseir::mpoxseir_date_as_date(date)),
           aes(x = date, y = deaths), colour = "grey50", fill = cols_region["sudkivu"],
           stat = "identity") +
  geom_line(aes(x = date, y = mean_counter, color = "No vaccination"),linewidth=1,alpha=0.7) +
  geom_line(aes(x = date, y = mean_value, color = Vaccination),linewidth=1,alpha=0.7) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b %y", limits = as.Date(c(mpoxseir::mpoxseir_date_as_date(371), mpoxseir::mpoxseir_date_as_date(1099))))+
  labs(title = "Weekly deaths",
       subtitle = "Maximum doses",
       x = NULL,
       y = NULL ) + 
  scale_color_manual(name = NULL, 
                     values=c("No vaccination" ="black",
                              scenarios_palette_slides[c(1:7)]))+
  theme_minimal()+
  guides(colour=guide_legend(ncol=3))+
  theme(legend.position = "bottom",
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=13),
        plot.title = element_text(size=14),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13.5))
line_plot_sudkivu_deaths


###  BAR CHARTS OF CASES AVERTED

big_df_cases_sk_deaths <- big_df_sudkivu |>  filter(Category %in% c("deaths_inc_cum"), TimeStep==max(big_df_sudkivu$TimeStep, na.rm = TRUE))


## CASES AVERTED SESNSITIVITY PLOTS

cases_averted_sensitivity_sk_deaths<-ggplot(data=big_df_cases_sk_deaths)+
  geom_point(aes(x=mean_cumulative_doses/province_pop$sudkivu, y=mean_percent, colour=scenario_name), alpha=0.7)+#,position = position_dodge(width = 0.5)
  geom_line(aes(x=mean_cumulative_doses/province_pop$sudkivu, y=mean_percent, colour=scenario_name), alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  geom_errorbar(aes(x=mean_cumulative_doses/province_pop$sudkivu, ymin=lower_95_percent, ymax=upper_95_percent, colour=scenario_name),width = 0)+ #,position = position_dodge(width = 0.5
  theme_minimal()+labs(x="Population vaccine coverage", y=NULL,subtitle="Percentage of deaths averted", title="Sud Kivu")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent, breaks=waiver(), n.breaks=6)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))



cases_averted_per_dose_sensitivity_sk_deaths<-ggplot(data=big_df_cases_sk_deaths)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name), alpha=0.7)+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", subtitle="Number of deaths averted per dose", y=NULL)+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=4)+
  scale_y_continuous(labels = scales::comma)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))

cases_averted_per_dose_sensitivity_sk_deaths_log<-ggplot(data=big_df_cases_sk_deaths)+
  geom_point(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name), alpha=0.7)+
  geom_line(aes(x=mean_cumulative_doses, y=mean_averted_per_dose, colour=scenario_name),alpha=0.7)+#,position = position_dodge(width = 0.5)
  scale_color_manual(name = NULL, values= scenarios_palette_slides)+
  theme_minimal()+labs(x="Total doses administered", subtitle="Number of deaths averted per dose (log scale)", y=NULL)+
  scale_x_continuous(labels = scales::comma, breaks=waiver(), n.breaks=4)+
  scale_y_log10(labels = scales::comma)+
  theme(legend.text = element_text(size=13),
        axis.text = element_text(size=13),
        plot.subtitle = element_text(size=14),
        plot.title = element_text(size=15),
        axis.title = element_text(size=13.5))+
  theme(legend.position = "none",plot.margin = unit(c(0.1, 1, 0.1, 1), "cm"))

# FIGURE PLOT

fig4_deaths<-cowplot::plot_grid(cases_averted_sensitivity_sk_deaths,cases_averted_per_dose_sensitivity_sk_deaths_log,line_plot_sudkivu_deaths, align = "v",
                         rel_heights=c(1.25,1,1.5),labels=c("A","B","C"), ncol=1)

fig4_deaths

ggsave(paste0("outputs/Figure_sudkivu_deaths.png"), fig4_deaths,
       width = 12.5, height = 15, scale=1.5, units="cm")


##### Tables of averted outcomes ####### -------------------------------------------------


# Equateur data 
data_averted_cases_eq<-big_df_cases_eq[,c("mean_averted",
                                                  "lower_95_averted","upper_95_averted",
                                                  "scenario_new",
                                                  "mean_cumulative_doses",
                                                 "lower_95_cumulative_doses",
                                                 "upper_95_cumulative_doses",
                                                  "scenario_dose")]

data_averted_deaths_eq<-big_df_cases_eq_deaths[,c("mean_averted",
                                                 "lower_95_averted","upper_95_averted",
                                                 "scenario_new",
                                                 "mean_cumulative_doses",
                                                 "lower_95_cumulative_doses",
                                                 "upper_95_cumulative_doses",
                                                 "scenario_dose")]


colnames(data_averted_cases_eq)<-c("Infections averted mean",
                                    "Infections averted lower 95% CrI","Infections averted upper 95% CrI",
                                    "Scenario", 
                                   "Doses administered mean","Doses administered lower 95% CrI","Doses administered upper 95% CrI",
                                   "Doses available" )

colnames(data_averted_deaths_eq)<-c("Deaths averted mean",
                                   "Deaths averted lower 95% CrI","Deaths averted upper 95% CrI",
                                   "Scenario",
                                   "Doses administered mean","Doses administered lower 95% CrI","Doses administered upper 95% CrI",
                                   "Doses available" )

data_averted_eq<-full_join(data_averted_deaths_eq, data_averted_cases_eq, 
                                 by=c("Scenario", "Doses available", "Doses administered mean", "Doses administered lower 95% CrI","Doses administered upper 95% CrI"))|>
  mutate(Province="Equateur")

# Sud Kivu data 
data_averted_cases_sk<-big_df_cases_sk[,c("mean_averted",
                                                 "lower_95_averted","upper_95_averted",
                                                 "scenario_new",
                                                "mean_cumulative_doses",
                                                "lower_95_cumulative_doses",
                                                "upper_95_cumulative_doses",
                                                 "scenario_dose")]

data_averted_deaths_sk<-big_df_cases_sk_deaths[,c("mean_averted",
                                                         "lower_95_averted","upper_95_averted",
                                                         "scenario_new",
                                                         "mean_cumulative_doses",
                                                         "lower_95_cumulative_doses",
                                                         "upper_95_cumulative_doses",
                                                         "scenario_dose")]


colnames(data_averted_cases_sk)<-c("Infections averted mean",
                                   "Infections averted lower 95% CrI","Infections averted upper 95% CrI",
                                   "Scenario", "Doses administered mean","Doses administered lower 95% CrI","Doses administered upper 95% CrI",
                                   "Doses available" )

colnames(data_averted_deaths_sk)<-c("Deaths averted mean",
                                    "Deaths averted lower 95% CrI","Deaths averted upper 95% CrI",
                                    "Scenario", "Doses administered mean","Doses administered lower 95% CrI","Doses administered upper 95% CrI",
                                    "Doses available" )

data_averted_sk<-full_join(data_averted_deaths_sk, data_averted_cases_sk, 
                           by=c("Scenario", "Doses administered mean","Doses administered lower 95% CrI","Doses administered upper 95% CrI",
                                "Doses available"))|>
  mutate(Province="Sud Kivu", `Doses administered mean`=floor(`Doses administered mean`))


### DRC tables 
data_averted_drc <- full_join(data_averted_eq,data_averted_sk)
data_averted_drc <- data_averted_drc[order(data_averted_drc$`Doses available`),c(12,4,8,5,6,7,9,10,11,1,2,3)]
data_averted_drc <- data_averted_drc[order(data_averted_drc$Scenario),]
data_averted_drc <- data_averted_drc[order(data_averted_drc$Province),]

# DRC scenario A table
table_averted_A<- gt::gt(filter(data_averted_drc, Scenario == unique(data_averted_drc$Scenario)[1]))|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(26)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(25)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(30),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(27),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**DRC outcomes averted**"),
    subtitle = gt::md(paste0("Scenario ",unique(data_averted_drc$Scenario)[1])))|>
  gt::fmt_number(n_sigfig = 3)

table_averted_A_html <- gt::gtsave(table_averted_A, "outputs/table_averted_A.html")
webshot::webshot(table_averted_A_html, "outputs/table_averted_A.png", vwidth=1800)

# DRC scenario B table

table_averted_B<- gt::gt(filter(data_averted_drc, Scenario == unique(data_averted_drc$Scenario)[2]))|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(26)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(25)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(30),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(27),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**DRC outcomes averted**"),
    subtitle = gt::md(paste0("Scenario ",unique(data_averted_drc$Scenario)[2])))|>
  gt::fmt_number(n_sigfig = 3)

table_averted_B_html <- gt::gtsave(table_averted_B, "outputs/table_averted_B.html")
webshot::webshot(table_averted_B_html, "outputs/table_averted_B.png",vwidth=1800)

# DRC scenario C table

table_averted_C<- gt::gt(filter(data_averted_drc, Scenario== unique(data_averted_drc$Scenario)[3]))|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(26)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(25)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(30),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(27),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**DRC outcomes averted**"),
    subtitle = gt::md(paste0("Scenario ",unique(data_averted_drc$Scenario)[3])))|>
  gt::fmt_number(n_sigfig = 3)


table_averted_C_html <- gt::gtsave(table_averted_C, "outputs/table_averted_C.html")
webshot::webshot(table_averted_C_html, "outputs/table_averted_C.png",vwidth=1800)


# DRC scenario D table

table_averted_D<- gt::gt(filter(data_averted_drc, Scenario== unique(data_averted_drc$Scenario)[4]))|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(26)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(25)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(30),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(27),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**DRC outcomes averted**"),
    subtitle = gt::md(paste0("Scenario ",unique(data_averted_drc$Scenario)[4])))|>
  gt::fmt_number(n_sigfig = 3)


table_averted_D_html <- gt::gtsave(table_averted_D, "outputs/table_averted_D.html")
webshot::webshot(table_averted_D_html, "outputs/table_averted_D.png",vwidth=1800)

# DRC scenario E table

table_averted_E<- gt::gt(filter(data_averted_drc, Scenario== unique(data_averted_drc$Scenario)[5]))|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(26)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(25)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(30),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(27),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**DRC outcomes averted**"),
    subtitle = gt::md(paste0("Scenario ",unique(data_averted_drc$Scenario)[5])))|>
  gt::fmt_number(n_sigfig = 3)


table_averted_E_html <- gt::gtsave(table_averted_E, "outputs/table_averted_E.html")
webshot::webshot(table_averted_E_html, "outputs/table_averted_E.png",vwidth=1800)

# DRC scenario F table

table_averted_F<- gt::gt(filter(data_averted_drc, Scenario== unique(data_averted_drc$Scenario)[6]))|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(26)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(25)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(30),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(27),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**DRC outcomes averted**"),
    subtitle = gt::md(paste0("Scenario ",unique(data_averted_drc$Scenario)[6])))|>
  gt::fmt_number(n_sigfig = 3)


table_averted_F_html <- gt::gtsave(table_averted_F, "outputs/table_averted_F.html")
webshot::webshot(table_averted_F_html, "outputs/table_averted_F.png",vwidth=1800)


# Bujumbura data 
data_averted_cases_bu<-big_df_cases_bu_app[,c("mean_averted",
                                                 "lower_95_averted","upper_95_averted",
                                                 "scenario_new",
                                                 "mean_cumulative_doses",                                    
                                                 "lower_95_cumulative_doses",
                                                 "upper_95_cumulative_doses",
                                                 "scenario_dose")]



colnames(data_averted_cases_bu)<-c("Infections averted mean",
                                   "Infections averted lower 95% CrI","Infections averted upper 95% CrI",
                                   "Scenario", "Doses administered mean",  
                                   "Doses administered lower 95% CrI",  "Doses administered upper 95% CrI", 
                                   "Doses available" )
### Burundi tables 
data_averted_cases_bu<-data_averted_cases_bu|>
  mutate(Province="Bujumbura", `Doses administered mean`=floor(`Doses administered mean`))

data_averted_cases_bu <- data_averted_cases_bu[order(data_averted_cases_bu$`Doses available`),c(9,4,8,5,6,7,1,2,3)]
data_averted_cases_bu <- data_averted_cases_bu[order(data_averted_cases_bu$Scenario),]
data_averted_cases_bu <- data_averted_cases_bu[order(data_averted_cases_bu$Province),]

table_averted_bu<- gt::gt(data_averted_cases_bu)|>
  gt::tab_style(
    style = gt::cell_text(
      font = gt::system_fonts(name = "industrial"),
      size = gt::px(52)
    ),
    locations = list(
      gt::cells_body(columns = Province),
      gt::cells_column_labels()
    )
  )|>
  gt::tab_style(
    style = gt::cell_text(
      size = gt::px(50)
    ),
    locations = list(
      gt::cells_body()
    )
  )|>
  gt::tab_options(  heading.title.font.size = gt::px(60),
                    heading.title.font.weight = "bold",
                    heading.subtitle.font.size = gt::px(54),
                    heading.subtitle.font.weight = "bold")|>
  gt::tab_header(
    title = gt::md("**Burundi outcomes averted**"),
    subtitle = gt::md("All Scenarios"))|>
  gt::fmt_number(n_sigfig = 3)

html_averted_bu <- gt::gtsave(table_averted_bu, "outputs/table_averted_bu.html")
webshot::webshot(html_averted_bu, "outputs/table_averted_bu.png", vwidth=4500)

