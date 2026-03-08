
orderly2::orderly_dependency(
  name = "pmcmc",
  query = "latest(parameter:region=='sudkivu')",
  files = c("inputs/samples_sudkivu.rds" = "outputs/samples.rds",
            "inputs/fitting_data_sudkivu.rds" = "outputs/fitting_data.rds"))
orderly2::orderly_dependency(
  name = "pmcmc",
  query = "latest(parameter:region=='equateur')",
  files = c("inputs/samples_equateur.rds" = "outputs/samples.rds",
            "inputs/fitting_data_equateur.rds" = "outputs/fitting_data.rds"))
orderly2::orderly_dependency(
  name = "pmcmc_burundi", 
  query =  "latest()", 
  files = c("inputs/samples_burundi.rds" = "outputs/samples.rds",
            "inputs/fitting_data_burundi.rds" = "outputs/fitting_data.rds"))

orderly2::orderly_resource("support.R")
orderly2::orderly_shared_resource("calc-Rt.R")
source("support.R")
source("calc-Rt.R")

library(tidyverse)
library(patchwork)
library(mpoxseir)
library(ggsflabel)
library(ggplot2)

cols_region <- c(MetBrewer::met.brewer("Cassatt1", 7)[c(5,1)],MetBrewer::met.brewer("Demuth", 3)[3], MetBrewer::met.brewer("Cassatt1", 7)[3]) 
names(cols_region) <- c("equateur", "sudkivu", "both", "burundi")
nms_region <- c(equateur = "DRC: Equateur", sudkivu = "DRC: Sud Kivu", burundi = "Burundi: Bujumbura")
cols_transmission <- MetBrewer::met.brewer("Hokusai3", 4)


# fitting plots ----------------------------------------------------------------
both <- "separate"

# data preparation
# read in data dependencies
pmcmc_results_e <- readRDS("inputs/samples_equateur.RDS")
fitting_data_e <- readRDS("inputs/fitting_data_equateur.RDS")
pmcmc_results_s <- readRDS("inputs/samples_sudkivu.RDS")
fitting_data_s <- readRDS("inputs/fitting_data_sudkivu.RDS")
pmcmc_results_b <- readRDS("inputs/samples_burundi.RDS")
fitting_data_b <- readRDS("inputs/fitting_data_burundi.RDS")

fitting_data_e <- fitting_data_e$data |> mutate(region="equateur") |> dplyr::select(date, cases, region) |> dplyr::rename(cases_total= cases) #cases_total
fitting_data_s <- fitting_data_s$data |> mutate(region="sudkivu") |> dplyr::select(date, cases, region) |> dplyr::rename(cases_total= cases) #cases_total
fitting_data_b <- fitting_data_b$data |> mutate(region="burundi") |> dplyr::select(date, cases, region) |> dplyr::rename(cases_total= cases) ## CHECK CASES IS CORRECT OUTPUT NOT CASES_BINOM
fitting_data <- rbind(fitting_data_e, fitting_data_s, fitting_data_b)

pars_tidy_e <- tidy_pars(pmcmc_results_e)
pars_tidy_s <- tidy_pars(pmcmc_results_s)
pars_tidy_b <- tidy_pars(pmcmc_results_b)

pars_tidy <- rbind(pars_tidy_e |> mutate(region="equateur"),
                   pars_tidy_s |> mutate(region="sudkivu"),
                   pars_tidy_b |> mutate(region="burundi"))

pars_unpacked_e <- apply(pmcmc_results_e$pars, 2, pmcmc_results_e$packer$unpack)
pars_unpacked_s <- apply(pmcmc_results_s$pars, 2, pmcmc_results_s$packer$unpack)
pars_unpacked_b <- apply(pmcmc_results_b$pars, 2, pmcmc_results_b$packer$unpack)

daily_zoonotic_summary_e <- extract_daily_zoonotic_region(pars_unpacked_e, region="equateur",both=both) 
# write.csv(daily_zoonotic_summary_e,
#          "outputs/zoonotic_eq.csv",row.names = FALSE)
daily_zoonotic_summary_s <- extract_daily_zoonotic_region(pars_unpacked_s, region="sudkivu",both=both) 
daily_zoonotic_summary_b <- extract_daily_zoonotic_region(pars_unpacked_b, region="burundi", both=both) 

daily_zoonotic_summary_tidy <- rbind(daily_zoonotic_summary_e |> mutate(region="equateur"),
                                     daily_zoonotic_summary_s |> mutate(region="sudkivu"),
                                     daily_zoonotic_summary_b |> mutate(region="burundi"))

cfr_prior_summary <- extract_cfr_by_age(pars_unpacked_e, region="equateur", both=both, name = "CFR_0") |> 
  mutate(region = "both")
cfr_summary_e <- extract_cfr_by_age(pars_unpacked_e, region="equateur", both=both, name = "CFR") |> 
  mutate(region = "equateur")
cfr_summary_s <- extract_cfr_by_age(pars_unpacked_s, region="sudkivu", both=both, name = "CFR") |> 
  mutate(region = "sudkivu")

cfr_summary_tidy <- rbind(cfr_summary_e, cfr_summary_s) |> bind_rows(cfr_prior_summary)

write.csv(cfr_summary_tidy %>% 
            mutate(region = ifelse(region=="both","prior",region)),
          "outputs/cfr.csv",row.names=FALSE)

state_e <- extract_state(pars_unpacked_e, region = "equateur", pmcmc_results=pmcmc_results_e) 
state_s <- extract_state(pars_unpacked_s, region = "sudkivu", pmcmc_results=pmcmc_results_s) 
state_b <- extract_state(pars_unpacked_b, region = "burundi", pmcmc_results=pmcmc_results_b) 

state_tidy <- rbind(state_e, state_s, state_b) |> 
  group_by(state, region) |> 
  summarise(value = mean(value)) |> 
  group_by(region) |> 
  mutate(prop = value / sum(value)) |>
  rowwise() |>
  mutate(region = nms_region[region]) |>
  tibble()

prop_transmission_with_uncertainty <- rbind(state_e,state_s,state_b) %>%
  group_by(region,particle) %>%
  mutate(prop = value/sum(value)) %>% group_by(region,state) %>%
  summarise(mean = round(mean(prop,na.rm=TRUE),2),
            lower = round(quantile(prop,0.025,na.rm=TRUE),2),
            upper = round(quantile(prop,0.975,na.rm=TRUE),2))
write.csv(prop_transmission_with_uncertainty,
          "outputs/prop_transmission_with_uncertainty.csv",row.names=FALSE)

trajectories_e <- tidy_trajectories(pmcmc_results_e, fitting_data_e, region = "equateur")
trajectories_s <- tidy_trajectories(pmcmc_results_s, fitting_data_s, region = "sudkivu")
trajectories_b <- tidy_trajectories(pmcmc_results_b, fitting_data_b, region = "burundi")
trajectories <- rbind(trajectories_e, trajectories_s, trajectories_b)

ascert_e <- tidy_ascertainment(pmcmc_results_e, region = "equateur")
ascert_s <- tidy_ascertainment(pmcmc_results_s, region = "sudkivu")
ascert_b <- tidy_ascertainment(pmcmc_results_b, region = "burundi")
ascert <- cbind(ascert_e, ascert_s, ascert_b)
write.csv(ascert, paste0("outputs/ascertainment_", both,".csv"))

r0_e <- calc_R0(pmcmc_results_e, "equateur", no_historic_vax = TRUE)
r0_s <- calc_R0(pmcmc_results_s, "sudkivu", no_historic_vax = TRUE)
r0_b <- calc_R0(pmcmc_results_b, "bujumbura", no_historic_vax = TRUE)
r0 <-rbind(r0_e|>as.data.frame()|>mutate(region="equateur"), 
           r0_s|>as.data.frame()|>mutate(region="sudkivu"),
           r0_b|>as.data.frame()|>mutate(region="burundi"))|>
  group_by(region)|>
  pivot_longer(names_to = "name", values_to = "value", cols=1:3)


rt_e <- tidy_Rt(pmcmc_results_e, fitting_data_e, region = "equateur")
rt_s <- tidy_Rt(pmcmc_results_s, fitting_data_s, region = "sudkivu")
rt_b <- tidy_Rt(pmcmc_results_b, fitting_data_b, region = "burundi")
rt <- rbind(rt_e, rt_s, rt_b)

rt |>
  filter(date == min(rt$date) | date == max(rt$date)) |>
  pivot_wider(names_from="state") -> rt_time

write.csv(rt_time, paste0("outputs/rt_start_end_", both,".csv"))
write.csv(state_tidy, paste0("outputs/transmission_route_values_", both,".csv"))


# plotting
g_zoonotic <- daily_zoonotic_summary_tidy |> 
  filter(str_detect(group, "[0-9]")) |> 
  ggplot(aes(x = group, group = region, colour = region, fill = region)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.7, col = "transparent") +
  geom_line(aes(y = mean)) +
  scale_fill_manual(values = cols_region[c(-3)], drop = FALSE,
                    labels = nms_region,
                    aesthetics = c("fill", "colour")) +
  theme_bw() +
  labs(x = "Age group (years)", y = "Expected daily zoonotic cases",
       fill = "Region", colour = "Region") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

g_zoonotic2 <- daily_zoonotic_summary_tidy |> 
  filter(str_detect(group, "[0-9]")) |> 
  filter(region=="equateur")|>
  ggplot(aes(x = group, group = region, colour = region)) +
  geom_segment(aes(y = q2.5, yend = q97.5)) +
  geom_point(aes(y = mean)) +
  scale_colour_manual(values = cols_region[c(-3)], 
                      labels = nms_region) +
  theme_bw() +
  labs(x = "Age group (years)", y = "Expected daily\nzoonotic cases",
       colour = "Region") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        text=element_text(size=12),
        #legend.position.inside = "inside",
        legend.position = "none")#c(0.8,0.8))


# new R0 plot
g_R0 <- r0 |> 
  mutate(name=case_when(name=="R0_gen_pop"~ "R0 general population",
                        name=="R0_sex"~ "R0 sexual networks",
                        name=="R0"~"R0")) |>
  ggplot(aes(x = value, fill = region)) +
  geom_density(alpha = 0.7, colour = "transparent", position = "identity",
               show.legend = TRUE) +
  scale_fill_manual(values = cols_region[c(-3)], 
                    labels = c(nms_region, both="Both"), 
                    drop = FALSE) +
  labs(y="Density", fill=NULL, x="R0")+
  facet_grid(rows = vars(name), scales = "free") +
  geom_vline(xintercept = 1, lty = 2) +
  theme_bw()+
  theme(strip.background = element_blank(),
        text=element_text(size=12))

g_SW <- pars_tidy |> 
  filter(name %in% c("prop_SW")) |> 
  mutate(name = "Sex worker population") |> 
  mutate(value=value*100)|>
  ggplot(aes(x = value, fill = region)) +
  geom_density(alpha = 0.7, colour = "transparent", position = "identity",
               show.legend = TRUE) +
  scale_fill_manual(values = cols_region[c(-3)], 
                    labels = c(nms_region, both="Both"), 
                    drop = FALSE) +
  labs(y="Density", fill="Region", x="Percentage %")+
  facet_grid(rows = vars(name), scales = "free") +
  theme_bw()+
  theme(strip.background = element_blank(),
        text=element_text(size=12),
        legend.position = "inside",
        legend.position.inside=c(0.8,0.8))

g_phi <- pars_tidy |> 
  filter(name %in% c("phi_05_14", "phi_15_plus")) |> 
  mutate(name=case_when(name=="phi_05_14"~ "Case ascertainment 5-14",
                        name=="phi_15_plus"~ "Case ascertainment 15+")) |>
  mutate(name=factor(name, levels=c("Case ascertainment 5-14", "Case ascertainment 15+")))|>
  mutate(value=value*100)|>
  ggplot(aes(x = value, fill = region)) +
  geom_density(alpha = 0.9, colour = "transparent", position = "identity",
               show.legend = TRUE) +
  scale_fill_manual(values = cols_region[c(-3)], 
                    labels = c(nms_region), 
                    drop = FALSE) +
  labs(y="Density", fill="Region", x="Percentage %")+
  facet_grid(rows = vars(name), scales = "free") +
  theme_bw()+
  theme(strip.background = element_blank(),
        text=element_text(size=12),
        legend.position = "none")

g_cfr1 <- cfr_summary_tidy |> 
  filter(str_detect(group, "[0-9]")) |> 
  ggplot(aes(x = group, group = region, colour = region)) +
  geom_segment(aes(y = q2.5, yend = q97.5)) +
  geom_point(aes(y = mean)) +
  scale_colour_manual(values = cols_region[-4],
                      labels =  c(nms_region, both="Historic clade I")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Age group (years)", y = "Case Fatality Ratio (CFR)",
       colour = "Region")

g_cfr2 <- cfr_summary_tidy |> 
  filter(str_detect(group, "[0-9]")) |> 
  ggplot(aes(x = group, group = region, colour = region, fill = region)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.7, col = "transparent") +
  geom_line(aes(y = mean)) +
  scale_fill_manual(values = cols_region[c(3:1)], 
                    labels =  c(nms_region[c(1:2)],both="Historic clade I"),
                    drop = FALSE,
                    aesthetics = c("fill", "colour")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Age group (years)", y = "Case Fatality Ratio (CFR)",
       colour = "", fill = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        text=element_text(size=12),
        #legend.position.inside = "inside",
        legend.position = "none")

g_transmission <- state_tidy |> 
  filter(state!="Nosocomial") |>
  ggplot(aes(x = prop, y = region, fill = state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols_transmission) +
  scale_x_continuous(label = scales::percent) +
  theme_bw() +
  labs(x = "\nProportion of cumulative infections",
       fill = "Transmission\nroute", y = "")

ggsave(paste0("outputs/transmission_route_", both, ".png"), g_transmission, 
       width = 10, height = 5,
       units = "cm", scale = 1.2)

ggsave(paste0("outputs/severity_", both,".png"), g_cfr2, 
       width = 10, height = 5,
       units = "cm", scale = 1.5)

g <- ((g_R0 | (g_zoonotic2)) /
        g_transmission) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", heights = c(3, 1)) 

ggsave(paste0("outputs/phi_", both, ".png"), g_phi,
       width = 17, height = 10, scale = 1.5,
       units = "cm")

# plots of trajectories for Figure 1
g_obs_under_cases_s <- trajectories |>
  pivot_wider() |>
  filter(state == "observed_cases_inc" | state == "cases_inc") |>
  filter(region=="sudkivu") |>
  filter(date>"2023-12-31") |>
  ggplot(aes(x = date)) +
  geom_bar(data = fitting_data |>
             filter(region=="sudkivu") |>
             mutate(date = mpoxseir_date_as_date(date)),
           aes(y = cases_total), fill=cols_region["sudkivu"],#fill = "grey50",
           stat = "identity") +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = state), alpha = 0.5,
              col = "transparent") +
  geom_line(aes(y = mean, col = state)) +
  scale_fill_manual(values = cols_transmission[c(2,3)],
                    labels = c("observed_cases_inc" = "Modelled cases",
                               "cases_inc" = "Modelled Infections"),
                    aesthetics = c("fill", "colour")) +
  scale_x_date(date_labels = "%b-%y") +
  theme_bw() +
  labs(x = "", y = "Weekly total",title="DRC: Sud Kivu", col=NULL, fill=NULL)+
  theme(legend.position="bottom", text=element_text(size=15),
        plot.title = element_text(hjust=0.5))

g_obs_under_cases_e <- trajectories |>
  pivot_wider() |>
  filter(state == "observed_cases_inc" | state == "cases_inc") |>
  filter(region=="equateur") |>
  filter(date>"2023-12-31") |>
  ggplot(aes(x = date)) +
  geom_bar(data = fitting_data |>
             filter(region=="equateur") |>
             mutate(date = mpoxseir_date_as_date(date)),
           aes(y = cases_total), fill=cols_region["equateur"],#fill = "grey50",
           stat = "identity") +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = state), alpha = 0.5,
              col = "transparent") +
  geom_line(aes(y = mean, col = state)) +
  scale_fill_manual(values = cols_transmission[c(2,3)],
                    labels = c("observed_cases_inc" = "Modelled cases",
                               "cases_inc" = "Modelled Infections"),
                    aesthetics = c("fill", "colour")) +
  scale_x_date(date_labels = "%b-%y") +
  theme_bw() +
  labs(x = "", y = "Weekly total",title="DRC: Equateur", col=NULL, fill=NULL)+
  theme(legend.position="bottom", text=element_text(size=15),
        plot.title = element_text(hjust=0.5))

g_obs_under_cases_b <- trajectories |>
  pivot_wider() |>
  filter(state == "observed_cases_inc" | state == "cases_inc") |>
  filter(region=="burundi") |>
  filter(date>"2023-12-31") |>
  ggplot(aes(x = date)) +
  geom_bar(data = fitting_data |>
             filter(region=="burundi") |>
             mutate(date = mpoxseir_date_as_date(date)),
           aes(y = cases_total), fill = cols_region["burundi"], # "grey50",
           stat = "identity") +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = state), alpha = 0.5,
              col = "transparent") +
  geom_line(aes(y = mean, col = state)) +
  scale_fill_manual(values = cols_transmission[c(2,3)],
                    labels = c("observed_cases_inc" = "Modelled cases",
                               "cases_inc" = "Modelled Infections"),
                    aesthetics = c("fill", "colour")) +
  scale_x_date(date_labels = "%b-%y") +
  theme_bw() +
  labs(x = "", y = "Weekly total", title="Burundi: Bujumbura",col=NULL, fill=NULL)+
  theme(legend.position="bottom", text=element_text(size=15),
        plot.title = element_text(hjust=0.5))+
  coord_cartesian(xlim=c(as.Date("2024-01-01"), max(trajectories_b$date)))


ggsave(paste0("outputs/cases_fits_underlying_sudkivu_", both,".png"), g_obs_under_cases_s,
       width = 10, height = 8, scale = 1.5,
       units = "cm")
ggsave(paste0("outputs/cases_fits_underlying_equateur_", both,".png"), g_obs_under_cases_e,
       width = 10, height = 8, scale = 1.5,
       units = "cm")
ggsave(paste0("outputs/cases_fits_underlying_burundi_", both,".png"), g_obs_under_cases_b,
       width = 10, height = 8, scale = 1.5,
       units = "cm")

g_Rt <- rt |>
  pivot_wider(names_from="state") |>
  ggplot(aes(x = date)) +
  facet_grid(vars(region), scales = "free_y",
             labeller = labeller(.rows = nms_region)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = region), alpha = 0.5,
              col = "transparent") +
  geom_line(aes(y = mean, col = region)) +
  geom_hline(aes(yintercept=1),linetype="dashed")+
  scale_fill_manual(values = cols_region,
                    labels = nms_region,
                    aesthetics = c("fill", "colour"), guide = "none") +
  scale_x_date(date_labels = "%b-%y") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  labs(x = "", y = "Effective reproduction number")

ggsave(paste0("outputs/rt_", both,".png"), g_Rt,
       width = 10, height = 10, scale = 1.5,
       units = "cm")


# map plot ---------------------------------------------------------------------
orderly2::orderly_shared_resource("DRC_shapefile")
orderly2::orderly_shared_resource("Tanzania_shapefile")
orderly2::orderly_shared_resource("Burundi_shapefile")
orderly2::orderly_shared_resource("Uganda_shapefile")
orderly2::orderly_shared_resource("Rwanda_shapefile")
sf_d_1 <- sf::st_read("DRC_shapefile/gadm41_COD_1.shp")
sf_d_0 <- sf::st_read("DRC_shapefile/gadm41_COD_0.shp")
sf_b_0 <- sf::st_read("Burundi_shapefile/gadm41_BDI_0.shp") 
sf_b_1 <- sf::st_read("Burundi_shapefile/gadm41_BDI_1.shp") 
sf_t <- sf::st_read("Tanzania_shapefile/gadm41_TZA_0.shp") |> mutate(loc="Tanzania")
sf_u <- sf::st_read("Uganda_shapefile/gadm41_UGA_0.shp") |> mutate(loc="Uganda")
sf_r <- sf::st_read("Rwanda_shapefile/gadm41_RWA_0.shp") |> mutate(loc="Rwanda") 

sf_d_1 |> mutate(loc=case_when(NAME_1=="Équateur" ~ "DRC: Equateur", 
                               NAME_1=="Sud-Kivu" ~ "DRC: Sud Kivu",
                               .default=NA)) -> sf_d_1
sf_b_1 |> mutate(loc=case_when(NAME_1=="Bujumbura Mairie" | NAME_1=="Bujumbura Rural" 
                               ~ "Burundi: Bujumbura",
                               .default=NA)) -> sf_b_1

map <- ggplot() +
  geom_sf(data=sf_d_1, aes(fill=loc, col=loc)) +
  geom_sf(data=sf_d_0, fill=NA, col="black") +
  geom_sf(data=sf_b_1, aes(fill=loc, col=loc)) +
  geom_sf(data=sf_b_0, fill=NA, col="black") +
  geom_sf(data=sf_t, fill="lightgrey", col="black") +
  geom_sf(data=sf_u, fill="lightgrey", col="black") +
  geom_sf(data=sf_r, fill="lightgrey", col="black") +
  
  scale_fill_manual(breaks=c("DRC: Equateur","DRC: Sud Kivu",  "Burundi: Bujumbura"),
                    values = c(cols_region[[1]], cols_region[[2]],cols_region[[4]]), 
                    na.value = "lightgrey")+
  
  scale_colour_manual(breaks=c("DRC: Equateur","DRC: Sud Kivu",  "Burundi: Bujumbura"),
                      values = c(cols_region[[1]], cols_region[[2]], cols_region[[4]]), 
                      na.value = "lightgrey")+
  
  geom_sf_label_repel(data=sf_b_1[sf_b_1$NAME_1!="Bujumbura Rural",], aes(label = loc), 
                      size = 5, fill=cols_region[[4]],
                      nudge_x = 8, nudge_y = 10) +
  geom_sf_label_repel(data=sf_d_1[sf_d_1$loc=="DRC: Equateur",], aes(label = loc), 
                      size = 5, fill=cols_region[[1]],
                      nudge_x = -4, nudge_y = 7) +
  geom_sf_label_repel(data=sf_d_1[sf_d_1$loc=="DRC: Sud Kivu",], aes(label = loc), 
                      size = 5, fill=cols_region[[2]],
                      nudge_x = -9, nudge_y = -9) +
  theme_void()+
  labs(fill=NULL, col=NULL)+
  theme(legend.position = "none")

ggsave("outputs/map.png", map,
       width = 10, height = 10, scale = 1.5,
       units = "cm")



# Figure 1 ---------------------------------------------------------------------
# Map
# Phi
# fits/cases/infections
legend_top <- cowplot::get_legend(g_obs_under_cases_b)
fig1_a <- cowplot::plot_grid(map, 
                             g_phi + theme(legend.position = "none"),
                             ncol=1, labels="AUTO")
fig1_b <- cowplot::plot_grid(legend_top,
                             g_obs_under_cases_e + theme(legend.position = "none"), 
                             g_obs_under_cases_s + theme(legend.position = "none"),
                             g_obs_under_cases_b + theme(legend.position = "none"),
                             ncol=1, labels=c("C","","D", "E"),
                             rel_heights = c(0.2,1,1,1))

fig1 <- cowplot::plot_grid(fig1_a, 
                           fig1_b, 
                           ncol=2, rel_widths = c(0.8,1)) 
fig1

ggsave("outputs/Figure_1.png", fig1,
       width = 15, height = 15, scale = 1.5,
       units = "cm")

# Figure 2 ---------------------------------------------------------------------
# Rt hist
# zoonotic
# bars
# CFR

g_transmission + theme(legend.position = "bottom",legend.text = element_text(size=12),
                       text=element_text(size=12))+
  labs(fill=NULL, 
       x="Proportion of cumulative infections\nby transmission route")+
  scale_y_discrete(labels=c('DRC\nSud Kivu', 'DRC\nEquateur', 'Burundi\nBujumbura')) -> g_transmission2


fig2_left <- cowplot::plot_grid(g_R0 + 
                                  theme(legend.text = element_text(size=10),
                                        legend.position=c(0.75,0.5)),
                                g_SW + labs(fill=NULL) + theme(legend.position = "none"),
                                
                                ncol=1, labels=c("A","B"),
                                rel_heights = c(0.9,0.325), align="v")

fig2_rightb <- cowplot::plot_grid(
  g_zoonotic2,
  g_cfr2 + theme(legend.text = element_text(size=10),
                 legend.position = c(0.7,0.7)),
  ncol=1,
  labels=c("D","E"), align="v",
  rel_heights = c(0.9,0.9))

fig2_right <- cowplot::plot_grid(g_transmission2,
                                 fig2_rightb,
                                 ncol=1,
                                 labels=c("C",""),
                                 rel_heights = c(0.9,1.8))


fig2 <- cowplot::plot_grid(fig2_left,
                           fig2_right,
                           rel_heights = c(1,0.9),
                           ncol=2,
                           labels=NULL, align="v") 
fig2

ggsave("outputs/Figure_2.png", fig2,
       width = 20, height = 16, scale = 1.5,
       units = "cm")

# R0 summaries
r0_summary <- r0 |> group_by(region, name) |> summarise(mean=mean(value), lower_95=quantile(value, 0.025), upper_95=quantile(value,0.975))

