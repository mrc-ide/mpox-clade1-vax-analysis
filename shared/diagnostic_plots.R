get_par_labels <- function(par_names) {
  
  regions <- character(length(par_names))
  regions[grepl("<equateur>", par_names)] <- '" (Equateur)"'
  regions[grepl("<sudkivu>", par_names)] <- '" (Sud-Kivu)"'
  
  labels <- c(alpha_cases = 'alpha^C',
              alpha_deaths = 'alpha^D',
              beta_z_max = 'zeta',
              cfr_00_04 = 'theta["[0,5)"]',
              lambda = 'mu',
              phi_05_14 = 'phi["[5,15)"]',
              phi_15_plus = 'phi["[15+)"]',
              prop_SW = 'p[SW]',
              rho = 'rho',
              R0_hh = 'R[0]^H',
              R0_sw_st = 'R[0]^"SW,PBS"'
  )
  
  par_names <- gsub("<.*", "", par_names)
  
  par_labels <- paste0(labels[par_names], regions)
  
}

pairs_plot <- function(samples) {
  
  par_names <- rownames(samples$pars)
  samples$pars <- array(samples$pars, c(dim(samples$pars), 1))
  rownames(samples$pars) <- par_names
  
  samples_df <- posterior::as_draws_df(samples)
  samples_df$.iteration <- NULL
  samples_df$.chain <- NULL
  samples_df$.draw <- NULL
  
  par_labels <- get_par_labels(par_names)
  
  GGally::ggpairs(samples_df,
                  progress = FALSE, 
                  columnLabels = par_labels,
                  labeller = "label_parsed") + 
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill="white"),
      axis.text.x = element_text(size = 11, angle = 270, 
                                 vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 11))
}

traceplots <- function(samples) {
  samples$pars <- samples$pars_full
  
  par_names <- rownames(samples$pars)
  par_labels <- get_par_labels(par_names)
  rownames(samples$pars) <- par_labels
  
  if ("density_full" %in% names(samples)) {
    samples$pars <- abind::abind(samples$pars, samples$density_full, along = 1)
    rownames(samples$pars)[nrow(samples$pars)] <- '"log posterior density"'
  }
  
  samples_df <- posterior::as_draws_df(samples)
  
  color_scheme <- unname(unlist(rev(bayesplot::color_scheme_get("viridis"))))
  color_scheme <- gsub("*FF", "", color_scheme)
  
  bayesplot::color_scheme_set(color_scheme)
  p <- bayesplot::mcmc_trace(samples_df,
                             n_warmup = samples$control$pmcmc$n_burnin,
                             facet_args = list(nrow = 3,
                                               labeller = label_parsed)) +
    bayesplot::theme_default(base_size = 16) +
    theme(legend.position = "none")
  p
}

rankplots <- function(samples) {
  samples$pars <- samples$pars_full
  rownames(samples$pars) <- 
    gsub(">", "", gsub("<", "_", rownames(samples$pars)))
  
  samples$pars <- samples$pars[, -seq_len(samples$control$pmcmc$n_burnin), ]
  
  par_names <- rownames(samples$pars)
  par_labels <- get_par_labels(par_names)
  rownames(samples$pars) <- par_labels
  
  samples_df <- posterior::as_draws_df(samples)
  
  color_scheme <- unname(unlist(rev(bayesplot::color_scheme_get("viridis"))))
  color_scheme <- gsub("*FF", "", color_scheme)
  
  bayesplot::color_scheme_set(color_scheme)
  p <- bayesplot::mcmc_rank_overlay(samples_df,
                                    facet_args = list(nrow = 3,
                                                      labeller = label_parsed)) +
    bayesplot::theme_default(base_size = 16) +
    theme(legend.position = "none")
  p
}
