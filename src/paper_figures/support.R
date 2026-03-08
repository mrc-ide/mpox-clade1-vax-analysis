mean_ci <- function(x, alpha = 0.05) {
  p <- c(alpha / 2, 1 - alpha / 2)
  q <- quantile(x, c(0.025, 0.975))
  names(q) <- paste0("q", p * 100)
  
  c(mean = mean(x), q)
}

tidy_pars <- function(pmcmc_results, region = NULL) {
  pars_tidy <- pmcmc_results$pars |> 
    t() |> 
    as.data.frame() |> 
    pivot_longer(cols = everything())
  
  is_combined <- length(dim(pmcmc_results$observations$trajectories)) == 4
  if (is_combined) {
    pattern <- "(.*)<(.*)>.*"
    pars_tidy <- pars_tidy |> 
      mutate(region = sub(pattern, "\\2", name),
             name = sub(pattern, "\\1", name),
             region = factor(if_else(region == name, "both", region)))
  } else {
    pars_tidy <- pars_tidy |> 
      mutate(region = region)
  }
  pars_tidy
}

plot_posteriors <- function(pmcmc_results, region = "both") {
  
  pars_tidy <- tidy_pars(pmcmc_results, region)
  
  
  if (region == "both") {
    
    g1 <- pars_tidy |> 
      filter(region != "both") |> 
      ggplot(aes(x = value, fill = region)) +
      geom_histogram(alpha = 0.7, colour = "transparent", position = "identity",
                     show.legend = TRUE) +
      scale_fill_manual(values = cols_region, drop = FALSE) +
      facet_wrap(vars(name), scales = "free") +
      theme_bw()
    
    g2 <- pars_tidy |> 
      filter(region == "both") |> 
      ggplot(aes(x = value, fill = region)) +
      geom_histogram(alpha = 0.7, colour = "transparent", position = "identity",
                     show.legend = TRUE) +
      scale_fill_manual(values = cols_region, drop = FALSE) +
      facet_wrap(vars(name), scales = "free") +
      theme_bw()
    
    g <- g1 + g2 +
      plot_layout(widths = c(3, 1), guides = "collect") +
      plot_annotation(tag_levels = "A") &
      theme(legend.position = "bottom",
            legend.direction = "horizontal")
    
  } else {
    
    g <- pars_tidy |> 
      ggplot(aes(x = value, fill = region)) +
      geom_histogram(alpha = 0.7, colour = "transparent", position = "identity",
                     show.legend = FALSE) +
      scale_fill_manual(values = cols_region, drop = FALSE) +
      facet_wrap(vars(name), scales = "free", nrow = 2) +
      theme_bw()
  }
  
  g
}

extract_par <- function(pars_unpacked, region, both, name, simplify = TRUE) {
  if(both=="joint"){ 
    sapply(pars_unpacked, function(x) x[[region]][[name]], simplify = simplify)
  }else{
    sapply(pars_unpacked, function(x) x[[name]], simplify = simplify)
  }
}

extract_daily_zoonotic_region <- function(pars_unpacked, region, both, summarise = TRUE) {
  N0 <- extract_par(pars_unpacked, region, both, "N0")[, 1]
  beta_z <- extract_par(pars_unpacked, region,both, "beta_z")
  ret <- beta_z * N0
  format_age_spline(ret, summarise, names(N0))
}

format_age_spline <- function(x, summarise, nms) {
  
  if (summarise) {
    ret <- t(apply(x, 1, mean_ci))
  } else {
    ret <- x
  }
  
  as.data.frame(ret) |> 
    mutate(group = factor(nms, nms, ordered = TRUE))
}


extract_pars_all_regions <- function(pars_unpacked, f) {
  nms_region <- c("Equateur" = "equateur", "Sud-Kivu" = "sudkivu")
  ret <- lapply(nms_region, f, pars_unpacked = pars_unpacked)
  names(ret) <- nms_region
  ret |> 
    bind_rows(.id = "region") 
}


extract_cfr_by_age <- function(pars_unpacked, region,both, summarise = TRUE, name = "CFR") {
  cfr_matrix <- extract_par(pars_unpacked, region,both, name = name, simplify = FALSE)
  idx <- mpoxseir::get_compartment_indices()
  ret <- sapply(cfr_matrix, "[", , idx$vax$unvaccinated)
  nms <- rownames(ret)
  rownames(ret) <- NULL
  format_age_spline(ret, summarise, nms)
}

tidy_trajectories <- function(pmcmc_results, fitting_data, region) {
  
  
  if (region == "both") {
    trajectories <- apply(pmcmc_results$observations$trajectories, c(1, 2, 3), mean_ci)
    grid <- dimnames(trajectories)
  
  } else {
    
    trajectories <- apply(pmcmc_results$observations$trajectories, c(1, 2), mean_ci)
    grid <- dimnames(trajectories)
    grid <- append(grid, values = region, after = 2)
    
  }

  names(grid) <- c("name", "state", "region", "date")
  grid$date <- mpoxseir_date_as_date(sort(unique(fitting_data$date))) 
  ret <- expand.grid(grid)
  ret$value <- c(trajectories)

  ret
}



tidy_Rt <- function(pmcmc_results, fitting_data, region) {
  
  if (region == "both") {
    rt <- apply(pmcmc_results$Rt, c(1,2), mean_ci)
    grid <- dimnames(rt)
    
  } else {
    
    rt <- apply(pmcmc_results$Rt, c(1), mean_ci)
    grid <- dimnames(rt)
    grid <- append(grid, values = region, after = 1)
    
  }
  
  names(grid) <- c( "state", "region", "date")
  grid$date <- mpoxseir_date_as_date(sort(unique(fitting_data$date))) 
  ret <- expand.grid(grid)
  ret$value <- c(rt)
  
  ret
}

extract_state <- function(pars_unpacked, region, pmcmc_results) {
  
  sys <- dust2::dust_system_create(
    mpoxseir::model_targeted_vax, pars_unpacked[[1]], time = 0,
    n_particles = 1, n_groups = ifelse(region == "both", 2, 1),
    seed = 1, dt = 1, deterministic = TRUE,
    preserve_particle_dimension = TRUE)
  full_index <- dust2::dust_unpack_index(sys)
  
  nms <- c("General" = "hh", "Zoonotic" = "z", "Sexual" = "s", "Nosocomial" = "hc")
  nms_state <- paste0("cases_cumulative_", nms)
  idx_state <- unlist(full_index[nms_state])
  
  if (region == "both") {
    state <- pmcmc_results$observations$state[idx_state, , ]
    nms_region <- c("equateur", "sudkivu")
  } else {
    if (region == "equateur") state <- pmcmc_results$observations$state[idx_state, ]
    if (region == "sudkivu") state <- pmcmc_results$observations$state[idx_state, ]
    if (region == "burundi") state <- pmcmc_results$observations$state[idx_state, ]
    
    nms_region <- region
  }
  
  n_particles <- tail(dim(state), 1)
  ret <- expand.grid(state = names(nms),
                     region = nms_region,
                     particle = seq_len(n_particles))
  
  ret$value <- c(state)
  ret 
}

tidy_ascertainment <- function(pmcmc_results, region) {
  
  if (region == "both") {
    
    traj <- apply(pmcmc_results$observations$trajectories, c(1, 2, 4), sum)
    traj2 <- apply(traj, c(2), function(x){x["observed_cases_inc",]/ x["cases_inc",]})
    ret <- apply(traj2, c(2), mean_ci)
    
  } else {
    
    traj <- apply(pmcmc_results$observations$trajectories, c(1, 3), sum)
    traj2 <- traj["observed_cases_inc",]/ traj["cases_inc",]
    ret <- mean_ci(traj2) |> data.frame() 
    colnames(ret) <- region
  }
  
  ret
}
