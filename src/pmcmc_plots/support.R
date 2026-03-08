mean_ci <- function(x, alpha = 0.05) {
  p <- c(alpha / 2, 1 - alpha / 2)
  q <- quantile(x, c(0.025, 0.975))
  names(q) <- paste0("q", p * 100)
  
  c(mean = mean(x), q)
}

pmcmc_tidy_trajectories <- function(object, data, summarise = FALSE, summary = NULL) {

  state <- object

  if (is.null(summary)) {
    summary <- function(x) {
      quantiles <- c(0.025, 0.5, 0.975)
      ret <- c(mean(x), quantile(x, quantiles, names = FALSE))
      set_names(ret, c("mean", paste0("q", quantiles * 100)))
    }
  }

  # create a named list with an entry for each dimension of `state`
  # the list entries will be dimnames, where they exist, and numeric otherwise
  grid <- list(state = rownames(state) %||% seq_len(nrow(state)))

  if (summarise) {
    # marginalise over particles
    margin <- if (object$nested) c(1, 2, 4) else c(1, 3)
    state_summary <- apply(state, margin, summary)
    # statistic moves from 1st to penultimate dimension
    perm <- if (object$nested) c(2, 3, 1, 4) else c(2, 1, 3)
    # replace state with summarised array with dimensions, where dimension
    # previously containing particles now contains summary statistics, i.e.
    # [state, time, (if nested) population, statistic]
    state <- aperm(state_summary, perm)

    # define summary_function and a vector of output names
    grid$statistic <- rownames(state_summary) %||% seq_len(nrow(state_summary))
  } else {
    grid$particle <- seq_len(nrow(object$pars))
  }
  grid$time <- mpoxseir::mpoxseir_date_as_date(data$data$date)

  ret <- do.call(expand.grid, grid)
  ret$value <- c(state)
  ret
}



pmcmc_tidy_chains <- function(object) {
  pars <- pmcmc_tidy_chain_one(object, type = "pars_full")
  probabilities <- pmcmc_tidy_chain_one(object, type = "probabilities_full")
  rbind(pars, probabilities)
}

pmcmc_tidy_chain_one <- function(object, type) {
  stopifnot(type %in% c("pars_full", "probabilities_full"))

  n_pars_full <- nrow(object$pars_full)
  n_chains <- max(object$chain %||% 1)
  n_iterations <- n_pars_full / n_chains
  chain <- rep(seq_len(n_chains), each = n_iterations)
  iteration <- rep(seq_len(n_iterations), n_chains)

  x <- object[[type]]
  grid <- list(particle = seq_len(nrow(object$pars_full)),
               type = type,
               name = colnames(x))
  if (object$nested) {
    grid$population <- dimnames(x)[[3]]
  }
  ret <- do.call(expand.grid, grid)
  ret$iteration <- iteration
  ret$chain <- chain
  ret$value <- c(x) # drops unwanted attributes
  ret
}


plot_traceplots <- function(tidy_chains, region = NULL) {
  if (!is.null(region)) {
    tidy_chains <- tidy_chains %>%
      filter(population == region)
  }
  
  tidy_chains %>%
    ggplot(aes(x = iteration, y = value, group = chain, colour = chain)) +
    geom_line() +
    theme_bw() +
    facet_wrap(vars(name), scales = "free")
}


plot_trajectories <- function(samples, data, type, region = NULL) {
  get_summary <- function(x) {
    quantiles <- c(0.025, 0.975)
    ret <- c(mean(x), quantile(x, quantiles, names = FALSE))
    set_names(ret, c("mean", "lb", "ub"))
  }
  
  traj <- samples$observations$trajectories
  group_names <- c("all", "00_04", "05_14", "15_plus")
  group_titles <- c("All", "0 to 4", "5 to 14", "15+")
  
  if (!is.null(region)) {
    traj <- traj[, region, , ]
    data$data <- data$data[data$data$region == region, ]
    data$mask <- data$mask[data$mask$region == region, ]
  }
  
  dates <- mpoxseir::mpoxseir_date_as_date(data$data$date)
  
  summarise_1 <- function(group) {
    suffix <- if (group == "all") "" else paste0("_", group)
    traj_name <- if (type == "cases") paste0("observed_cases_inc", suffix) else
      paste0("deaths_inc", suffix)
    group_summary <- t(apply(traj[traj_name, , ], 1, get_summary))
    data_name <- paste0(type, suffix)
    fitted_data <- ifelse(data$mask[[data_name]], data$data[[data_name]], NA)
    unfitted_data <- ifelse(!data$mask[[data_name]], data$data[[data_name]], NA)
    
    data.frame(date = dates, 
               group_summary,
               fitted_data = fitted_data,
               unfitted_data = unfitted_data)
  }
  
  res <- lapply(group_names, summarise_1)
  names(res) <- group_titles
  
  res <- dplyr::bind_rows(res, .id = "group")
  
  ylab <- paste0("Weekly ", type)
  
  ggplot(res, aes(x = date)) +
    geom_line(aes(y = mean), col = "blue") +
    geom_ribbon(alpha = 0.2, aes(ymin = lb, ymax = ub), fill = "blue") +
    geom_point(aes(y = fitted_data), col = "green4", na.rm = TRUE)  +
    geom_point(aes(y = unfitted_data), col = "red3", na.rm = TRUE)  +
    scale_x_date(date_breaks = "2 months", date_minor_breaks = "1 months", 
                 date_labels = "%b %y") +
    facet_wrap(group ~ ., scales = "free") +
    labs(y = ylab, x = "") +
    theme_bw(base_size = 16) +
    theme(plot.margin = margin(5.5, 17.5, 5.5, 5.5),
          strip.background = element_rect(fill = "white"))
  
}


plot_trajectories_KPs <- function(samples, data, region = NULL) {

  get_summary <- function(x) {
    quantiles <- c(0.025, 0.975)
    ret <- c(mean(x, na.rm = TRUE), 
             quantile(x, quantiles, na.rm = TRUE, names = FALSE))
    set_names(ret, c("mean", "lb", "ub"))
  }
  
  get_data_summary <- function(group, fitted = TRUE) {
    data_name <- paste0("cases_", group)
    pos <- ifelse(data$mask[[data_name]] == fitted, data$data[[data_name]], NA)
    tot <- ifelse(data$mask$cases_total == fitted, data$data$cases_total, NA)
    na_data <- is.na(pos) | is.na(tot)
    pos[na_data] <- 0
    tot[na_data] <- 0
    
    ret <- Hmisc::binconf(pos, tot, alpha = 0.05, return.df = TRUE) * 100
    data_name <- if (fitted) "fitted_data" else "unfitted_data"
    names(ret) <- paste0(data_name, "_", c("mean", "lb", "ub"))
    ret
  }
  
  traj <- samples$observations$trajectories
  group_names <- c("HCW", "SW")
  
  if (!is.null(region)) {
    traj <- traj[, region, , ]
    data$data <- data$data[data$data$region == region, ]
    data$mask <- data$mask[data$mask$region == region, ]
  }
  
  dates <- mpoxseir::mpoxseir_date_as_date(data$data$date)
  
  summarise_1 <- function(group) {
    suffix <- paste0("_", group)
    prop <- traj[paste0("observed_cases_inc_", group), , ] /
      traj["observed_cases_inc", , ] * 100
    group_summary <- t(apply(prop, 1, get_summary))
  
    fitted_data_summary <- get_data_summary(group, TRUE)
    unfitted_data_summary <- get_data_summary(group, FALSE)
    
    data.frame(date = dates, 
               group_summary,
               fitted_data_summary,
               unfitted_data_summary)
  }
  
  res <- summarise_1("SW")
  
  ylab <- paste0("SW cases as percentage of weekly cases")
  
  ggplot(res, aes(x = date)) +
    geom_line(aes(y = mean), col = "blue") +
    geom_ribbon(alpha = 0.2, aes(ymin = lb, ymax = ub), fill = "blue") +
    geom_pointrange(aes(y = fitted_data_mean, ymin = fitted_data_lb, 
                        ymax = fitted_data_ub), col = "green4", na.rm = TRUE) +
    geom_pointrange(aes(y = unfitted_data_mean, ymin = unfitted_data_lb, 
                        ymax = unfitted_data_ub), col = "red3", na.rm = TRUE) +
    scale_x_date(date_breaks = "2 months", date_minor_breaks = "1 months",
                 date_labels = "%b %y") +
    labs(y = ylab, x = "") +
    theme_bw(base_size = 16)
  
}


plot_fitted_by_group <- function(samples,data, type, region = NULL) {
  
  traj <- samples$observations$trajectories
  dates <- mpoxseir::mpoxseir_date_as_date(data$data$date)
  
  if (!is.null(region)) {
    traj <- traj[, region, , ]
    }
  
  # select non age groupings
  groups <- c("SW", "HCW", "PBS", "all")
  
  name_lookup <- data.frame(group = c("SW","HCW", "PBS", "gen_pop"),
                            name = c("sex workers", "healthcare workers","people who buy sex", "general population"))
  
  
  summarise_1 <- function(group) {
    suffix <- if (group == "all") "" else paste0("_", group)
    group_summary <- rowMeans(traj[paste0(type, "_inc", suffix), , ])
    
    data.frame(date = dates, 
               mean = group_summary)
  }
  
  res <- lapply(groups, summarise_1)
  names(res) <- groups
  
  res <- dplyr::bind_rows(res, .id = "group")
  res %>% filter(group!="all") %>% group_by(date) %>% reframe(minus=sum(mean)) %>%
    left_join(res%>%filter(group=="all") %>% select(!group))%>%
    mutate(mean = mean - minus, group="gen_pop") %>%
    select(!minus) %>%
    rbind(res) %>%
    filter(group!="all")%>%
    left_join(name_lookup)->res
  
  
  ylab <- paste0("Weekly ", type)
  
  ggplot(res, aes(x = date)) +
    geom_col(aes(y=mean,fill=name)) +
    scale_x_date(date_breaks = "1 months", date_labels = "%b %y") +
    labs(y = ylab, x = "")+
    theme(legend.position="bottom", 
          legend.title = element_blank(),
          text=element_text(size=15))
  
}



plot_cfr <- function(pmcmc_results, data, region = NULL) {
  
  trajectories <- pmcmc_results$observations$trajectories

  if (!is.null(region)) {
    trajectories <- trajectories[, region, , ]
    data$data <- data$data[data$data$region == region, ]
    data$mask <- data$mask[data$mask$region == region, ]
  }
  
  nms <- rownames(trajectories)
  nms_cases <- grep("cases_cum", nms, value = TRUE)
  nms_deaths <- grep("deaths_cum", nms, value = TRUE)
  idx_max_t <- ncol(trajectories)
  
  cases <- trajectories[nms_cases, idx_max_t, ]
  deaths <- trajectories[nms_deaths, idx_max_t, ]
  cfr <- deaths / cases
  
  res <- data.frame(name = gsub("cases_cumulative_", "", nms_cases),
                    t(apply(cfr, 1, mean_ci)))
  res <- res %>%
    dplyr::mutate(name = case_when(name == "00_04" ~ "0 to 4",
                                   name == "05_14" ~ "5 to 14",
                                   name == "15_plus" ~ "15+"))
  res$name <- factor(res$name, c("0 to 4", "5 to 14", "15+"))
  
  res |> 
    ggplot(aes(x = name, y = mean)) +
    geom_point(col = "blue") +
    geom_segment(alpha = 0.5, aes(y = q2.5, yend = q97.5), col = "blue") +
    geom_point() +
    geom_point(
      data = data$data |> 
        dplyr::select(contains("cfr")) |> 
        dplyr::filter(!is.na(cfr_00_04)) |> 
        pivot_longer(everything()) |> 
        mutate(name = case_when(name == "cfr_00_04" ~ "0 to 4",
                                name == "cfr_05_14" ~ "5 to 14",
                                name == "cfr_15_plus" ~ "15+")),
      aes(y = value),
      colour = "green4") +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Case Fatality Ratio (CFR)", x = "Age Group") +
    theme_bw(base_size = 16)
}

plot_Rt <- function(samples, data, region = NULL) {

  get_summary <- function(x) {
    quantiles <- c(0.025, 0.975)
    ret <- c(mean(x), quantile(x, quantiles, names = FALSE))
    set_names(ret, c("mean", "lb", "ub"))
  }
  
  rt <- samples$Rt
  
  if (!is.null(region)) {
    rt <- rt[region, , ]
    data$data <- data$data[data$data$region == region, ]
  }
  
  
  dates <- mpoxseir::mpoxseir_date_as_date(data$data$date)
  
  res <- data.frame(date = dates,
                    t(apply(rt, 1, get_summary)))
  
  ylab <- "Rt"
  
  ggplot(res, aes(x = date)) +
    geom_line(aes(y = mean), col = "blue") +
    geom_ribbon(alpha = 0.2, aes(ymin = lb, ymax = ub), fill = "blue") +
    scale_x_date(date_breaks = "1 months", date_labels = "%b %y") +
    labs(y = ylab, x = "")
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


plot_SW_cases_and_cfr <- function(samples, data, region = NULL) {
  
  SW <- plot_trajectories_KPs(samples, data, region)
  CFR <- plot_cfr(samples, data, region)
  
  (SW | CFR) +
    plot_annotation(tag_levels = list(c('A', 'B')))
  
}
