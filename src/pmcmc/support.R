
get_data <- function(x, x_linelist, x_who_shiny, start_date, end_date,
                     fit_by_age, fit_KPs, r) {
  
  # Use data scraped from WHO shiny when age-disaggregation not available
  raw_data <- bind_rows(x, filter(x_who_shiny, week_start > max(x$week_start)))

  if (fit_KPs) {
    raw_data <- dplyr::left_join(raw_data, 
        dplyr::select(x_linelist,
                      region, epi_week, cases_total, cases_HCW, cases_SW))
  }

  if (r == "both") {
    data <- raw_data %>%
      dplyr::filter(region %in% c("equateur", "sudkivu"))
  } else {
    data <- raw_data %>%
      dplyr::filter(region == r)
    
    data$region <- NULL
  }
  
  end_date <- mpoxseir::mpoxseir_date(end_date)
  SW_end_date <- mpoxseir::mpoxseir_date("2024-04-01")
  data <- as.data.frame(data) %>%
    dplyr::mutate(date = mpoxseir::mpoxseir_date(week_start + 6)) %>%
    dplyr::mutate(cases_SW = ifelse(date < SW_end_date, cases_SW, NA)) %>%
    dplyr::mutate(cases_HCW = NA) %>%
    dplyr::filter(date <= end_date) %>%
    dplyr::select(-c(week_start, epi_week)) %>%
    dplyr::relocate(date)
  
  
  ## Add rows of NA between start_date and beginning of data
  mpoxseir_start_date <- mpoxseir::mpoxseir_date(start_date) 
  dates_na <- 
    seq(ceiling((mpoxseir_start_date + 1) / 7) * 7, min(data$date) - 1, by = 7)
  if (r == "both") {
    data <- data %>% 
      dplyr::rows_append(data.frame(date = dates_na, region = "equateur")) %>%
      dplyr::rows_append(data.frame(date = dates_na, region = "sudkivu")) %>%
      dplyr::arrange(region, date)
  } else {
    data <- data %>% 
      dplyr::rows_append(data.frame(date = dates_na)) %>%
      dplyr::arrange(date)
  }

  if (any(data$date %% 7 != 0)) {
    stop("Weekly data not correctly aligned")
  }
  
  mask <- data
  mask[, !(names(mask) %in% c("date", "region"))] <- TRUE

  if (fit_by_age) {
    idx_mask <- !is.na(data$cases_00_04)
    mask[idx_mask, c("cases", "deaths")] <- FALSE
  } else {
    mask[, c("cases_00_04", "cases_05_14", "cases_15_plus",
             "deaths_00_04", "deaths_05_14", "deaths_15_plus")] <- FALSE
  }
  
  if (!fit_KPs) {
    mask$cases_total <- mask$cases_HCW <- mask$cases_SW <- FALSE
  }
  

  
  ## We are saving:
  ## 1. data: a data frame with a date column (in terms of mpoxseir dates),
  ##          a region column if running a multiregion fit, and then columns
  ##          for each datastream
  ## 2. mask: a data frame with the same columns as above, but the columns
  ##          for each datastream identify whether each datapoint is fitted
  ##          to (TRUE) or not (FALSE)
  list(data = data,
       mask = mask)
}


get_prior <- function(region, names){
  
  ## Priors for prop_SW have a = 2 and then b is set to ensure mean of 0.03 for
  ## Sud-Kivu and 0.007 for Equateur
  if (region == "both") {
    domain <- array(0, c(length(names), 2))
    rownames(domain) <- names
    domain[, 2] <- ifelse(grepl("^(alpha|cfr|phi|prop)", names), 1, 1000)
    monty_model(
      list(
        parameters = names,
        density = function(x) {
          names(x) <- names
          sum(dbeta(x[grepl("^alpha", names(x))], 1, 9, log = TRUE)) +
            sum(dbeta(x[grepl("^(cfr|phi)", names(x))], 1, 1, log = TRUE)) +
            sum(dunif(x[!grepl("^(alpha|cfr|phi|prop)", names(x))], 0, 1000, log = TRUE)) +
            dbeta(x["prop_SW<equateur>"], 2, 2 * (1 / 0.007 - 1), log = TRUE) +
            dbeta(x["prop_SW<sudkivu>"], 2, 2 * (1 / 0.03 - 1), log = TRUE)},
        domain = domain
      ))
  } else if (region == "equateur") {
    monty::monty_dsl({
      alpha_cases ~ Beta(1, 9)
      alpha_deaths ~ Beta(1, 9)
      cfr_00_04 ~ Beta(1, 1)
      beta_z_max ~ Uniform(0, 1000)
      lambda ~ Uniform(0, 1000)
      phi_05_14 ~ Beta(1, 1)
      phi_15_plus ~ Beta(1, 1)
      prop_SW ~ Beta(2, 2 * (1 / 0.007 - 1))
      R0_hh ~ Uniform(0, 1000)
      R0_sw_st ~ Uniform(0, 1000)
    })
  } else {
    monty::monty_dsl({
      alpha_cases ~ Beta(1, 9)
      alpha_deaths ~ Beta(1, 9)
      cfr_00_04 ~ Beta(1, 1)
      lambda ~ Uniform(0, 1000)
      phi_05_14 ~ Beta(1, 1)
      phi_15_plus ~ Beta(1, 1)
      prop_SW ~ Beta(2, 2 * (1 / 0.03 - 1))
      R0_hh ~ Uniform(0, 1000)
      R0_sw_st ~ Uniform(0, 1000)
    })
  }
}


create_mcmc_pars <- function(region, deterministic, assumptions) {
  type <- if (deterministic) "deterministic" else "stochastic"
  pars_folder <- paste("parameters", region, type, assumptions, sep = "/")
  info_file <- paste(pars_folder, "info.csv", sep = "/")
  proposal_file <- paste(pars_folder, "proposal.csv", sep = "/")
  info <- read.csv(info_file)
  proposal <- read.csv(proposal_file)
  proposal_names <- proposal[, 1]
  proposal <- as.matrix(proposal[, -1])
  row.names(proposal) <- colnames(proposal) <- proposal_names

  
  if (region == "both") {
    if (!identical(gsub("<both>", "", paste0(info$name, "<", info$region, ">")), 
                   row.names(proposal))) {
      stop("info and proposal files do not match")
    }
    
    initial <- info$initial
    names(initial) <- row.names(proposal)
  
  } else {
    
    if (!identical(info$name, row.names(proposal))) {
      stop("info and proposal files do not match")
    }
    initial <- info$initial
    names(initial) <- info$name
    
  }
  
  if (!deterministic) {
    proposal <- 2.38^2 / nrow(proposal) * proposal
  }

  list(initial = initial,
       vcv = proposal)
}


create_filter <- function(data, deterministic, region, start_date, control) {
  ## Here we construct the data for fitting by taking the full data
  ## and apply the mask to it 
  fit_data <- data$data
  data_names <- names(fit_data)[!(names(fit_data) %in% c("date", "region"))]
  for (nm in data_names) {
    fit_data[[nm]] <- ifelse(data$mask[[nm]], fit_data[[nm]], NA)
  }
  
  if (region == "both") {
    fit_data <- dust2::dust_filter_data(fit_data, time = "date",
                                        group = "region")
  } else {
    fit_data <- dust2::dust_filter_data(fit_data, time = "date")
  }
  
  sys <- mpoxseir::model_targeted_vax()
  time_start <- mpoxseir::mpoxseir_date(start_date)
  
  if (deterministic) {
    dust2::dust_unfilter_create(sys, time_start, fit_data, dt = 1)

  } else {
    dust2::dust_filter_create(sys, time_start, fit_data, dt = 1,
                              n_particles = control$filter$n_particles,
                              n_threads = control$filter$n_threads)
  }
}

fit_control <- function(deterministic, n_steps, n_burnin, n_sample, n_chains,
                        n_particles) {
  adaptive_proposal <- deterministic

  thinning_factor <- floor((n_steps - n_burnin) / (n_sample / n_chains))

  n_threads_total <- control_cores()
  if (deterministic) {
    n_workers <- min(n_chains, n_threads_total)
  } else {
    max_workers <- 4
    pos <- seq_len(max_workers)
    n_workers <- max(pos[n_threads_total %% pos == 0 & pos <= n_chains])
  }
  n_threads <- n_threads_total / n_workers
  
  rerun_every <- if (deterministic) Inf else 100

  filter <- list(n_particles = n_particles,
                 n_threads = n_threads)

  pmcmc <- list(n_steps = n_steps,
                n_chains = n_chains,
                n_burnin = n_burnin,
                thinning_factor = thinning_factor,
                n_threads_total = n_threads_total,
                n_workers = n_workers,
                save_trajectories = TRUE,
                rerun_every = rerun_every,
                rerun_random = TRUE, 
                adaptive_proposal = adaptive_proposal,
                progress = TRUE)

  list(
    filter = filter,
    pmcmc = pmcmc
    )
}


control_cores <- function() {
  as.integer(Sys.getenv("CONTEXT_CORES",
                        Sys.getenv("MC_CORES",
                                   getOption("mc.cores", 1))))
}


save_pars_inputs <- function(samples, mcmc_pars, region) {
  if (region == "both") {
    i <- which.max(samples$density)
    
    initial <- samples$pars[, i]
    region <- ifelse(grepl(".*<(.+)>", names(initial)),
                     gsub(".*<(.+)>", "\\1", names(initial)),
                     "both")
    info <- data.frame(name = gsub("<.*", "", names(initial)),
                       region = region,
                       initial = initial)
  } else {
    i <- which.max(samples$density)
    
    initial <- samples$pars[, i]
    info <- data.frame(name = names(initial),
                       initial = initial)
  }
  
  proposal <- cov(t(samples$pars))
  
  write.csv(info, "outputs/info.csv", row.names = FALSE)
  write.csv(proposal, "outputs/proposal.csv")
}


run_pmcmc <- function(filter, packer, mcmc_pars, control, deterministic, region,
                      snapshots) {
  
  index <- save_index(rt = TRUE)
  
  snapshots <- mpoxseir::mpoxseir_date(snapshots)
  
  likelihood <- dust2::dust_likelihood_monty(filter, packer,
                                             save_state = TRUE,
                                             save_trajectories = index$idx,
                                             save_snapshots = snapshots)
  
  prior <- get_prior(region, packer$names())
  
  density <- likelihood + prior 
  
  vcv <- mcmc_pars$vcv
  initial <- mcmc_pars$initial
  
  n_steps <- control$n_steps
  n_chains <- control$n_chains
  
  rerun_every <- control$rerun_every
  rerun_random <- control$rerun_random
  
  if (deterministic) {
    sampler <- monty::monty_sampler_adaptive(vcv, initial_vcv_weight = 100,
                                             min_scaling = 0.9)
  } else {
    sampler <- monty::monty_sampler_random_walk(vcv, rerun_every = rerun_every,
                                                rerun_random = rerun_random)
  }
  
  runner <- monty::monty_runner_callr(control$n_workers)
  samples <- monty::monty_sample(density, sampler, n_steps,
                                 initial = initial, n_chains = n_chains,
                                 runner = runner)
  
  rownames(samples$observations$trajectories) <- index$names
  if (region == "both") {
    colnames(samples$observations$trajectories) <- filter$groups
    colnames(samples$observations$state) <- filter$groups
    colnames(samples$observations$snapshots) <- filter$groups
  }
  ## save the packer for downstream use
  samples$packer <- packer
  samples
}


thin_samples <- function(samples, control) {
  burnin <- control$pmcmc$n_burnin
  thinning_factor <- control$pmcmc$thinning_factor
  
  traj_names <- dimnames(samples$observations$trajectories)
  state_names <- dimnames(samples$observations$state)
  snapshots_names <- dimnames(samples$observations$snapshots)
  
  pars_full <- samples$pars
  density_full <- samples$density
  
  samples <- monty::monty_samples_thin(samples, thinning_factor, burnin)
  dimnames(samples$observations$trajectories) <- traj_names
  dimnames(samples$observations$state) <- state_names
  dimnames(samples$observations$snapshots) <- snapshots_names
  
  
  samples$pars <- array_flatten(samples$pars, c(2, 3))
  
  samples$density <- c(samples$density)
  
  if (region == "both") {
    samples$observations$state <- 
      array_flatten(samples$observations$state, c(3, 4))
    
    samples$observations$trajectories <- 
      array_flatten(samples$observations$trajectories, c(4, 5))
    
    samples$observations$snapshots <- 
      array_flatten(samples$observations$snapshots, c(4, 5))
  } else {
    samples$observations$state <- 
      array_flatten(samples$observations$state, c(2, 3))
    
    samples$observations$trajectories <- 
      array_flatten(samples$observations$trajectories, c(3, 4))
    
    samples$observations$snapshots <- 
      array_flatten(samples$observations$snapshots, c(3, 4))
  }
  
  samples$pars_full <- pars_full
  samples$density_full <- density_full
  samples$control <- control
  
  samples
}


array_flatten <- function(x, i) {
  dx <- dim(x)
  if (any(i < 1 | i > length(dx))) {
    stop(sprintf("Values of 'i' must be in [1, %d]", length(dx)))
  }
  if (length(i) < 2) {
    stop("i must be vector of at least length 2")
  }
  if (any(diff(i) != 1)) {
    stop("All values of 'i' must be consecutive integers")
  }
  dx[[i[[1L]]]] <- prod(dx[i])
  dx_new <- dx[seq_along(dx)[-i[-1]]]
  if (length(dx_new) == 1L) {
    dx_new <- NULL
  }
  nms_x <- dimnames(x)
  if (!is.null(nms_x)) {
    nms_x[i[1]] <- list(NULL)
    nms_x[i[-1]] <- NULL
  }
  dim(x) <- dx_new
  dimnames(x) <- nms_x
  x
}
