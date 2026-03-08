
get_data <- function(x, start_date, end_date) {
  end_date <- mpoxseir::mpoxseir_date(end_date)

  data <- as.data.frame(x) %>%
    dplyr::rename(cases_binom = cases,
                  cases = cases_total) %>% 
    dplyr::mutate(date = mpoxseir_date(week_ending),
                  cases_00_04_binom = cases_00_04,
                  cases_00_14_binom = cases_00_04 + cases_05_14) %>%
    dplyr::mutate(cases_00_04 = NA,
                  cases_05_14 = NA,
                  cases_15_plus = NA) %>%
    dplyr::filter(date <= end_date) %>%
    dplyr::relocate(date) %>%
    dplyr::select(-c(week_ending, cases_check))
  
  ## Add rows of NA between start_date and beginning of data
  mpoxseir_start_date <- mpoxseir::mpoxseir_date(start_date) 
  dates_na <- 
    seq(ceiling((mpoxseir_start_date + 1) / 7) * 7, min(data$date) - 1, by = 7)
  data <- data %>% 
    dplyr::rows_append(data.frame(date = dates_na)) %>%
    dplyr::arrange(date)

  if (any(data$date %% 7 != 0)) {
    stop("Weekly data not correctly aligned")
  }
  
  mask <- data
  # include everything to start with
  mask[, !(names(mask) %in% c("date", "region"))] <- TRUE
  
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


get_prior <- function(){
  monty::monty_dsl({
    alpha_cases ~ Beta(1, 9)
    phi_05_14 ~ Beta(1, 1)
    phi_15_plus ~ Beta(1, 1)
    prop_SW ~ Beta(2, 2 * (1 / 0.018 - 1))
    R0_hh ~ Gamma(shape = 18.01306, scale = 0.02295737)
    R0_sw_st ~ Uniform(0, 10)
    rho ~ Beta(1, 1)
  }, gradient = FALSE)
}


create_mcmc_pars <- function(deterministic, region) {
  type <- if (deterministic) "deterministic" else "stochastic"
  info_file <- paste("parameters", region, type, "info.csv", sep = "/")
  proposal_file <- paste("parameters", region, type, "proposal.csv", sep = "/")
  info <- read.csv(info_file)

  proposal <- read.csv(proposal_file)
  proposal_names <- proposal[, 1]
  proposal <- as.matrix(proposal[, -1])
  row.names(proposal) <- colnames(proposal) <- proposal_names
  
  rownames(info) <- info$name
  
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


create_filter <- function(data, deterministic, start_date, control) {
  ## Here we construct the data for fitting by taking the full data
  ## and apply the mask to it 
  fit_data <- data$data
  data_names <- names(fit_data)[!(names(fit_data) %in% c("date", "region"))]
  for (nm in data_names) {
    fit_data[[nm]] <- ifelse(data$mask[[nm]], fit_data[[nm]], NA)
  }
  
  fit_data <- dust2::dust_filter_data(fit_data, time = "date")
  
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


run_pmcmc <- function(filter, packer, mcmc_pars, control, deterministic,
                      snapshots) {
  index <- save_index(rt = TRUE)
  
  snapshots <- mpoxseir::mpoxseir_date(snapshots)
  
  likelihood <- dust2::dust_likelihood_monty(filter, packer,
                                             save_state = TRUE,
                                             save_trajectories = index$idx,
                                             save_snapshots = snapshots)
  
  prior <- get_prior()
  
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
  
  samples$observations$state <- 
    array_flatten(samples$observations$state, c(2, 3))
  
  samples$observations$trajectories <- 
    array_flatten(samples$observations$trajectories, c(3, 4))
  
  samples$observations$snapshots <- 
    array_flatten(samples$observations$snapshots, c(3, 4))
  
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
