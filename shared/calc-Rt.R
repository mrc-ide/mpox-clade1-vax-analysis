calc_Rt_fits <- function(samples, region) {
  pars <- lapply(seq_len(ncol(samples$pars)),
                 function (i) samples$packer$unpack(samples$pars[, i]))
  
  samples$Rt <- calc_Rt(samples$observations$trajectories, pars, region)
  
  samples$observations$trajectories <-
    reduce_trajectories(samples$observations$trajectories)
  
  samples             
}


calc_Rt_scenarios <- function(trajectories, pars, region) {
  traj <- aperm(trajectories, c(1, 3, 2))
  
  Rt <- t(calc_Rt(traj, pars, region, no_vax = FALSE))
  Rt_no_vax <- t(calc_Rt(traj, pars, region, no_vax = TRUE))
  
  out <- abind::abind(trajectories, Rt, Rt_no_vax, along = 1, make.names = TRUE)
  out <- reduce_trajectories(out)
  
  out
}


## This function takes a 3d array of trajectories (n_states x n_time x 
## n_samples) and a length n_samples list of lists of pars. It outputs a 2d
## array of Rt (n_time x n_samples)
calc_Rt <- function(trajectories, pars, region, no_vax = FALSE) {
  
  S_names <- get_S_names()
  
  if (region == "both") {
    S <- trajectories[S_names, , , ]
    regs <- dimnames(S)[[2]]
    
    Rt <-
      vapply(seq_along(pars),
             function (i) {
               vapply(regs, function (r) calc_Rt_1(pars[[i]][[r]], S[, r, , i]),
                      numeric(dim(S)[3]))},
             array(0, dim(S)[c(3, 2)]))
    Rt <- aperm(Rt, c(2, 1, 3))
    rownames(Rt) <- regs
    
  } else {
    S <- trajectories[S_names, , ]
    
    Rt <- vapply(seq_along(pars),
                 function (i) calc_Rt_1(pars[[i]], S[, , i]),
                 numeric(dim(S)[2]))
  }
  
  Rt
}


## This function takes in a single list of parameters and a 2d array of S
## trajectories (n_S x n_time) and outputs a vector of Rt (length n_time)
calc_Rt_1 <- function(pars, S, dt = 1, no_vax = FALSE) {
  N <- pars$N0
  
  ## get the individual-to-individual contact matrix for gen pop/sex contacts
  M_gen_pop <- pars$beta_h * 
    pars$m_gen_pop / rep(N, each = ncol(pars$m_gen_pop))
  M_sex <- pars$beta_s * pars$m_sex / rep(N, each = ncol(pars$m_gen_pop))
  
  idx <- mpoxseir::get_compartment_indices()
  
  ## Add together gen_pop and sex contacts
  M <- M_gen_pop + M_sex
  ## Add the extra contacts to HCW (not symmetric!)
  M[, idx$group$HCW] <- M[, idx$group$HCW] + pars$beta_hcw / sum(N)
  ## Expand to a 4x4 block matrix where each block is the same but corresponds
  ## to vaccine classes
  M <- block_expand(M, 4)
  
  ve_I <- pars$ve_I
  ve_T <- pars$ve_T
  
  if (no_vax) {
    # set VE for vaccinated to 0
    ve_I[, 3:4] <- 0
    ve_T[3:4] <- 0
  }
  
  k <- 1 # number of compartments per infectious class
  mean_duration <- pars$CFR * ((k * dt) / (1 - exp(-pars$gamma_Id * dt))) +
    (1 - pars$CFR) * ((k * dt) / (1 - exp(-pars$gamma_Ir * dt)))  * 
    rep(1 - ve_T, each = pars$n_group)
  
  ## This function will calculate Rt at one time-point
  compute_Rt <- function(x) {
    ## NGM[i, j] = M[i, j] * mean_duration[i] * S[j] * (1 - ve_I[j])
    ## (tcrossprod(x, y) for vectors creates a matrix with entries x[i] * y[j])
    NGM <- M * tcrossprod(c(mean_duration), x * c(1 - ve_I))
    j <- which(x != 0)
    eigen1::eigen1(NGM[j, j], max_iterations = 1e5, tolerance = 1e-6,
                   method = "power_iteration")
  }
  
  apply(S, 2, compute_Rt)
}


block_expand <- function(m, n) {
  if (n == 1L) {
    return(m)
  }
  len <- nrow(m) * n
  matrix(t(matrix(m, nrow(m), len)), len, len, byrow = TRUE)
}


get_S_names <- function() {
  compartment_indices <- mpoxseir::get_compartment_indices()
  groups <- seq_along(compartment_indices$group)
  vax <- seq_along(compartment_indices$vax)
  paste("S", groups, rep(vax, each = length(groups)), sep = "_")
}


## This function is for removing trajectories once used
## (e.g. S once Rt calculated)
reduce_trajectories <- function(trajectories) {
  S_names <- get_S_names()
  
  keep <- setdiff(rownames(trajectories), S_names)
  
  if (length(dim(trajectories)) == 4) {
    trajectories <- trajectories[keep, , , ]
  } else {
    trajectories <- trajectories[keep, , ]
  }
  
}

calc_R0 <- function(samples, region, no_historic_vax = FALSE) {
  pars <- lapply(seq_len(ncol(samples$pars)),
                 function (i) samples$packer$unpack(samples$pars[, i]))
  
  if (region == "both") {
    regs <- names(pars[[1]])
    
    R0 <-
      lapply(regs,
             function (r) {
               t(vapply(seq_along(pars), 
                        function (i) calc_R0_1(pars[[i]][[r]], 
                                               no_historic_vax = no_historic_vax),
                        numeric(3)))})
  } else {
    R0 <- t(vapply(seq_along(pars), 
                   function (i) calc_R0_1(pars[[i]], 
                                          no_historic_vax = no_historic_vax),
                   numeric(3)))
  }
  
  R0
}


calc_R0_1 <- function(pars, dt = 1, no_historic_vax = FALSE) {
  N <- pars$N0
  
  ## get the individual-to-individual contact matrix for gen pop/sex contacts
  M_gen_pop <- pars$beta_h * 
    pars$m_gen_pop / rep(N, each = ncol(pars$m_gen_pop))
  M_gen_pop[, N == 0] <- 0
  M_sex <- pars$beta_s * pars$m_sex / rep(N, each = ncol(pars$m_gen_pop))
  M_sex[, N == 0] <- 0
  
  idx <- mpoxseir::get_compartment_indices()
  
  ## Add together gen_pop and sex contacts
  M <- M_gen_pop + M_sex
  ## Add the extra contacts to HCW (not symmetric!)
  M[, idx$group$HCW] <- M[, idx$group$HCW] + pars$beta_hcw / sum(N)
  ## Expand to a 4x4 block matrix where each block is the same but corresponds
  ## to vaccine classes
  
  ve_I <- pars$ve_I
  ve_T <- pars$ve_T
  
  ## Should not be anyone in vax strata 3/4 but set their VE to 0 just in case
  ve_I[, 3:4] <- 0
  ve_T[3:4] <- 0
  
  if (no_historic_vax) {
    # set VE for historic vaccinated to 0
    ve_I[, 1] <- 0
    ve_T[1] <- 0
  }
  
  k <- 1 # number of compartments per infectious class
  mean_duration <- pars$CFR * ((k * dt) / (1 - exp(-pars$gamma_Id * dt))) +
    (1 - pars$CFR) * ((k * dt) / (1 - exp(-pars$gamma_Ir * dt)))  * 
    rep(1 - ve_T, each = pars$n_group)
  
  S0 <- c(pars$S0)
  
  ## This function will calculate Rt at one time-point
  compute_R0 <- function(MM) {
    ## NGM[i, j] = M[i, j] * mean_duration[i] * S[j] * (1 - ve_I[j])
    ## (tcrossprod(x, y) for vectors creates a matrix with entries x[i] * y[j])
    MM <- block_expand(MM, 4)
    NGM <- MM * tcrossprod(c(mean_duration), S0 * c(1 - ve_I))
    j <- which(colSums(NGM) != 0)
    max(Re(eigen(NGM[j, j])$values))
  }
  
  c(R0 = compute_R0(M), 
    R0_gen_pop = compute_R0(M_gen_pop),
    R0_sex = compute_R0(M_sex))
}