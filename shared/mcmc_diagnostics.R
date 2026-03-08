calc_diagnostics <- function(samples) {
  samples$pars <- samples$pars_full
  samples$pars <- samples$pars[, -seq_len(samples$control$pmcmc$n_burnin), ]
  samples_df <- posterior::as_draws_df(samples)
  
  posterior::summarise_draws(samples_df)
}
