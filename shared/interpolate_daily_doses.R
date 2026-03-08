interpolate_daily_doses <- function(daily_doses_time, daily_doses_value) {
  t(vapply(seq_len(nrow(daily_doses_value)),
           function (i) {
             approx(daily_doses_time, 
                    daily_doses_value[i, ],
                    xout = seq(daily_doses_time[1],
                               daily_doses_time[length(daily_doses_time)]), 
                    method = "constant")$y},
           numeric(max(daily_doses_time) - min(daily_doses_time) + 1)))
}

uninterpolate_daily_doses <- function(daily_doses) {
  if (ncol(daily_doses) > 1) {
    i_change <- which(rowSums(diff(t(daily_doses)) != 0) > 0)
    i_change <- c(1, i_change + 1)
  } else {
    i_change <- 1
  }
  
  list(time = i_change,
       value = daily_doses[, i_change, drop = FALSE])
}
