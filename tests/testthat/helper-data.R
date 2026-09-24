example_trial <- function(n = 80, seed = 1) {
  set.seed(seed)
  arm <- rep(0:1, each = n)
  event_time <- rexp(2 * n, rate = ifelse(arm == 1, 0.12, 0.18))
  censor_time <- runif(2 * n, 2, 7)
  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    arm = arm
  )
}
