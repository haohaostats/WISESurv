#' Generate a reproducible example two-arm survival trial
#'
#' Creates a small synthetic right-censored trial for learning the package API
#' and checking an installation. The generator is not part of the scientific
#' simulation design used to evaluate WISE-Surv.
#'
#' @param n_per_arm Number of participants in each randomized arm.
#' @param seed Random-number seed.
#' @param administrative_censoring Maximum follow-up time in years.
#'
#' @return A data frame with columns `time`, `status`, and `arm`.
#' @export
wise_example_trial <- function(n_per_arm = 300L, seed = 20260924L,
                               administrative_censoring = 5) {
  n_per_arm <- as.integer(n_per_arm)
  if (length(n_per_arm) != 1L || is.na(n_per_arm) || n_per_arm < 10L) {
    stop("`n_per_arm` must be an integer of at least 10.", call. = FALSE)
  }
  .wise_assert_scalar(administrative_censoring, "administrative_censoring",
                      0, Inf, closed = FALSE)
  set.seed(seed)
  arm <- rep(0:1, each = n_per_arm)
  event_rate <- ifelse(arm == 1L, 0.11, 0.18)
  event_time <- stats::rexp(2L * n_per_arm, rate = event_rate)
  dropout_time <- stats::rexp(2L * n_per_arm, rate = 0.03)
  censor_time <- pmin(dropout_time, administrative_censoring)
  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    arm = arm
  )
}
