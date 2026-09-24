.wise_km_process <- function(time, status, grid) {
  time <- as.numeric(time)
  status <- as.integer(status)
  n <- length(time)
  m <- length(grid)
  event_times <- sort(unique(time[status == 1L]))
  estimate <- rep(1, m)
  raw_influence <- matrix(0, nrow = n, ncol = m)
  cloglog_influence <- matrix(0, nrow = n, ncol = m)
  cumulative <- numeric(n)
  survival <- 1
  cursor <- 1L

  for (j in seq_along(grid)) {
    while (cursor <= length(event_times) && event_times[cursor] <= grid[j] + 1e-12) {
      u <- event_times[cursor]
      at_risk <- time >= u - 1e-12
      died <- status == 1L & abs(time - u) < 1e-10
      y <- sum(at_risk)
      d <- sum(died)
      if (y > 0L && d > 0L) {
        cumulative <- cumulative +
          (as.numeric(died) - as.numeric(at_risk) * d / y) / (y / n)
        survival <- survival * (1 - d / y)
      }
      cursor <- cursor + 1L
    }
    estimate[j] <- survival
    if (survival > 1e-10 && survival < 1 - 1e-10) {
      raw_influence[, j] <- -survival * cumulative
      cloglog_influence[, j] <-
        raw_influence[, j] / (survival * log(survival))
    }
  }
  raw_influence <- sweep(raw_influence, 2L, colMeans(raw_influence), "-")
  cloglog_influence <- sweep(cloglog_influence, 2L,
                             colMeans(cloglog_influence), "-")
  sigma <- sqrt(colMeans(cloglog_influence^2))
  active <- estimate > 1e-10 & estimate < 1 - 1e-10 & sigma > 1e-12
  list(estimate = estimate, raw = raw_influence,
       cloglog = cloglog_influence, sigma = sigma,
       active = active, grid = grid, n = n)
}

.wise_multiplier_process <- function(component, multiplier) {
  active <- component$active
  if (!any(active)) return(matrix(0, nrow(multiplier), 0L))
  z <- multiplier %*% component$cloglog[, active, drop = FALSE]
  sweep(z / sqrt(component$n), 2L, component$sigma[active], "/")
}

.wise_survival_band <- function(component, critical) {
  lower <- upper <- component$estimate
  active <- component$active
  if (any(active)) {
    q <- log(-log(component$estimate[active]))
    half <- critical * component$sigma[active] / sqrt(component$n)
    lower[active] <- exp(-exp(q + half))
    upper[active] <- exp(-exp(q - half))
  }
  list(lower = pmax(0, lower), upper = pmin(1, upper),
       estimate = component$estimate, active = active,
       grid = component$grid)
}

#' Construct the joint simultaneous WISE-Surv input region
#'
#' @param object A fitted `wise_surv_fit` object.
#' @param level Simultaneous confidence level.
#' @param multipliers Number of Gaussian multiplier repetitions.
#' @param seed Random-number seed.
#'
#' @return An object of class `wise_input_region` containing arm-specific and
#'   external survival bands, randomized-arm difference constraints, and the
#'   observed-period RMST interval.
#' @export
wise_input_region <- function(object, level = 0.95, multipliers = 999L,
                              seed = 20260918L) {
  if (!inherits(object, "wise_surv_fit")) {
    stop("`object` must be a wise_surv_fit.", call. = FALSE)
  }
  .wise_assert_scalar(level, "level", 0, 1, closed = FALSE)
  multipliers <- as.integer(multipliers)
  if (length(multipliers) != 1L || is.na(multipliers) || multipliers < 2L) {
    stop("`multipliers` must be an integer of at least 2.", call. = FALSE)
  }
  set.seed(seed)
  pre_grid <- object$observed_grid
  trial_data <- lapply(0:1, function(a) object$trial[object$trial$arm == a, ])
  names(trial_data) <- c("control", "intervention")
  components <- lapply(trial_data, function(z) {
    .wise_km_process(z$time, z$status, pre_grid)
  })

  external_ipd <- Filter(function(z) identical(z$kind, "ipd"), object$external)
  external_components <- lapply(external_ipd, function(z) {
    .wise_km_process(z$time, z$status, object$grid)
  })
  external_names <- vapply(external_ipd, `[[`, character(1), "name")
  if (length(external_components)) names(external_components) <- external_names

  all_components <- c(components, external_components)
  multiplier_draws <- lapply(all_components, function(component) {
    matrix(stats::rnorm(multipliers * component$n),
           nrow = multipliers, ncol = component$n)
  })
  maxima <- numeric(multipliers)
  for (k in seq_along(all_components)) {
    z <- .wise_multiplier_process(all_components[[k]], multiplier_draws[[k]])
    if (ncol(z)) maxima <- pmax(maxima, apply(abs(z), 1L, max))
  }

  control <- components$control
  intervention <- components$intervention
  n0 <- control$n
  n1 <- intervention$n
  difference <- intervention$estimate - control$estimate
  difference_se <- sqrt(colMeans(intervention$raw^2) / n1 +
                          colMeans(control$raw^2) / n0)
  difference_active <- difference_se > 1e-12
  if (any(difference_active)) {
    z_difference <-
      (multiplier_draws$intervention %*%
         intervention$raw[, difference_active, drop = FALSE] / n1 -
       multiplier_draws$control %*%
         control$raw[, difference_active, drop = FALSE] / n0)
    z_difference <- sweep(z_difference, 2L,
                          difference_se[difference_active], "/")
    maxima <- pmax(maxima, apply(abs(z_difference), 1L, max))
  }

  pre_weights <- .wise_weights(pre_grid, object$cutoff)
  intervention_area_if <- as.numeric(intervention$raw %*% pre_weights)
  control_area_if <- as.numeric(control$raw %*% pre_weights)
  observed_area <- sum(pre_weights * difference)
  observed_area_se <- sqrt(mean(intervention_area_if^2) / n1 +
                             mean(control_area_if^2) / n0)
  if (observed_area_se > 1e-12) {
    z_area <-
      (as.numeric(multiplier_draws$intervention %*% intervention_area_if) / n1 -
       as.numeric(multiplier_draws$control %*% control_area_if) / n0) /
      observed_area_se
    maxima <- pmax(maxima, abs(z_area))
  }

  critical <- as.numeric(stats::quantile(maxima, probs = level,
                                         names = FALSE, type = 7))
  trial_bands <- lapply(components, .wise_survival_band, critical = critical)
  external_bands <- lapply(external_components, .wise_survival_band,
                           critical = critical)

  structure(list(
    level = level,
    repetitions = multipliers,
    seed = seed,
    critical = critical,
    trial = trial_bands,
    external = external_bands,
    difference = list(
      estimate = difference,
      lower = difference - critical * difference_se,
      upper = difference + critical * difference_se,
      se = difference_se,
      active = difference_active,
      grid = pre_grid
    ),
    observed_area = c(
      estimate = observed_area,
      lower = observed_area - critical * observed_area_se,
      upper = observed_area + critical * observed_area_se,
      se = observed_area_se
    )
  ), class = "wise_input_region")
}

#' @export
print.wise_input_region <- function(x, digits = 3, ...) {
  cat("WISE-Surv joint simultaneous input region\n")
  cat(sprintf("  level: %.1f%%\n", 100 * x$level))
  cat(sprintf("  Gaussian multipliers: %d\n", x$repetitions))
  cat(sprintf("  joint critical value: %.*f\n", digits, x$critical))
  invisible(x)
}
