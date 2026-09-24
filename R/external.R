#' Define interval-valued external survival evidence
#'
#' @param time Time coordinates in years.
#' @param survival External survival estimates.
#' @param arm Trial arm constrained by the source (0 or 1).
#' @param tolerance Absolute transportability tolerance. May be a scalar,
#'   length-two vector (lower, upper), or a function returning one or two
#'   columns.
#' @param active_range Time range over which the source is imposed.
#' @param name Source label used in diagnostics.
#'
#' @return An object of class `wise_external`.
#' @export
wise_external <- function(time, survival, arm = 0L, tolerance = 0,
                          active_range = range(time), name = "external") {
  time <- as.numeric(time)
  survival <- as.numeric(survival)
  if (length(time) != length(survival) || length(time) < 1L ||
      any(!is.finite(time)) || any(!is.finite(survival)) ||
      any(survival < 0 | survival > 1)) {
    stop("`time` and `survival` must be equal-length finite vectors, with survival in [0, 1].",
         call. = FALSE)
  }
  arm <- as.integer(arm)
  if (length(arm) != 1L || !arm %in% c(0L, 1L)) {
    stop("`arm` must be 0 or 1.", call. = FALSE)
  }
  ord <- order(time)
  time <- time[ord]
  survival <- survival[ord]
  if (any(diff(survival) > 1e-8)) {
    stop("External survival must be nonincreasing.", call. = FALSE)
  }
  if (!is.numeric(active_range) || length(active_range) != 2L ||
      active_range[1L] > active_range[2L]) {
    stop("`active_range` must be a length-two increasing vector.", call. = FALSE)
  }
  if (!(is.function(tolerance) || is.numeric(tolerance))) {
    stop("`tolerance` must be numeric or a function.", call. = FALSE)
  }
  structure(list(time = time, survival = survival, arm = arm,
                 tolerance = tolerance, active_range = active_range,
                 name = as.character(name)[1L], kind = "curve"),
            class = "wise_external")
}

#' Define external individual-level survival evidence
#'
#' @param data External right-censored survival data.
#' @param time,status Column names for follow-up time and event indicator.
#' @param arm Trial arm constrained by the source (0 or 1).
#' @param tolerance Absolute transportability tolerance; see [wise_external()].
#' @param active_range Time range over which the source is imposed.
#' @param name Source label used in diagnostics.
#'
#' @return An object of class `wise_external`.
#' @export
wise_external_ipd <- function(data, time = "time", status = "status", arm = 0L,
                              tolerance = 0, active_range = c(0, Inf),
                              name = "external") {
  if (!is.data.frame(data) || !all(c(time, status) %in% names(data))) {
    stop("`data` must contain the requested time and status columns.",
         call. = FALSE)
  }
  observed_time <- as.numeric(data[[time]])
  observed_status <- as.integer(data[[status]])
  if (any(!is.finite(observed_time)) || any(observed_time < 0) ||
      any(!observed_status %in% c(0L, 1L))) {
    stop("External times must be nonnegative and status must be 0/1.",
         call. = FALSE)
  }
  arm <- as.integer(arm)
  if (length(arm) != 1L || !arm %in% c(0L, 1L)) {
    stop("`arm` must be 0 or 1.", call. = FALSE)
  }
  if (!(is.function(tolerance) || is.numeric(tolerance))) {
    stop("`tolerance` must be numeric or a function.", call. = FALSE)
  }
  structure(list(time = observed_time, status = observed_status, arm = arm,
                 tolerance = tolerance, active_range = active_range,
                 name = as.character(name)[1L], kind = "ipd"),
            class = "wise_external")
}

.wise_external_interval <- function(source, grid) {
  estimate <- if (identical(source$kind, "ipd")) {
    .wise_km(source$time, source$status, grid)
  } else {
    .wise_step_at(source$time, source$survival, grid,
                  left = source$survival[1L])
  }
  tolerance <- if (is.function(source$tolerance)) source$tolerance(grid) else
    source$tolerance
  if (is.matrix(tolerance) || is.data.frame(tolerance)) {
    tolerance <- as.matrix(tolerance)
    if (nrow(tolerance) != length(grid) || ncol(tolerance) != 2L) {
      stop("A tolerance function returning a matrix must have two columns.",
           call. = FALSE)
    }
    lower_tolerance <- tolerance[, 1L]
    upper_tolerance <- tolerance[, 2L]
  } else {
    tolerance <- as.numeric(tolerance)
    if (length(tolerance) == 1L) tolerance <- rep(tolerance, 2L)
    if (length(tolerance) == 2L) {
      lower_tolerance <- rep(tolerance[1L], length(grid))
      upper_tolerance <- rep(tolerance[2L], length(grid))
    } else if (length(tolerance) == length(grid)) {
      lower_tolerance <- upper_tolerance <- tolerance
    } else {
      stop("Tolerance must be scalar, length two, one value per node, or a two-column matrix.",
           call. = FALSE)
    }
  }
  if (any(!is.finite(c(lower_tolerance, upper_tolerance))) ||
      any(lower_tolerance < 0 | upper_tolerance < 0)) {
    stop("Transportability tolerances must be finite and nonnegative.", call. = FALSE)
  }
  active <- grid >= source$active_range[1L] - 1e-10 &
    grid <= source$active_range[2L] + 1e-10
  lower <- upper <- rep(NA_real_, length(grid))
  lower[active] <- pmax(0, estimate[active] - lower_tolerance[active])
  upper[active] <- pmin(1, estimate[active] + upper_tolerance[active])
  list(arm = source$arm, lower = lower, upper = upper, name = source$name,
       estimate = estimate)
}
