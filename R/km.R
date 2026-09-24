.wise_km <- function(time, status, grid) {
  time <- as.numeric(time)
  status <- as.integer(status)
  event_times <- sort(unique(time[status == 1L]))
  survival <- 1
  estimate <- rep(1, length(grid))
  cursor <- 1L
  for (j in seq_along(grid)) {
    while (cursor <= length(event_times) && event_times[cursor] <= grid[j] + 1e-12) {
      u <- event_times[cursor]
      risk <- sum(time >= u - 1e-12)
      deaths <- sum(status == 1L & abs(time - u) < 1e-10)
      if (risk > 0L && deaths > 0L) survival <- survival * (1 - deaths / risk)
      cursor <- cursor + 1L
    }
    estimate[j] <- survival
  }
  estimate
}

.wise_prepare_trial <- function(data, time, status, arm, cutoff) {
  if (!is.data.frame(data)) stop("`trial` must be a data frame.", call. = FALSE)
  missing_names <- setdiff(c(time, status, arm), names(data))
  if (length(missing_names)) {
    stop("Missing trial columns: ", paste(missing_names, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    time = as.numeric(data[[time]]),
    status = as.integer(data[[status]]),
    arm = .wise_arm(data[[arm]])
  )
  if (any(!is.finite(out$time)) || any(out$time < 0) ||
      any(!out$status %in% c(0L, 1L))) {
    stop("Trial times must be nonnegative and status must be 0/1.", call. = FALSE)
  }
  original_time <- out$time
  out$time <- pmin(out$time, cutoff)
  out$status <- as.integer(out$status == 1L & original_time <= cutoff + 1e-12)
  out
}
