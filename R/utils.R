.wise_assert_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                closed = TRUE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) {
    ok <- if (closed) x >= lower && x <= upper else x > lower && x < upper
  }
  if (!ok) {
    stop(sprintf("`%s` must be one finite number in %s%g, %g%s.",
                 name, if (closed) "[" else "(", lower, upper,
                 if (closed) "]" else ")"), call. = FALSE)
  }
  invisible(x)
}

.wise_eval_parameter <- function(x, grid, name) {
  value <- if (is.function(x)) x(grid) else x
  value <- as.numeric(value)
  if (length(value) == 1L) value <- rep(value, length(grid))
  if (length(value) != length(grid) || any(!is.finite(value))) {
    stop(sprintf("`%s` must be a scalar, a function, or one value per grid node.",
                 name), call. = FALSE)
  }
  value
}

.wise_trapezoid <- function(y, x) {
  if (length(x) < 2L) return(0)
  n <- length(y)
  sum(diff(x) * (y[seq_len(n - 1L)] + y[seq.int(2L, n)]) / 2)
}

.wise_step_at <- function(time, value, xout, left = 1) {
  ord <- order(time)
  time <- as.numeric(time[ord])
  value <- as.numeric(value[ord])
  keep <- !duplicated(time, fromLast = TRUE)
  time <- time[keep]
  value <- value[keep]
  index <- findInterval(xout, time)
  out <- rep(left, length(xout))
  active <- index > 0L
  out[active] <- value[index[active]]
  out
}

.wise_arm <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    values <- unique(x)
    if (length(values) != 2L) stop("The trial must contain exactly two arms.", call. = FALSE)
    x <- match(x, values) - 1L
  }
  x <- as.integer(x)
  if (!all(x %in% c(0L, 1L))) {
    stop("The arm variable must contain only 0 and 1 (or two labels).", call. = FALSE)
  }
  x
}
