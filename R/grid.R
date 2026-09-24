#' Construct the post-trial analysis grid
#'
#' @param cutoff Trial follow-up cutoff in years.
#' @param horizon Decision horizon in years.
#' @param step Regular grid spacing in years.
#' @param exact_nodes Additional nodes to insert exactly.
#'
#' @return A strictly increasing numeric vector.
#' @export
wise_grid <- function(cutoff, horizon = 15, step = 1 / 12,
                      exact_nodes = c(7, 10, 15)) {
  .wise_assert_scalar(cutoff, "cutoff", 0, Inf)
  .wise_assert_scalar(horizon, "horizon", cutoff, Inf, closed = FALSE)
  .wise_assert_scalar(step, "step", 0, Inf, closed = FALSE)
  n_regular <- floor((horizon - cutoff) / step + 1e-10)
  regular <- cutoff + step * seq.int(0, n_regular)
  exact_nodes <- exact_nodes[is.finite(exact_nodes) &
                               exact_nodes >= cutoff & exact_nodes <= horizon]
  grid <- sort(unique(round(c(regular, exact_nodes, horizon), 12L)))
  grid[grid >= cutoff - 1e-10 & grid <= horizon + 1e-10]
}

.wise_weights <- function(grid, endpoint, node_weight = NULL) {
  if (!any(abs(grid - endpoint) < 1e-9)) {
    stop("The endpoint must be present as an exact grid node.", call. = FALSE)
  }
  z <- if (is.null(node_weight)) rep(1, length(grid)) else node_weight(grid)
  z <- as.numeric(z)
  if (length(z) == 1L) z <- rep(z, length(grid))
  if (length(z) != length(grid) || any(!is.finite(z))) {
    stop("The weight function must return one finite value per grid node.",
         call. = FALSE)
  }
  out <- numeric(length(grid))
  for (j in seq.int(2L, length(grid))) {
    left <- grid[j - 1L]
    right <- min(grid[j], endpoint)
    if (right <= left) next
    width <- right - left
    out[j - 1L] <- out[j - 1L] + width * z[j - 1L] / 2
    out[j] <- out[j] + width * z[j] / 2
    if (grid[j] >= endpoint - 1e-10) break
  }
  out
}
