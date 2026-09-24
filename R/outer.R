.wise_interp_band <- function(band, grid) {
  list(
    lower = .wise_step_at(band$grid, band$lower, grid,
                          left = band$lower[1L]),
    upper = .wise_step_at(band$grid, band$upper, grid,
                          left = band$upper[1L])
  )
}

.wise_external_outer_intervals <- function(object, input, grid) {
  out <- vector("list", length(object$external))
  for (k in seq_along(object$external)) {
    source <- object$external[[k]]
    if (identical(source$kind, "ipd")) {
      band <- input$external[[source$name]]
      source_lower <- .wise_step_at(band$grid, band$lower, grid,
                                    left = band$lower[1L])
      source_upper <- .wise_step_at(band$grid, band$upper, grid,
                                    left = band$upper[1L])
    } else {
      source_estimate <- .wise_step_at(source$time, source$survival, grid,
                                       left = source$survival[1L])
      source_lower <- source_upper <- source_estimate
    }
    tolerance <- if (is.function(source$tolerance)) source$tolerance(grid) else
      source$tolerance
    if (is.matrix(tolerance) || is.data.frame(tolerance)) {
      tolerance <- as.matrix(tolerance)
      lower_tolerance <- tolerance[, 1L]
      upper_tolerance <- tolerance[, 2L]
    } else {
      tolerance <- as.numeric(tolerance)
      if (length(tolerance) == 1L) tolerance <- rep(tolerance, 2L)
      if (length(tolerance) == 2L) {
        lower_tolerance <- rep(tolerance[1L], length(grid))
        upper_tolerance <- rep(tolerance[2L], length(grid))
      } else {
        lower_tolerance <- upper_tolerance <- tolerance
      }
    }
    active <- grid >= source$active_range[1L] - 1e-10 &
      grid <= source$active_range[2L] + 1e-10
    lower <- upper <- rep(NA_real_, length(grid))
    lower[active] <- pmax(0, source_lower[active] - lower_tolerance[active])
    upper[active] <- pmin(1, source_upper[active] + upper_tolerance[active])
    out[[k]] <- list(arm = source$arm, lower = lower, upper = upper,
                     name = source$name)
  }
  out
}

.wise_build_outer_problem <- function(object, input) {
  grid <- sort(unique(c(object$observed_grid, object$grid)))
  n <- length(grid)
  i0 <- function(j) j
  i1 <- function(j) n + j
  rows <- .wise_new_rows(2L * n)
  cutoff_index <- which.min(abs(grid - object$cutoff))
  post <- which(grid >= object$cutoff - 1e-10)

  for (offset in c(0L, n)) {
    for (j in seq.int(2L, n)) {
      .wise_add_row(rows, c(offset + j, offset + j - 1L), c(1, -1), 0,
                    group = "survival")
    }
  }

  specification <- object$specification
  if (specification$non_harm) {
    for (j in post) {
      .wise_add_row(rows, c(i0(j), i1(j)), c(1, -1), 0,
                    group = "non-harm")
    }
  }
  envelope <- .wise_eval_parameter(specification$effect_upper, grid,
                                   "effect_upper")
  for (j in post) {
    .wise_add_row(rows, c(i1(j), i0(j)), c(1, -1), envelope[j],
                  group = "envelope")
  }
  end_index <- length(grid)
  .wise_add_row(
    rows,
    c(i1(end_index), i0(end_index), i1(cutoff_index), i0(cutoff_index)),
    c(1, -1, -specification$residual_upper,
      specification$residual_upper), 0, group = "residual"
  )
  .wise_add_row(
    rows,
    c(i1(end_index), i0(end_index), i1(cutoff_index), i0(cutoff_index)),
    c(-1, 1, specification$residual_lower,
      -specification$residual_lower), 0, group = "residual"
  )
  for (j in post[post > cutoff_index]) {
    width <- grid[j] - grid[j - 1L]
    columns <- c(i1(j), i0(j), i1(j - 1L), i0(j - 1L))
    .wise_add_row(rows, columns, c(1, -1, -1, 1),
                  specification$rebound_rate * width,
                  group = "rate", scale = width)
    .wise_add_row(rows, columns, c(-1, 1, 1, -1),
                  specification$decline_rate * width,
                  group = "rate", scale = width)
  }

  trial_bands <- lapply(input$trial, .wise_interp_band, grid = grid)
  pre <- which(grid <= object$cutoff + 1e-10)
  for (j in pre) {
    .wise_add_row(rows, i0(j), 1, trial_bands$control$upper[j],
                  group = "trial-band")
    .wise_add_row(rows, i0(j), -1, -trial_bands$control$lower[j],
                  group = "trial-band")
    .wise_add_row(rows, i1(j), 1, trial_bands$intervention$upper[j],
                  group = "trial-band")
    .wise_add_row(rows, i1(j), -1, -trial_bands$intervention$lower[j],
                  group = "trial-band")
  }
  difference_lower <- .wise_step_at(input$difference$grid,
                                    input$difference$lower, grid,
                                    left = input$difference$lower[1L])
  difference_upper <- .wise_step_at(input$difference$grid,
                                    input$difference$upper, grid,
                                    left = input$difference$upper[1L])
  for (j in pre) {
    .wise_add_row(rows, c(i1(j), i0(j)), c(1, -1),
                  difference_upper[j], group = "trial-difference")
    .wise_add_row(rows, c(i1(j), i0(j)), c(-1, 1),
                  -difference_lower[j], group = "trial-difference")
  }

  pre_weights <- .wise_weights(grid, object$cutoff)
  active <- which(pre_weights != 0)
  .wise_add_row(rows, c(n + active, active),
                c(pre_weights[active], -pre_weights[active]),
                input$observed_area["upper"], group = "trial-area")
  .wise_add_row(rows, c(n + active, active),
                c(-pre_weights[active], pre_weights[active]),
                -input$observed_area["lower"], group = "trial-area")

  external_intervals <- .wise_external_outer_intervals(object, input, grid)
  for (source in external_intervals) {
    offset <- source$arm * n
    active <- which(is.finite(source$lower) & is.finite(source$upper) &
                      grid >= object$cutoff - 1e-10)
    for (j in active) {
      group <- paste0("external:", source$name)
      .wise_add_row(rows, offset + j, 1, source$upper[j], group = group)
      .wise_add_row(rows, offset + j, -1, -source$lower[j], group = group)
    }
  }

  list(
    A = .wise_rows_matrix(rows), rhs = rows$rhs,
    lower = rep(0, 2L * n), upper = rep(1, 2L * n),
    group = rows$group, scale = rows$scale,
    grid = grid, external = external_intervals
  )
}

#' Compute sampling-error-protected outer bounds
#'
#' @param object A fitted `wise_surv_fit` object.
#' @param estimand Either `"rmst"` or `"qaly"`.
#' @param horizons Evaluation horizons.
#' @param utility Utility while alive for QALYs.
#' @param discount_rate Continuous annual discount rate.
#' @param level Simultaneous confidence level.
#' @param multipliers Number of Gaussian multiplier repetitions.
#' @param seed Random-number seed.
#'
#' @return A `wise_outer_bounds` object, also inheriting from `wise_bounds`.
#' @export
wise_outer_bounds <- function(object, estimand = c("rmst", "qaly"),
                              horizons = object$horizon, utility = 1,
                              discount_rate = 0, level = 0.95,
                              multipliers = 999L, seed = 20260918L) {
  if (!inherits(object, "wise_surv_fit")) {
    stop("`object` must be a wise_surv_fit.", call. = FALSE)
  }
  estimand <- match.arg(estimand)
  input <- wise_input_region(object, level, multipliers, seed)
  problem <- .wise_build_outer_problem(object, input)
  utility_function <- if (is.function(utility)) utility else {
    value <- as.numeric(utility)
    function(t) rep(value, length(t))
  }
  node_weight <- if (estimand == "rmst") NULL else function(t) {
    utility_function(t) * exp(-discount_rate * t)
  }
  rows <- vector("list", length(horizons))
  paths <- vector("list", length(horizons))
  for (k in seq_along(horizons)) {
    h <- horizons[k]
    if (!any(abs(problem$grid - h) < 1e-9)) {
      stop("Every horizon must be an exact node of the complete grid.", call. = FALSE)
    }
    weights <- .wise_weights(problem$grid, h, node_weight)
    solution <- .wise_solve_bounds(problem, c(-weights, weights), constant = 0)
    if (is.null(solution)) {
      rows[[k]] <- data.frame(estimand = estimand, horizon = h,
                              status = "conflict", lower = NA_real_,
                              upper = NA_real_, width = NA_real_)
    } else {
      rows[[k]] <- data.frame(estimand = estimand, horizon = h,
                              status = "feasible", lower = solution$lower,
                              upper = solution$upper,
                              width = solution$upper - solution$lower)
      paths[[k]] <- solution
    }
  }
  out <- list(summary = do.call(rbind, rows), paths = paths,
              fit = object, estimand = estimand, utility = utility,
              discount_rate = discount_rate, input_region = input,
              problem = problem)
  class(out) <- c("wise_outer_bounds", "wise_bounds")
  out
}

#' @export
print.wise_outer_bounds <- function(x, digits = 3, ...) {
  cat("WISE-Surv sampling-error-protected outer bounds\n")
  cat(sprintf("  simultaneous level: %.1f%%\n", 100 * x$input_region$level))
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}
