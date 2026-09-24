.wise_observed_contribution <- function(object, node_weight = NULL) {
  grid <- object$observed_grid
  weight <- if (is.null(node_weight)) rep(1, length(grid)) else node_weight(grid)
  d <- object$observed_survival[[2L]] - object$observed_survival[[1L]]
  .wise_trapezoid(weight * d, grid)
}

#' Compute sharp structural bounds
#'
#' @param object A fitted `wise_surv_fit` object.
#' @param estimand Either `"rmst"` or `"qaly"`.
#' @param horizons Evaluation horizons. Each must exceed the trial cutoff.
#' @param utility Utility while alive, as a scalar or function of time.
#' @param discount_rate Continuous annual discount rate used for QALYs.
#'
#' @return A `wise_bounds` object containing a summary table and extremizing
#'   trajectories.
#' @export
wise_bounds <- function(object, estimand = c("rmst", "qaly"),
                        horizons = object$horizon, utility = 1,
                        discount_rate = 0) {
  if (!inherits(object, "wise_surv_fit")) {
    stop("`object` must be a wise_surv_fit.", call. = FALSE)
  }
  estimand <- match.arg(estimand)
  horizons <- as.numeric(horizons)
  if (any(!is.finite(horizons)) || any(horizons <= object$cutoff) ||
      any(horizons > object$horizon)) {
    stop("All horizons must be greater than cutoff and no greater than the fitted horizon.",
         call. = FALSE)
  }
  missing_nodes <- horizons[!vapply(horizons, function(h)
    any(abs(object$grid - h) < 1e-9), logical(1))]
  if (length(missing_nodes)) {
    stop("Refit with these horizons in `exact_nodes`: ",
         paste(missing_nodes, collapse = ", "), call. = FALSE)
  }

  utility_function <- if (is.function(utility)) utility else {
    value <- as.numeric(utility)
    if (length(value) != 1L || !is.finite(value)) {
      stop("`utility` must be one finite number or a function.", call. = FALSE)
    }
    function(t) rep(value, length(t))
  }
  node_weight <- if (estimand == "rmst") NULL else function(t) {
    utility_function(t) * exp(-discount_rate * t)
  }

  rows <- vector("list", length(horizons))
  paths <- vector("list", length(horizons))
  for (k in seq_along(horizons)) {
    h <- horizons[k]
    if (!object$feasible) {
      rows[[k]] <- data.frame(estimand = estimand, horizon = h,
                              status = "conflict", lower = NA_real_,
                              upper = NA_real_, width = NA_real_)
      next
    }
    weights <- .wise_weights(object$grid, h, node_weight)
    objective <- c(-weights, weights)
    constant <- .wise_observed_contribution(object, node_weight)
    solution <- .wise_solve_bounds(object$problem, objective, constant)
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
  structure(list(summary = do.call(rbind, rows), paths = paths,
                 fit = object, estimand = estimand,
                 utility = utility, discount_rate = discount_rate),
            class = "wise_bounds")
}

#' @export
print.wise_bounds <- function(x, digits = 3, ...) {
  cat("WISE-Surv sharp structural bounds\n")
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
summary.wise_bounds <- function(object, ...) object$summary

#' Plot extremizing trajectories
#'
#' @param x A `wise_bounds` object.
#' @param which Row of the bounds summary to plot.
#' @param type Plot marginal survival-effect trajectories (`"effect"`) or the
#'   four arm-specific extremizing survival trajectories (`"survival"`).
#' @param main Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param ... Additional graphical arguments passed to `plot`.
#' @export
plot.wise_bounds <- function(x, which = 1L,
                             type = c("effect", "survival"), main = NULL,
                             xlab = "Time (years)",
                             ylab = NULL, ...) {
  which <- as.integer(which)
  type <- match.arg(type)
  if (length(which) != 1L || which < 1L || which > nrow(x$summary)) {
    stop("`which` must identify one row of the bounds summary.", call. = FALSE)
  }
  solution <- x$paths[[which]]
  if (is.null(solution)) stop("No trajectories are available for a conflict.", call. = FALSE)
  grid <- if (!is.null(x$problem$grid)) x$problem$grid else x$fit$grid
  n <- length(grid)
  horizon <- x$summary$horizon[which]
  active <- which(grid <= horizon + 1e-10)
  if (is.null(main)) {
    main <- sprintf("%g-year %s extremizing %s trajectories",
                    horizon, toupper(x$estimand),
                    if (type == "effect") "benefit" else "survival")
  }
  if (type == "effect") {
    lower <- solution$lower_path[n + active] -
      solution$lower_path[active]
    upper <- solution$upper_path[n + active] -
      solution$upper_path[active]
    ylim <- range(c(lower, upper), finite = TRUE)
    padding <- max(diff(ylim) * 0.12, 0.005)
    ylim <- c(ylim[1L] - padding, ylim[2L] + padding)
    if (is.null(ylab)) ylab <- "Marginal survival benefit"
    graphics::plot(grid[active], lower, type = "l", col = "#007C83",
                   lwd = 2.4, ylim = ylim, xlab = xlab, ylab = ylab,
                   main = main, bty = "l", las = 1, xaxs = "i", ...)
    graphics::lines(grid[active], upper, col = "#D97706", lwd = 2.4,
                    lty = 2)
    if (ylim[1L] <= 0 && ylim[2L] >= 0) {
      graphics::abline(h = 0, col = "#AAB2BD", lty = 3)
    }
    graphics::legend(
      "bottomleft", legend = c("Lower-bound path", "Upper-bound path"),
      col = c("#007C83", "#D97706"), lty = c(1, 2), lwd = 2.4,
      bty = "n"
    )
  } else {
    if (is.null(ylab)) ylab <- "Survival probability"
    graphics::plot(grid[active], solution$lower_path[active], type = "l",
                   col = "#007C83", lwd = 2, ylim = c(0, 1),
                   xlab = xlab, ylab = ylab, main = main,
                   bty = "l", las = 1, xaxs = "i", yaxs = "i", ...)
    graphics::lines(grid[active], solution$lower_path[n + active],
                    col = "#007C83", lwd = 2, lty = 2)
    graphics::lines(grid[active], solution$upper_path[active],
                    col = "#D97706", lwd = 2)
    graphics::lines(grid[active], solution$upper_path[n + active],
                    col = "#D97706", lwd = 2, lty = 2)
    graphics::legend(
      "bottomleft",
      legend = c("Lower: control", "Lower: intervention",
                 "Upper: control", "Upper: intervention"),
      col = c("#007C83", "#007C83", "#D97706", "#D97706"),
      lty = c(1, 2, 1, 2), lwd = 2, bty = "n", ncol = 2
    )
  }
  if (min(grid[active]) < x$fit$cutoff - 1e-10) {
    graphics::abline(v = x$fit$cutoff, col = "#7A8490", lty = 3)
  }
  invisible(x)
}
