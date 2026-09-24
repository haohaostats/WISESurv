#' Fit a structural WISE-Surv model
#'
#' @param trial Two-arm trial data frame.
#' @param time,status,arm Column names for follow-up time, event indicator, and
#'   randomized arm.
#' @param cutoff Trial information cutoff in years.
#' @param horizon Maximum decision horizon in years.
#' @param external Optional `wise_external` object or list of such objects.
#' @param specification A `wise_spec` object.
#' @param grid_step Post-trial grid spacing in years.
#' @param exact_nodes Evaluation horizons inserted into the grid.
#'
#' @return An object of class `wise_surv_fit`.
#' @export
wise_surv <- function(trial, time = "time", status = "status", arm = "arm",
                      cutoff, horizon = 15, external = NULL,
                      specification = wise_spec(), grid_step = 1 / 12,
                      exact_nodes = c(7, 10, 15)) {
  if (!inherits(specification, "wise_spec")) {
    stop("`specification` must be created by wise_spec().", call. = FALSE)
  }
  trial <- .wise_prepare_trial(trial, time, status, arm, cutoff)
  if (is.null(external)) external <- list()
  if (inherits(external, "wise_external")) external <- list(external)
  if (!is.list(external) ||
      any(!vapply(external, inherits, logical(1), what = "wise_external"))) {
    stop("`external` must be a wise_external object or a list of them.",
         call. = FALSE)
  }
  grid <- wise_grid(cutoff, horizon, grid_step, exact_nodes)
  cutoff_survival <- vapply(0:1, function(a) {
    z <- trial[trial$arm == a, , drop = FALSE]
    .wise_km(z$time, z$status, cutoff)
  }, numeric(1))

  event_times <- sort(unique(trial$time[trial$status == 1L & trial$time <= cutoff]))
  observed_grid <- sort(unique(c(0, event_times, cutoff)))
  observed_survival <- lapply(0:1, function(a) {
    z <- trial[trial$arm == a, , drop = FALSE]
    .wise_km(z$time, z$status, observed_grid)
  })
  problem <- .wise_build_problem(grid, cutoff_survival, specification, external)
  feasible <- .wise_feasible(problem)

  structure(list(
    call = match.call(),
    trial = trial,
    cutoff = cutoff,
    horizon = horizon,
    grid = grid,
    cutoff_survival = cutoff_survival,
    observed_grid = observed_grid,
    observed_survival = observed_survival,
    external = external,
    specification = specification,
    problem = problem,
    feasible = feasible
  ), class = "wise_surv_fit")
}

#' @export
print.wise_surv_fit <- function(x, ...) {
  cat("WISE-Surv structural fit\n")
  cat(sprintf("  trial cutoff / horizon: %g / %g years\n", x$cutoff, x$horizon))
  cat(sprintf("  trial sample sizes: control=%d, intervention=%d\n",
              sum(x$trial$arm == 0L), sum(x$trial$arm == 1L)))
  cat(sprintf("  post-trial grid nodes: %d\n", length(x$grid)))
  cat(sprintf("  external sources: %d\n", length(x$external)))
  cat(sprintf("  feasible class: %s\n", if (x$feasible) "yes" else "no"))
  invisible(x)
}
