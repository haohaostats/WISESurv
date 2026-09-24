#' Diagnose an incompatible WISE-Surv class
#'
#' @param object A fitted `wise_surv_fit` object.
#'
#' @return A list containing the total and group-specific minimum relaxation.
#' @export
wise_diagnose <- function(object) {
  if (!inherits(object, "wise_surv_fit")) {
    stop("`object` must be a wise_surv_fit.", call. = FALSE)
  }
  problem <- object$problem
  relaxed_groups <- setdiff(unique(problem$group), "survival")
  if (!length(relaxed_groups)) {
    return(structure(list(status = if (object$feasible) "feasible" else "conflict",
                          total = 0, by_group = numeric()),
                     class = "wise_diagnosis"))
  }
  n_original <- length(problem$lower)
  group_index <- structure(seq_along(relaxed_groups), names = relaxed_groups)
  slack_column <- n_original + group_index[problem$group]
  relaxable <- problem$group != "survival"
  slack_i <- which(relaxable)
  slack_j <- slack_column[relaxable]
  slack_x <- -problem$scale[relaxable]
  A <- cbind(problem$A, Matrix::Matrix(0, nrow(problem$A), length(relaxed_groups),
                                      sparse = TRUE))
  A[cbind(slack_i, slack_j)] <- slack_x
  objective <- c(rep(0, n_original), rep(1, length(relaxed_groups)))
  fit <- highs::highs_solve(
    L = objective,
    lower = c(problem$lower, rep(0, length(relaxed_groups))),
    upper = c(problem$upper, rep(Inf, length(relaxed_groups))),
    A = A,
    lhs = rep(-Inf, length(problem$rhs)),
    rhs = problem$rhs,
    control = highs::highs_control(threads = 1L, log_to_console = FALSE,
                                   solver = "simplex")
  )
  if (!identical(fit$status_message, "Optimal")) {
    stop("Minimum-relaxation problem failed: ", fit$status_message,
         call. = FALSE)
  }
  by_group <- fit$primal_solution[n_original + seq_along(relaxed_groups)]
  names(by_group) <- relaxed_groups
  structure(list(status = if (object$feasible) "feasible" else "conflict",
                 total = sum(by_group), by_group = by_group),
            class = "wise_diagnosis")
}

#' @export
print.wise_diagnosis <- function(x, digits = 4, ...) {
  cat("WISE-Surv compatibility diagnosis\n")
  cat("  fitted class: ", x$status, "\n", sep = "")
  cat("  minimum total relaxation: ", format(round(x$total, digits)), "\n",
      sep = "")
  if (length(x$by_group)) print(round(x$by_group, digits))
  invisible(x)
}
