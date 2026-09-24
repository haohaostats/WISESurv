.wise_solve <- function(problem, objective, maximum = FALSE) {
  fit <- highs::highs_solve(
    L = as.numeric(objective),
    lower = problem$lower,
    upper = problem$upper,
    A = problem$A,
    lhs = rep(-Inf, length(problem$rhs)),
    rhs = problem$rhs,
    maximum = maximum,
    control = highs::highs_control(
      threads = 1L,
      log_to_console = FALSE,
      primal_feasibility_tolerance = 1e-9,
      dual_feasibility_tolerance = 1e-9,
      solver = "simplex",
      simplex_strategy = 1L
    )
  )
  optimal <- identical(fit$status_message, "Optimal") &&
    isTRUE(fit$solver_msg$value_valid)
  if (!optimal) return(NULL)
  x <- fit$primal_solution
  violation <- max(c(0, as.numeric(problem$A %*% x - problem$rhs),
                     problem$lower - x, x - problem$upper), na.rm = TRUE)
  if (violation > 1e-7) {
    stop(sprintf("HiGHS solution violates a constraint by %.3e.", violation),
         call. = FALSE)
  }
  list(value = as.numeric(fit$objective_value), x = x,
       violation = violation, raw = fit)
}

.wise_feasible <- function(problem) {
  !is.null(.wise_solve(problem, numeric(length(problem$lower))))
}

.wise_solve_bounds <- function(problem, objective, constant = 0) {
  lower <- .wise_solve(problem, objective, maximum = FALSE)
  upper <- .wise_solve(problem, objective, maximum = TRUE)
  if (is.null(lower) || is.null(upper)) return(NULL)
  list(
    lower = constant + lower$value,
    upper = constant + upper$value,
    lower_path = lower$x,
    upper_path = upper$x,
    max_violation = max(lower$violation, upper$violation)
  )
}
