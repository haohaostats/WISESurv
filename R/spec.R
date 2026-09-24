#' Define a WISE-Surv scientific specification
#'
#' @param effect_upper Scalar, vector, or function giving the upper envelope
#'   for the marginal survival difference.
#' @param residual Length-two vector giving the lower and upper fractions of
#'   cutoff benefit allowed to remain at the decision horizon.
#' @param decline_rate Maximum annual decline in marginal survival difference.
#' @param rebound_rate Maximum annual increase in marginal survival difference.
#' @param non_harm Whether to impose nonnegative marginal survival benefit.
#'
#' @return An object of class `wise_spec`.
#' @export
wise_spec <- function(effect_upper = 0.20, residual = c(0, 0.80),
                      decline_rate = 0.035, rebound_rate = 0.045,
                      non_harm = TRUE) {
  if (!is.numeric(residual) || length(residual) != 2L ||
      any(!is.finite(residual)) || residual[1L] < 0 ||
      residual[2L] > 1 || residual[1L] > residual[2L]) {
    stop("`residual` must satisfy 0 <= residual[1] <= residual[2] <= 1.",
         call. = FALSE)
  }
  .wise_assert_scalar(decline_rate, "decline_rate", 0, Inf)
  .wise_assert_scalar(rebound_rate, "rebound_rate", 0, Inf)
  if (!is.logical(non_harm) || length(non_harm) != 1L || is.na(non_harm)) {
    stop("`non_harm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!(is.function(effect_upper) || is.numeric(effect_upper))) {
    stop("`effect_upper` must be numeric or a function.", call. = FALSE)
  }
  structure(list(
    effect_upper = effect_upper,
    residual_lower = residual[1L],
    residual_upper = residual[2L],
    decline_rate = decline_rate,
    rebound_rate = rebound_rate,
    non_harm = non_harm
  ), class = "wise_spec")
}

#' @export
print.wise_spec <- function(x, ...) {
  cat("WISE-Surv scientific specification\n")
  cat(sprintf("  marginal non-harm: %s\n", x$non_harm))
  cat(sprintf("  residual fraction: [%g, %g]\n",
              x$residual_lower, x$residual_upper))
  cat(sprintf("  decline/rebound rates: %g / %g per year\n",
              x$decline_rate, x$rebound_rate))
  cat("  effect upper envelope: ",
      if (is.function(x$effect_upper)) "function\n" else
        paste(x$effect_upper, collapse = ", "), "\n", sep = "")
  invisible(x)
}
