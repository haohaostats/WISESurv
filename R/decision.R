#' Compute robust net monetary benefit decisions
#'
#' @param object A QALY `wise_bounds` object or a `wise_surv_fit` object.
#' @param cost Fixed incremental cost.
#' @param willingness_to_pay One or more willingness-to-pay values.
#' @param horizon QALY horizon when `object` is a fit.
#' @param utility Utility while alive when `object` is a fit.
#' @param discount_rate Continuous annual discount rate.
#'
#' @return A data frame with NMB bounds and decision classifications.
#' @export
wise_nmb <- function(object, cost, willingness_to_pay,
                     horizon = NULL, utility = 1, discount_rate = 0) {
  .wise_assert_scalar(cost, "cost")
  willingness_to_pay <- as.numeric(willingness_to_pay)
  if (!length(willingness_to_pay) || any(!is.finite(willingness_to_pay)) ||
      any(willingness_to_pay < 0)) {
    stop("`willingness_to_pay` must contain finite nonnegative values.",
         call. = FALSE)
  }
  if (inherits(object, "wise_surv_fit")) {
    if (is.null(horizon)) horizon <- object$horizon
    object <- wise_bounds(object, "qaly", horizon, utility, discount_rate)
  }
  if (!inherits(object, "wise_bounds") || object$estimand != "qaly") {
    stop("`object` must be a wise_surv_fit or QALY wise_bounds object.",
         call. = FALSE)
  }
  if (nrow(object$summary) != 1L) {
    stop("Supply QALY bounds for exactly one horizon.", call. = FALSE)
  }
  lower <- object$summary$lower
  upper <- object$summary$upper
  nmb_lower <- willingness_to_pay * lower - cost
  nmb_upper <- willingness_to_pay * upper - cost
  decision <- ifelse(
    is.na(nmb_lower), "conflict",
    ifelse(nmb_lower > 0, "robust adoption",
           ifelse(nmb_upper < 0, "robust rejection",
                  "decision not identified"))
  )
  data.frame(
    horizon = object$summary$horizon,
    willingness_to_pay = willingness_to_pay,
    incremental_cost = cost,
    nmb_lower = nmb_lower,
    nmb_upper = nmb_upper,
    decision = decision
  )
}
