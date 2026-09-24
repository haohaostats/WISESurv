.wise_panel_title <- function(label, title) {
  graphics::mtext(label, side = 3, adj = 0, line = 0.7,
                  font = 2, cex = 1.05)
  graphics::title(main = title, font.main = 2, cex.main = 1.02,
                  line = 0.5)
}

.wise_band_polygon <- function(time, lower, upper, colour) {
  keep <- is.finite(time) & is.finite(lower) & is.finite(upper)
  if (sum(keep) < 2L) return(invisible(NULL))
  time <- time[keep]
  lower <- lower[keep]
  upper <- upper[keep]
  graphics::polygon(c(time, rev(time)), c(lower, rev(upper)),
                    col = colour, border = NA)
  invisible(NULL)
}

#' Plot the evidence entering a WISE-Surv analysis
#'
#' Produces a publication-style three-panel diagnostic showing randomized
#' trial survival, external comparator evidence, and the trial-identified
#' marginal survival difference with its simultaneous confidence band.
#'
#' @param x A fitted `wise_surv_fit` object.
#' @param input_region Optional result from [wise_input_region()]. If omitted,
#'   it is computed from `x`.
#' @param level Simultaneous confidence level used when `input_region` is not
#'   supplied.
#' @param multipliers Number of Gaussian multiplier repetitions.
#' @param seed Random-number seed for the multiplier calculation.
#' @param external Index of the external source to display.
#' @param ... Reserved for future graphical options.
#'
#' @return `x`, invisibly.
#' @export
plot.wise_surv_fit <- function(x, input_region = NULL, level = 0.95,
                               multipliers = 999L, seed = 20260918L,
                               external = 1L, ...) {
  if (is.null(input_region)) {
    input_region <- wise_input_region(x, level = level,
                                      multipliers = multipliers, seed = seed)
  }
  if (!inherits(input_region, "wise_input_region")) {
    stop("`input_region` must be created by wise_input_region().",
         call. = FALSE)
  }
  external <- as.integer(external)
  if (length(external) != 1L || is.na(external) || external < 1L) {
    stop("`external` must be one positive integer.", call. = FALSE)
  }

  teal <- "#008C8C"
  blue <- "#4771B2"
  brown <- "#A66A23"
  dark <- "#34495E"
  sampling_colour <- grDevices::adjustcolor("#AFC9DC", alpha.f = 0.55)
  transport_colour <- grDevices::adjustcolor("#EFCF72", alpha.f = 0.55)
  difference_colour <- grDevices::adjustcolor("#8FD0C8", alpha.f = 0.48)

  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2L, byrow = TRUE),
                   heights = c(1, 0.18))
  graphics::par(mar = c(4.2, 4.2, 3.0, 0.9), mgp = c(2.35, 0.7, 0),
                tcl = -0.25, las = 1, family = "sans")

  trial_grid <- x$observed_grid
  trial_values <- unlist(x$observed_survival, use.names = FALSE)
  trial_ylim <- range(c(trial_values, 1), finite = TRUE)
  trial_ylim[1L] <- max(0, trial_ylim[1L] - 0.05)
  graphics::plot(NA, xlim = c(0, x$cutoff), ylim = trial_ylim,
                 xlab = "Time since randomization, years",
                 ylab = "Overall survival", bty = "l", xaxs = "i")
  graphics::lines(trial_grid, x$observed_survival[[1L]], type = "s",
                  col = blue, lwd = 2)
  graphics::lines(trial_grid, x$observed_survival[[2L]], type = "s",
                  col = teal, lwd = 2)
  graphics::abline(v = x$cutoff, col = "#6F7A86", lty = 3)
  .wise_panel_title("A", "Randomized trial evidence")

  if (length(x$external) >= external) {
    source <- x$external[[external]]
    upper_time <- min(x$horizon, max(source$time))
    external_grid <- sort(unique(c(0, source$time[source$time <= upper_time],
                                   x$cutoff, upper_time)))
    interval <- .wise_external_interval(source, external_grid)
    external_ylim <- range(c(interval$lower, interval$upper,
                             interval$estimate), na.rm = TRUE)
    external_ylim <- pmax(0, pmin(1, external_ylim + c(-0.04, 0.04)))
    graphics::plot(NA, xlim = c(0, upper_time), ylim = external_ylim,
                   xlab = "Time, years", ylab = "Overall survival",
                   bty = "l", xaxs = "i")
    .wise_band_polygon(external_grid, interval$lower, interval$upper,
                       transport_colour)
    if (identical(source$kind, "ipd") &&
        source$name %in% names(input_region$external)) {
      band <- input_region$external[[source$name]]
      .wise_band_polygon(band$grid, band$lower, band$upper,
                         sampling_colour)
    }
    graphics::lines(external_grid, interval$estimate, type = "s",
                    col = brown, lwd = 2)
    graphics::abline(v = x$cutoff, col = "#6F7A86", lty = 3)
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No external evidence supplied",
                   col = "#6F7A86")
  }
  .wise_panel_title("B", "External comparator evidence")

  difference <- input_region$difference
  difference_ylim <- range(c(0, difference$lower, difference$upper),
                           finite = TRUE)
  difference_ylim <- difference_ylim + c(-0.03, 0.03)
  graphics::plot(NA, xlim = c(0, x$cutoff), ylim = difference_ylim,
                 xlab = "Time, years", ylab = expression(D(t) == S[1](t)-S[0](t)),
                 bty = "l", xaxs = "i")
  .wise_band_polygon(difference$grid, difference$lower, difference$upper,
                     difference_colour)
  graphics::lines(difference$grid, difference$estimate, type = "s",
                  col = dark, lwd = 2)
  graphics::abline(h = 0, col = "#6F7A86", lty = 3)
  graphics::abline(v = x$cutoff, col = "#6F7A86", lty = 3)
  .wise_panel_title("C", "Trial-identified difference")

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = c("Intervention KM", "Control KM", "External KM",
               "External sampling band", "Transport interval",
               "Simultaneous difference band"),
    col = c(teal, blue, brown, sampling_colour, transport_colour,
            difference_colour),
    lty = c(1, 1, 1, NA, NA, NA),
    lwd = c(2, 2, 2, NA, NA, NA),
    pch = c(NA, NA, NA, 15, 15, 15),
    pt.cex = c(NA, NA, NA, 2, 2, 2),
    ncol = 3, bty = "n", cex = 0.86, xpd = NA
  )
  invisible(x)
}
