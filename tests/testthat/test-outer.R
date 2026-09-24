test_that("joint input region includes trial contrast and observed area", {
  trial <- example_trial(n = 50, seed = 10)
  fit <- wise_surv(trial, cutoff = 2)
  region <- wise_input_region(fit, multipliers = 99, seed = 10)
  expect_s3_class(region, "wise_input_region")
  expect_true(is.finite(region$critical))
  expect_equal(length(region$difference$estimate), length(fit$observed_grid))
  expect_true(all(is.finite(region$observed_area)))
})

test_that("outer bounds contain their optimized endpoints in order", {
  trial <- example_trial(n = 60, seed = 11)
  fit <- wise_surv(trial, cutoff = 2)
  structural <- wise_bounds(fit, "rmst", horizons = c(7, 15))
  outer <- wise_outer_bounds(fit, "rmst", horizons = c(7, 15),
                             multipliers = 99, seed = 11)
  expect_s3_class(outer, "wise_outer_bounds")
  expect_true(all(outer$summary$status == "feasible"))
  expect_true(all(outer$summary$lower <= outer$summary$upper))
  expect_true(all(outer$summary$lower <= structural$summary$lower + 1e-8))
  expect_true(all(outer$summary$upper >= structural$summary$upper - 1e-8))
})

test_that("external IPD can constrain one trial arm", {
  trial <- example_trial(n = 50, seed = 12)
  external <- example_trial(n = 100, seed = 13)
  external <- external[external$arm == 0, c("time", "status")]
  source <- wise_external_ipd(external, arm = 0, tolerance = .15,
                              active_range = c(2, 7), name = "control source")
  fit <- wise_surv(trial, cutoff = 2, external = source)
  expect_s3_class(fit, "wise_surv_fit")
  region <- wise_input_region(fit, multipliers = 49, seed = 12)
  expect_true("control source" %in% names(region$external))
})
