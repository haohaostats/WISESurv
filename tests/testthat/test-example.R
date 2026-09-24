test_that("README example data are immediately usable", {
  trial <- wise_example_trial(n_per_arm = 50, seed = 20)
  expect_equal(names(trial), c("time", "status", "arm"))
  expect_equal(nrow(trial), 100)
  expect_true(all(trial$status %in% 0:1))

  fit <- wise_surv(
    trial,
    cutoff = 2,
    horizon = 15,
    specification = wise_spec(
      effect_upper = 1,
      residual = c(0, 1),
      decline_rate = 1,
      rebound_rate = 1
    )
  )
  expect_s3_class(fit, "wise_surv_fit")
  expect_true(fit$feasible)
  expect_s3_class(wise_bounds(fit, "rmst", c(7, 10, 15)), "wise_bounds")
})

test_that("README example gives a clear robust decision", {
  trial <- wise_example_trial()
  specification <- wise_spec(
    effect_upper = 0.35,
    residual = c(0.50, 0.80),
    decline_rate = 0.005,
    rebound_rate = 0.045
  )
  fit <- wise_surv(
    trial = trial,
    cutoff = 2,
    horizon = 15,
    specification = specification
  )
  qaly <- wise_outer_bounds(
    fit,
    estimand = "qaly",
    utility = 0.75,
    discount_rate = log(1.035),
    multipliers = 199,
    seed = 20260918
  )
  decision <- wise_nmb(
    qaly,
    cost = 20000,
    willingness_to_pay = c(20000, 30000, 50000)
  )

  expect_true(all(qaly$summary$status == "feasible"))
  expect_true(all(qaly$summary$lower > 1))
  expect_true(all(decision$decision == "robust adoption"))

  figure <- tempfile(fileext = ".pdf")
  grDevices::pdf(figure)
  expect_silent(plot(qaly))
  grDevices::dev.off()
  expect_true(file.exists(figure))
})

test_that("evidence overview plot accepts external IPD", {
  trial <- wise_example_trial(n_per_arm = 80, seed = 31)
  external_data <- subset(
    wise_example_trial(n_per_arm = 100, seed = 32,
                       administrative_censoring = 10),
    arm == 0,
    select = c(time, status)
  )
  source <- wise_external_ipd(
    external_data,
    arm = 0,
    tolerance = 0.20,
    active_range = c(2, 10),
    name = "External control"
  )
  fit <- wise_surv(
    trial,
    cutoff = 2,
    horizon = 15,
    external = source,
    specification = wise_spec(effect_upper = 1, residual = c(0, 1),
                              decline_rate = 1, rebound_rate = 1)
  )
  input <- wise_input_region(fit, multipliers = 49, seed = 33)
  figure <- tempfile(fileext = ".pdf")
  grDevices::pdf(figure, width = 10, height = 4)
  expect_silent(plot(fit, input_region = input))
  grDevices::dev.off()
  expect_true(file.exists(figure))
})
