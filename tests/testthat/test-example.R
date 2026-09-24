test_that("README example data are immediately usable", {
  trial <- wise_example_trial(n_per_arm = 50, seed = 20)
  expect_equal(names(trial), c("time", "status", "arm"))
  expect_equal(nrow(trial), 100)
  expect_true(all(trial$status %in% 0:1))

  fit <- wise_surv(trial, cutoff = 2, horizon = 15)
  expect_s3_class(fit, "wise_surv_fit")
  expect_true(fit$feasible)
  expect_s3_class(wise_bounds(fit, "rmst", c(7, 10, 15)), "wise_bounds")
})
