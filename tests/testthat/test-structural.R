test_that("structural bounds are feasible and ordered", {
  trial <- example_trial()
  fit <- wise_surv(trial, cutoff = 2, exact_nodes = c(7, 10, 15))
  expect_s3_class(fit, "wise_surv_fit")
  expect_true(fit$feasible)

  out <- wise_bounds(fit, "rmst", c(7, 10, 15))
  expect_s3_class(out, "wise_bounds")
  expect_true(all(out$summary$status == "feasible"))
  expect_true(all(out$summary$lower <= out$summary$upper))
  expect_true(all(out$summary$width >= 0))
})

test_that("stronger restrictions cannot widen structural bounds", {
  trial <- example_trial(seed = 2)
  broad <- wise_surv(trial, cutoff = 2,
                     specification = wise_spec(residual = c(0, .8)))
  narrow <- wise_surv(trial, cutoff = 2,
                      specification = wise_spec(residual = c(.5, .8),
                                                rebound_rate = 0))
  b <- wise_bounds(broad, "rmst", 15)$summary
  n <- wise_bounds(narrow, "rmst", 15)$summary
  if (n$status == "feasible") expect_lte(n$width, b$width + 1e-8)
})

test_that("fixed-cost NMB rule returns the three-way classification", {
  trial <- example_trial(seed = 3)
  fit <- wise_surv(trial, cutoff = 2)
  q <- wise_bounds(fit, "qaly", 15, utility = .75,
                   discount_rate = log(1.035))
  decision <- wise_nmb(q, cost = 20000, willingness_to_pay = c(0, 30000))
  expect_equal(nrow(decision), 2)
  expect_true(all(decision$decision %in%
                    c("robust adoption", "robust rejection",
                      "decision not identified", "conflict")))
})
