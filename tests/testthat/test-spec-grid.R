test_that("specification validates scientific parameters", {
  expect_s3_class(wise_spec(), "wise_spec")
  expect_error(wise_spec(residual = c(.8, .2)))
  expect_error(wise_spec(decline_rate = -1))
})

test_that("grid includes requested horizons exactly", {
  grid <- wise_grid(2.1, 15, exact_nodes = c(7, 10, 15))
  expect_true(all(c(7, 10, 15) %in% grid))
  expect_equal(grid[1], 2.1)
  expect_equal(tail(grid, 1), 15)
})
