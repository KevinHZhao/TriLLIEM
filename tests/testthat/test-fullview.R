test_that("expected example_dat4R", {
  set.seed(1)
  expect_identical(fullview(simulateData()), fullview(example_dat4R))
})

test_that("expected control and environmental data", {
  set.seed(2)
  expect_snapshot(
    simulateData(
      nCases = 1234,
      nControl = 2345,
      mtCoef = c(1.2, 1.2, 1.2),
      propE = 0.3
    ) |>
      fullview()
  )
})
