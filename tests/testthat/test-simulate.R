test_that("expected example_dat4R", {
  set.seed(1)
  expect_identical(simulateData(), example_dat4R)
})

test_that("expected control and environmental data", {
  set.seed(2)
  expect_snapshot(
    simulateData(
      nCases = 1234,
      nControl = 2345,
      mtCoef = c(1.2, 1.2, 1.2),
      propE = 0.3
    )
  )
})

test_that("expected population stratification data", {
  set.seed(3)
  dat <-
    simulateData(
      nCases = 1234,
      nControl = 2345,
      propE = c(0.3, 0.5, 0.7),
      nPop = 3,
      maf = c(0.1, 0.2, 0.3),
      prev.byPop = c(0.1, 0.15, 0.2),
      prop.byPop = c(0.3, 0.3, 0.4)
    )
  expect_snapshot(dat)
})
