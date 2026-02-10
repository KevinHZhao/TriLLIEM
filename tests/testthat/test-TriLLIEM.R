test_that("default (C, M) model works", {
  res <- TriLLIEM(dat = example_dat4R)
  expect_snapshot(res)
})

test_that("control and environment data model works",{
  set.seed(2)
  dat <-
    simulateData(
      nCases = 1234,
      nControl = 2345,
      mtCoef = c(1.2, 1.2, 1.2),
      propE = 0.3
    )
  res <-
    TriLLIEM::TriLLIEM(
      mtmodel = "MaS",
      effects = c("C", "M", "Im", "E:M"),
      dat = dat,
      includeE = TRUE,
      includeD = TRUE
    )
  resStrat <-
    TriLLIEM(
      mtmodel = "MaS",
      effects = c("C", "M", "Im", "E:M"),
      dat = dat,
      includeE = TRUE,
      Estrat = TRUE,
      includeD = TRUE
    )
  expect_snapshot(res)
  expect_snapshot(resStrat)
})
