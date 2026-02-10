test_that("summary of default (C, M) model works", {
  res <- TriLLIEM(dat = example_dat4R)
  expect_snapshot(summary(res))
})

test_that("summary of control and environment data model works", {
  set.seed(2)
  dat <-
    simulateData(
      nCases = 1234,
      nControl = 2345,
      mtCoef = c(1.2, 1.2, 1.2),
      propE = 0.3
    )
  res <-
    TriLLIEM(
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
  expect_snapshot(summary(res))
  expect_snapshot(summary(resStrat))
})

test_that("population stratification model works", {
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
  resHWE <-
    TriLLIEM(
      mtmodel = "HWE",
      effects = c("C", "M", "Im", "E:Im"),
      dat = dat,
      includeE = TRUE,
      includeD = TRUE
    )
  resMS <-
    TriLLIEM(
      mtmodel = "MS",
      effects = c("C", "M", "Im", "E:Im"),
      dat = dat,
      includeE = TRUE,
      includeD = TRUE
    )
  resMaS <-
    TriLLIEM(
      mtmodel = "MaS",
      effects = c("C", "M", "Im", "E:Im"),
      dat = dat,
      includeE = TRUE,
      includeD = TRUE
    )
  expect_snapshot(summary(resHWE))
  expect_snapshot(summary(resMS))
  expect_snapshot(summary(resMaS))
})
