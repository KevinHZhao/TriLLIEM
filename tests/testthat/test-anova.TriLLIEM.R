test_that("simple example_dat4R anova works", {
  model_1 <- TriLLIEM::TriLLIEM(dat = example_dat4R, effects = c("C", "M"))
  model_2 <- TriLLIEM::TriLLIEM(dat = example_dat4R, effects = c("C", "M", "Im"))
  expect_snapshot(anova(model_1, model_2))
})

test_that("Eanova works", {
  set.seed(2)
  dat <-
    simulateData(
      nCases = 1234,
      nControl = 2345,
      mtCoef = c(1.2, 1.2, 1.2),
      propE = 0.3
    )
  res1 <-
    TriLLIEM(
      mtmodel = "HWE",
      effects = c("C", "M", "Im", "E:Im"),
      dat = dat,
      includeE = TRUE,
      includeD = TRUE
    )
  res2 <-
    TriLLIEM(
      mtmodel = "HWE",
      effects = c("C", "M", "Im"),
      dat = dat,
      Eanova = TRUE,
      includeD = TRUE
    )
  expect_snapshot(anova(res1, res2))
})
