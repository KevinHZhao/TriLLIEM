## code to prepare `example_dat4R` dataset goes here

example_dat4R <- withr::with_seed(
  seed = 1,
  simulateData()$dat4R
)

usethis::use_data(example_dat4R, overwrite = TRUE)
