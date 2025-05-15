## code to prepare `PoO_df` dataset goes here

type <- 1:15
M <- c(rep(0, 3), rep(1, 2), 0, 2, rep(1, 5), rep(2, 3))
F <- c(0, rep(1, 2), rep(0, 2), 2, 0, rep(1, 3), rep(2, 2), rep(1, 2), 2)
C <- c(rep(0, 2), 1, 0, rep(1, 3), 0, 1, 2, 1, 2, 1, rep(2, 2))

## Indicator variables for parent-of-origin
matOrg <- c(rep(0, 4), 1, 0, 1, 0, NA, 1, 0, 1, 1, 1, 1)
patOrg <- c(rep(0, 2), 1, rep(0, 2), 1, 0, 0, NA, 1, 1, 1, 0, 1, 1)

PoO_df <- data.frame(type, M, F, C, matOrg, patOrg)

usethis::use_data(PoO_df, overwrite = TRUE)
