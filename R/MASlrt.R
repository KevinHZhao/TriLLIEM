#' Mating asymmetry likelihood ratio test
#'
#' `MaSlrt()` runs a Likelihood Ratio Test for a trios model with Mating Asymmetry.
#'
#' @param dat A dataframe formatted with the same columns as [example_dat4R]
#'
#' @return A p-value for the Likelihood Ratio Test.
#' @noRd
#'
#' @examples
#' MaSlrt(example_dat4R)
MaSlrt <- function(dat) {
  subdat <- dat[c(2:7, 11:14), ]
  subdat <- aggregate(subdat$count, by = list(M = subdat$M, F = subdat$F), sum)
  subdat <- subdat[c(1, 3, 2, 5, 4, 6), ]
  colnames(subdat)[3] <- "count"

  num <- 3 * (subdat$count[5] + subdat$count[6]) + 2 * (subdat$count[3] + subdat$count[4]) +
    (subdat$count[1] + subdat$count[2])
  denom <- 4 * sum(subdat$count)
  p <- num / denom

  ishet.M <- subdat$M == 1
  ishet.F <- subdat$F == 1
  factor.M <- rep(1, nrow(subdat))
  factor.M[ishet.M] <- 2
  factor.F <- rep(1, nrow(subdat))
  factor.F[ishet.F] <- 2
  HWfreqs <- factor.M * p^subdat$M * (1 - p)^(2 - subdat$M) * factor.F * p^subdat$F * (1 - p)^(2 - subdat$F)


  c4 <- 2 * subdat$count[5] / (subdat$count[5] + subdat$count[6])
  c2 <- 2 * subdat$count[3] / (subdat$count[3] + subdat$count[4])
  c1 <- 2 * subdat$count[1] / (subdat$count[1] + subdat$count[2])
  Cfactor <- c(c1, 2 - c1, c2, 2 - c2, c4, 2 - c4)

  l1 <- sum(subdat$count * log(HWfreqs * Cfactor))
  l0 <- sum(subdat$count * log(HWfreqs))
  LRT <- 2 * (l1 - l0)
  return(pchisq(LRT, df = 3, lower = F))
}
