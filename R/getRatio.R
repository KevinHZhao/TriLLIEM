## Get ratio of cases/controls by population. Used only under
## population stratification
getRatio <- function(alpha, maf, omega, R = c(1, 1, 1), S = c(1, 1, 1), V = c(1, 1, 1),
                     mtcoef = c(1, 1, 1)) {
  numPop <- length(alpha)
  K <- vector(length = numPop)

  if ((length(maf) != numPop) || (length(omega) != numPop)) {
    stop("Number of elements of maf, alpha and omega must equal the
         number of subpopulations\n")
  }
  if (!is.list(mtcoef)) {
    mtcoef <- list(mtcoef)
    mtcoef <- rep(mtcoef, numPop)
  }

  browser()
  for (i in 1:numPop) {
    genomat <- mtmat(maf[i], C = mtcoef[[i]])
    diseasefactor <- R[(genomat$C + 1)] * S[(genomat$M + 1)]
    diseaseprob <- alpha[i] * genomat$prMFC * diseasefactor
    K[i] <- sum(diseaseprob)
  }

  prPopGivenCase <- K * omega / sum(K * omega)

  prPopGivenControl <- (1 - K) * omega / sum((1 - K) * omega)


  return(list(
    prevalence = K, prPopGivenCase = prPopGivenCase,
    prPopGivenControl = prPopGivenControl
  ))
}
