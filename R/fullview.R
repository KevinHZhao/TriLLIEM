#' View underlying simulation details for `TriLLIEM.sim` objects
#'
#' @param sim an object of class `TriLLIEM.sim`
#'
#' @returns A data frame with the same counts in each category as `sim`, but
#' with underlying simulation probabilities given as columns:
#' \describe{
#'  \item{"`typeOrig`"}{Index for each genetically distinct row,
#'  based on "`M`", "`F`", "`C`", "`matOrg`", and "`patOrg`".
#'  Ranges from 1 to 64, unlike the `type` column in `sim` which ranges from 1
#'  to 60, as the inclusion of "`matOrg`" and "`patOrg`" lead to potentially
#'  four additional rows.}
#'  \item{"`matOrg`" and "`patOrg`"}{Binary variables for when the minor allele
#'  is maternally/paternally inherited, so that counts when "`M`", "`F`", "`C`"
#'  are all 1 are not ambiguous when determining parent-of-origin effects.}
#'  \item{"`prMF`"}{Probability of the mother ("`M`") and father ("`F`") pairing
#'  in the simulated population, based on the minor allele frequency and mating
#'  type coefficients given during simulation.}
#'  \item{"`prCGivenMFOrg`"}{Probability of the child's genotype ("`C`") and the
#'  allele parent of origin ("`matOrg`" and "`patOrg`") conditional on the
#'  genotypes of the mother and father.}
#'  \item{"`prMFCOrg`"}{Probability of the triad based on genotypes and parent
#'  of origin, product of "`prMF`" and "`prCGivenMFOrg`".}
#'  \item{"`PrMFCOrg`"}{Probability of the triad based on genotypes, parent of
#'  origin, and environmental exposure ("`E`") conditional on disease status
#'  ("`D`"), obtained by scaling "`prMFCOrg`" by the risk of disease (equal to
#'  "`prMFCOrg`" when "`D`" is 0).}
#'  \item{"`pop`"}{Sub-population each row is simulated from (all `1` if
#' simulated without population stratification)}
#' }
#'
#' The "`count`" column in "`sim`" is the sum of the rows of the "`count`"
#' column in the returned data frame when grouped by "`M`", "`F`", "`C`", "`E`",
#' "`D`".
#' @export
#'
#' @examples
#'
#' ## View the underlying distributions behind some of the example models given in simulateData()
#' dat1 <- simulateData(S = c(1, 2, 4), If = 3)
#' fullview(dat1)
#'
#' dat2 <- simulateData(
#'   nControl = 1000,
#'   propE = c(0.1, 0.4),
#'   propE.control = c(0.2, 0.2),
#'   nPop = 2,
#'   maf = c(0.3, 0.4),
#'   prev.byPop = c(0.2, 0.3),
#'   prop.byPop = c(0.6, 0.4)
#'  )
#' fullview(dat2)
fullview <- function (sim) {
  if (!("TriLLIEM.sim" %in% class(sim))){
    stop("sim is not an object of class TriLLIEM.sim.")
  }
  attr(sim, "ALL")
}
