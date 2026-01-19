#' Simulate data for TriLLIEM
#'
#' Function to simulate maternal, paternal, and child genotype counts under
#' different genetic effect models.
#'
#' @param nCases number of case trios to simulate
#' @param nControl number of control trios
#' @param R vector of 3 elements representing child effects for 0, 1, and 2
#' copies of the risk allele, respectively
#' @param S vector of 3 elements representing maternal effects for 0, 1, and 2
#' copies of the risk allele, respectively
#' @param V vector of 3 elements representing gene-environment effects
#' for 0, 1, and 2 copies of the risk allele, respectively
#' @param mtCoef vector of 3 elements representing the mating type coefficients,
#' see details for more information
#' @param Im maternal imprinting effect
#' @param If paternal imprinting effect
#' @param propE scalar or vector of proportion of case trios in the
#' environmental exposure group in a single or multiple populations (when
#' `nPop > 1`)
#' @param propE.control scalar or vector of proportion of control trios in the
#' environmental exposure group in a single or multiple populations (when
#' `nPop > 1`)
#' @param Einteraction string indicating what variable environmental effects
#' interact with (one of "`Im`", "`If`", "`C`", or "`M`")
#' @param nPop number of populations to simulate for population
#' stratification, see details for information on
#' specifying other parameters when there are multiple strata
#' @param maf scalar or vector of the minor allele frequency proportion
#' in a single or multiple populations (when `nPop > 1`)
#' @param prev.byPop prevalence of cases in each sub population (vector with
#' length equal to `nPop`, ignored if `nPop = 1`)
#' @param prop.byPop proportion of each sub population, must sum to 1 (vector
#' with length equal to `nPop`, ignored if `nPop = 1`)
#' @param Fst F parameter for estimation of `maf` when sub-population level
#' `maf`'s are unknown (ignored when `nPop = 1` or when `maf` has the same
#' length as `nPop`), see details for more information
#'
#' @details
#' To simulate the counts, first the total number of case trios
#' (`nCases`) and control trios (`nControl`) are partitioned into the `nPop`
#' sub-populations by randomly sampling from respective multinomial
#' distributions, with probabilities calculated using `prop.byPop` and
#' `prev.byPop` as follows:
#' \deqn{\Pr(X = i, E = 1 \mid D = 1) = \frac{\alpha_i \varepsilon_{i, 1} \omega_i}{\alpha_1 \omega_1 + \alpha_2 \omega_2},}
#' \deqn{\Pr(X = i, E = 0 \mid D = 1) = \frac{\alpha_i (1 - \varepsilon_{i, 1}) \omega_i}{\alpha_1 \omega_1 + \alpha_2 \omega_2},}
#' \deqn{\Pr(X = i, E = 1 \mid D = 0) = \frac{(1 - \alpha_i) \varepsilon_{i, 0} \omega_i}{(1 - \alpha_1) \omega_1 + (1 - \alpha_2) \omega_2},}
#' \deqn{\Pr(X = i, E = 0 \mid D = 0) = \frac{(1 - \alpha_i) (1 - \varepsilon_{i, 0}) \omega_i}{(1 - \alpha_1) \omega_1 + (1 - \alpha_2) \omega_2},}
#' where \eqn{X = i} is the \eqn{i}-th sub-population, \eqn{\alpha_i} is the
#' \eqn{i}-th element of `prev.byPop`, \eqn{\varepsilon_{i,j}} \eqn{i}-th element
#' of `propE` (when \eqn{j=1}) or `propE.control` (when \eqn{j=0}), and
#' \eqn{\omega_i} is the \eqn{i}-th element of `prop.byPop.`
#'
#' If only a single aggregated `maf` value is provided when there are multiple
#' sub-populations, `maf` may be estimated using the model from
#' \insertCite{BaldNich95;textual}{TriLLIEM}.  This requires an \eqn{F_{ST}} parameter to be provided,
#' with estimates of `maf` being obtained by sampling from a beta distribution
#' with shape parameters
#' \deqn{\alpha = \text{MAF} \times \left(\frac{1}{F_{ST}} - 1\right),}
#' \deqn{\beta = (1 - \text{MAF}) \times \left(\frac{1}{F_{ST}} - 1\right),}
#' where MAF is the provided `maf`.
#'
#'  If `nPop = 1`, this step is skipped as everyone will be in the
#'  same sub-population.  Note that this method of simulating population
#'  stratification assumes \eqn{D} and \eqn{E} are conditionally independent
#'  (i.e., `V = c(1,1,1)`).
#'
#'  Distributions of counts in each partition based on the different genotype,
#'  exposure, and disease status categories are then independently sampled
#'  based on the risk of disease (from the given `R`, `S`, `V`, `Im`, and `If`
#'  values) and mating type probabilities (from `mtCoef`).
#'
#'  `R`, `S`, and `V` are vectors of 3 elements, where the first element should
#'  be 1, acting as a baseline for when \eqn{M}, \eqn{C}, or \eqn{M\times E} are
#'  0 respectively.  The following two numbers represent the multiplicative
#'  increase to risk when the respective genotype is 1 or 2.  `Im` and `If` are
#'  the imprinting parameters, representing the multiplicative increase to risk
#'  when the minor allele is maternally or paternally inherited.
#'
#'  Mating type probabilities are simulated using the mating asymmetry
#'  parameterization of \insertCite{Bour+11;textual}{TriLLIEM}.
#'  Here, `mtCoef` is a vector of 3 elements representing \eqn{(C_1, C_2, C_4)}.
#'  We begin by assigning each
#'  \eqn{\Pr(M, F)} assuming the population is in Hardy-Weinberg equilibrium,
#'  with probabilities obtained from the given `maf`.  Letting
#'  \deqn{\mu_1 = \Pr(M = 2, F = 1) = \Pr(M = 1, F = 2),}
#'  \deqn{\mu_2 = \Pr(M = 2, F = 0) = \Pr(M = 0, F = 2),}
#'  \deqn{\mu_4 = \Pr(M = 1, F = 0) = \Pr(M = 0, F = 1),}
#'  we reassign these particular mating pairs by setting
#'  \deqn{\Pr(M = 2, F = 1) = \mu_1C_1,}
#'  \deqn{\Pr(M = 1, F = 2) = \mu_1\left(2 - C_1\right),}
#'  \deqn{\Pr(M = 2, F = 0) = \mu_2C_2,}
#'  \deqn{\Pr(M = 0, F = 2) = \mu_2\left(2 - C_2\right),}
#'  \deqn{\Pr(M = 1, F = 0) = \mu_4C_4,}
#'  \deqn{\Pr(M = 0, F = 1) = \mu_4\left(2 - C_4\right).}
#'
#'  To view the data frame with underlying simulation probabilities, counts for
#'  either parent of origin case when mother, father, and child are all,
#'  heterozygous, and by sub-population when simulating population
#'  stratification, use [fullview()].
#'
#' @return An object of class `TriLLIEM.sim`, a data frame with columns for:
#' \describe{
#'    \item{"`type`"}{Index for the category corresponding to maternal ("`M`"),
#'    paternal ("`F`"), and child ("`C`") genotypes.}
#'    \item{"`mt_MS`"}{Mating type category for maternal and paternal
#'    genotypes under a mating symmetry model.}
#'    \item{"`mt_MaS`"}{Mating type category for maternal and paternal
#'    genotypes under a model that does not assume mating symmetry.}
#'    \item{"`M`"}{Maternal genotype.}
#'    \item{"`F`"}{Paternal genotype.}
#'    \item{"`C`"}{Child genotype.}
#'    \item{"`E`"}{Binary variable for if environmental effects are present (1) or
#'    not present (0).}
#'    \item{"`D`"}{Case child (1)/Control child (0).}
#'    \item{"`count`"}{Number of triads falling under the specified category based
#'    on "`M`", "`F`", "`C`", "`E`", "`D`".}
#' }
#'
#' @references{
#' \insertAllCited{}
#' }
#'
#' @seealso [example_dat4R], [fullview()]
#'
#' @export
#'
#' @examples
#' ## Maternal effect of 2, and paternal imprinting of 3.
#' simulateData(S = c(1, 2, 4), If = 3)
#'
#' ## Paternal imprinting by environment interaction of 1.5.
#' simulateData(V = c(1, 1.5, 1.5^2), propE = 0.3, Einteraction = "If")
#'
#' ## Maternal gene environment interaction of 1.6 with controls
#' simulateData(nControl = 1000, V = c(1, 1.6, 1.6^2), propE = 0.3, Einteraction = "M")
#'
#' ## Null model with 3 different sub-populations
#' simulateData(nPop = 3, maf = c(0.1, 0.2, 0.3), prev.byPop = c(0.2, 0.1, 0.4), prop.byPop = c(0.3, 0.3, 0.4))
#'
#'## Null model with 2 different sub-populations, environmental exposures, and controls
#' simulateData(
#'   nControl = 1000,
#'   propE = c(0.1, 0.4),
#'   propE.control = c(0.2, 0.2),
#'   nPop = 2,
#'   maf = c(0.3, 0.4),
#'   prev.byPop = c(0.2, 0.3),
#'   prop.byPop = c(0.6, 0.4)
#'  )
#'
#' @importFrom Rdpack reprompt
simulateData <- function(nCases = 1000, nControl = 0,
                         R = c(1, 1, 1), S = c(1, 1, 1), V = c(1, 1, 1),
                         mtCoef = c(1, 1, 1), Im = 1, If = 1,
                         propE = 0, propE.control = propE, Einteraction = "M",
                         nPop = 1, maf = 0.3, prev.byPop = NULL, prop.byPop = NULL,
                         Fst = 0.005) {
  includePopStrat = nPop > 1

  # prE.control should be a random binomial parameter, NOT a ratio
  # Start off small with tests
  if(Im != 1 && If != 1 && !identical(R, c(1,1,1))){
    warning("Maternal and paternal imprinting included with child effects,
            resulting data cannot have all parameters simultaneously modelled.")
  }

  if (includePopStrat == TRUE) {
    # Ensure that the proportion of the cases from each population is provided.
    if (is.null(prev.byPop) && is.null(prop.byPop)){
      stop("Must include prev.byPop and prop.byPop if simulating population
           stratification.\n")
    }
    else {
      if (length(prev.byPop) != nPop) {
        stop("The length of prev.byPop must be nPop. Each element is the prevalence
                of disease from each subpopulation.\n")
      }

      if (round(sum(prop.byPop), 8) != 1) {
        stop("The sum of prop.byPop must equal 1. Each element is the final proportion
                of individuals from each subpopulation.\n")
      }
    }
  } else {
    prop.byPop <- 1
    prev.byPop <- 0.5 # Prevalence doesn't matter in the non popstrat scenario
  }

  if (any(V != 1) && nPop != 1){
    warning("Population stratification simulation is invalid when true gene
              environment interactions are present.")
  }

  # Ensure we know proportion of families that are control trios
  if (nControl != 0) {
    if ((nControl < 0)) {
      stop("nControl must be non-negative.\n")
    }
  }

  # If only a single frequency of environmental variable is given,
  # make all populations have the same frequency
  if (includePopStrat == TRUE) {
    if (length(propE) == 1){
      propE <- rep(propE, nPop)
    }
    if (length(propE.control) == 1){
      propE.control <- rep(propE.control, nPop)
    }
    if (length(propE) != nPop || length(propE.control) != nPop){
      stop(paste0(length(propE), " elements in propE and ", length(propE.control),
      " elements propE.control, please enter either 1 or ", nPop, " elements for both.\n"))
    }
  }


  # Check proportion of environmental variable
  if (any(propE != 0)) {
    if (any(propE <= 0) || any(propE > 1)) {
      stop("Environmental exposure proportions for cases must be between 0 and 1.\n")
    }
    if (any(propE.control <= 0) || any(propE.control > 1)) {
      stop("Environmental exposure proportions for controls must be between 0 and 1.\n")
    }
  } else {
    if (any(propE.control != 0)){
      warning("propE.control is non-zero when propE is 0, setting propE.control
              to 0. No environmental exposures will be simulated.\n")
    }
    propE <- 0 # In case these were not 0 by accident
    propE.control <- propE
  }


  # Use Balding-Nichols for MAF in subpopulations
  # MAF is supposed to be very dependent on population
  # Fst can be obtained somewhat easily (e.g., Europe, Africa & Europe)
  # under this sim model, observed Fst approx the set Fst
  # Can put the Balding-Nichols model as an example in the doc
  if ((includePopStrat == TRUE) && (length(maf) == 1)) {
    warning(paste0("MAF of subpopulations not provided. Will use Balding-Nichols
                   model with Fst = ", Fst, ".\n"))
    alpha <- maf * (1 - Fst) / Fst
    beta <- (1 - maf) * (1 - Fst) / Fst

    # Generate allele frequencies for the nPop populations
    q <- rbeta(nPop, shape1 = alpha, shape2 = beta)
  } else {
    q <- maf
  }


  # Count tables across all populations
  caseE0.all <- NULL
  caseE1.all <- NULL
  controlE0.all <- NULL
  controlE1.all <- NULL
  full.tables <- NULL

  prControl.E.byPop = (1 - prev.byPop) * propE.control * prop.byPop / sum((1 - prev.byPop) * prop.byPop)
  prCase.E.byPop = prev.byPop * propE * prop.byPop / sum(prev.byPop * prop.byPop)
  prControl.notE.byPop = (1 - prev.byPop) * (1 - propE.control) * prop.byPop / sum((1 - prev.byPop) * prop.byPop)
  prCase.notE.byPop = prev.byPop * (1 - propE) * prop.byPop / sum(prev.byPop * prop.byPop)
  ## THESE EQUATIONS ARE VALID IFF WE MAKE THE ASSUMPTION THAT D AND E ARE CONDITIONALLY INDEPENDENT (V = 1)

  prCase.byPop = c(prCase.E.byPop, prCase.notE.byPop)
  prControl.byPop = c(prControl.E.byPop, prControl.notE.byPop)
  ## These will be 2*npop vectors, first npop is for E = 1, last npop is for E = 0

  sampled_pop_E_counts <-
    rmultinom(
      n = 1,
      size = nCases,
      prob = c(
        prCase.byPop
      )
    ) %>%
    cbind(
      rmultinom(
        n = 1,
        size = nControl,
        prob = c(
          prControl.byPop
        )
      )
    )
  ## Gives us a two column, 2 * nPop row matrix where col1 and col2 are the counts
  ## for number of cases and controls, respectively,
  ## where first nPop rows represent the subpopulations with E = 1, and last
  ## nPop rows represent subpops with E = 0.

  ## REMEMBER THAT THIS ONLY WORKS IF WE ASSUME D AND E ARE INDEPENDENT

  # Create all the datasets separately in each population
  for (i in 1:nPop) {
    ## NEW CODE
    ntrios.pop.E.case <- sampled_pop_E_counts[i, 1]
    ntrios.pop.notE.case <- sampled_pop_E_counts[i + nPop, 1]
    ntrios.pop.E.control <- sampled_pop_E_counts[i, 2]
    ntrios.pop.notE.control <- sampled_pop_E_counts[i + nPop, 2]

    ## OLD CODE
    # ntrios.pop.case <- round(ntrios * (1 - prControl) * prCase.byPop[i])
    # ntrios.pop.control <- round(ntrios * prControl * prControl.byPop[i])
    #
    # ntrios.pop.notE.case <- round(ntrios.pop.case * (1 - prE[i]))
    # ntrios.pop.E.case <- round(ntrios.pop.case * prE[i])
    # ntrios.pop.notE.control <- round(ntrios.pop.control * (1 - prE.control[i]))
    # ntrios.pop.E.control <- round(ntrios.pop.control * prE.control[i])


    # Simulate "case" trios,
    caseE0 <- simulateDataSubset(
      ntrios = ntrios.pop.notE.case, maf = q[i], R = R, S = S, mtCoef = mtCoef,
      includeE = FALSE, Im = Im, If = If
    )
    if (i == 1) {
      caseE0.all <- caseE0
      full.tables <- caseE0$datFull |> dplyr::mutate(pop = i)
    } else {
      caseE0.all <- mergeCounts(caseE0.all, caseE0)
      full.tables <- rbind(full.tables, caseE0$datFull |> dplyr::mutate(pop = i))
    }

    if (any(propE != 0)) {
      caseE1 <- simulateDataSubset(
        ntrios = ntrios.pop.E.case, maf = q[i], R = R, S = S, mtCoef = mtCoef,
        V = V, includeE = TRUE, Einteraction = Einteraction, Im = Im, If = If
      )
      if (i == 1) {
        caseE1.all <- caseE1
      } else {
        caseE1.all <- mergeCounts(caseE1.all, caseE1)
      }
      full.tables <- rbind(full.tables, caseE1$datFull |> dplyr::mutate(pop = i))
    }

    # Simulate control trios
    if (nControl > 0) {
      controlE0 <- simulateDataSubset(
        ntrios = ntrios.pop.notE.control, maf = q[i],
        R = c(1, 1, 1), S = c(1, 1, 1),
        mtCoef = mtCoef, includeE = FALSE, includeControl = TRUE
      )
      if (i == 1) {
        controlE0.all <- controlE0
      } else {
        controlE0.all <- mergeCounts(controlE0.all, controlE0)
      }
      full.tables <- rbind(full.tables, controlE0$datFull |> dplyr::mutate(pop = i))

      if (any(propE != 0)) {
        controlE1 <- simulateDataSubset(
          ntrios = ntrios.pop.E.control, maf = q[i],
          R = c(1, 1, 1), S = c(1, 1, 1), mtCoef = mtCoef,
          V = c(1, 1, 1), includeE = TRUE, Einteraction = Einteraction,
          includeControl = TRUE
        )
        if (i == 1) {
          controlE1.all <- controlE1
        } else {
          controlE1.all <- mergeCounts(controlE1.all, controlE1)
        }
        full.tables <- rbind(full.tables, controlE1$datFull |> dplyr::mutate(pop = i))
      }
    }
  }

  # Create final dataset by stacking all the individual tables (E=T/F, control=T/F)
  finaldat <- list(dat4R = caseE0.all$dat4R, dat4haplin = caseE0.all$dat4haplin)

  if (any(propE != 0)) {
    finaldat <- stackCounts(finaldat, caseE1.all)
  }
  if (nControl > 0) {
    finaldat <- stackCounts(finaldat, controlE0.all)
    if (any(propE != 0)) {
      finaldat <- stackCounts(finaldat, controlE1.all)
    }
  }

  peddat <- createPed(finaldat$dat4R)

  # return(list(
  #   dat4R = finaldat$dat4R, dat4haplin = finaldat$dat4haplin, dat4EMIM = peddat,
  #   datAll = full.tables
  # ))
  attr(finaldat$dat4R, "ALL") <- full.tables
  class(finaldat$dat4R) <- c("TriLLIEM.sim", "data.frame")
  return(finaldat$dat4R)
}
