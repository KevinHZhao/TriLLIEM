#' Simulate data for TriLLIEM
#'
#' @param ntrios Number of trios to simulate.
#' @param maf Minor allele frequency in the population, a proportion between 0 and 1.
#' @param R Vector of 3 elements representing child effects for 0, 1, and 2
#' copies of the risk allele, respectively.
#' @param S Vector of 3 elements representing maternal effects for 0, 1, and 2
#' copies of the risk allele, respectively.
#' @param V Vector of 3 elements representing gene-environment effects
#' for 0, 1, and 2 copies of the risk allele, respectively.
#' @param mtCoef Mating type coefficients.
#' @param Im Maternal imprinting effect.
#' @param If Paternal imprinting effect.
#' @param includeE A logical value indicating whether environmental effects
#' should be included in the simulation.
#' @param Einteraction A string indicating what variable environmental effects
#' interact with.  Can be "`Im`", "`If`", "`C`", or "`M`".
#' @param propE Proportion of case trios in the environmental exposure group by
#' population, assumed equal to case population proportion.
#' @param includeControl A logical value indicating whether controls should
#' be included in the simulations.
#' @param nControl Number of control trios, must be less than ntrios.
#' @param propE.control Proportion of control trios in the environmental exposure
#' group by population, assumed equal to case population proportion.
#' @param includePopStrat A logical value indicating whether to include
#' population stratification in the simulation.
#' @param numPop
#' @param Fst
#' @param prev.byPop Prevalence of cases in each sub population.
#' @param prop.byPop Proportion of each sub population, must sum to 1.
#'
#' @return A data frame of the same format as [example_dat4R]
#' @export
#'
#' @examples
#' ## Simulating data with multiplicative maternal effect of 2, and paternal imprinting of 3.
#' simulateData(S = c(1, 2, 4), If = 3)
simulateData <- function(ntrios = 1000, maf = 0.3,
                         R = c(1, 1, 1), S = c(1, 1, 1), V = c(1, 1, 1),
                         mtCoef = c(1, 1, 1), Im = 1, If = 1,
                         includeE = FALSE, Einteraction = "M", propE = 0,
                         includeControl = FALSE, nControl = 0, propE.control = propE,
                         includePopStrat = FALSE, numPop = 1, Fst = 0.005,
                         prev.byPop = NULL, prop.byPop = NULL) {
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
      if (length(prev.byPop) != numPop) {
        stop("The length of prev.byPop must be numPop. Each element is the prevalence
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

  # Ensure we know proportion of families that are control trios
  if (includeControl == TRUE) {
    if ((nControl <= 0) || (nControl >= ntrios)) {
      stop("nControl must be between 0 and ntrios exclusive.\n")
    }
  } else {
    nControl <- 0 # In case this was not 0 by accident
  }


  # If only a single frequency of environmental variable is given,
  # make all populations have the same frequency
  if (includePopStrat == TRUE) {
    if (length(propE) == 1){
      propE <- rep(propE, numPop)
    }
    if (length(propE.control) == 1){
      propE.control <- rep(propE.control, numPop)
    }
    if (length(propE) != numPop || length(propE.control) != numPop){
      stop(paste0(length(propE), " elements in propE and ", length(propE.control),
      " elements propE.control, please enter either 1 or ", numPop, " elements for both.\n"))
    }
  }


  # Check proportion of environmental variable
  if (includeE == TRUE) {
    if (any(propE <= 0) || any(propE > 1)) {
      stop("Environmental exposure proportions for cases must be between 0 and 1.\n")
    }
    if (any(propE.control <= 0) || any(propE.control > 1)) {
      stop("Environmental exposure proportions for controls must be between 0 and 1.\n")
    }
  } else {
    propE <- 0 # In case these were not 0 by accident
    propE.control <- propE
  }


  # Use Balding-Nichols for MAF in subpopulations
  if ((includePopStrat == TRUE) && (length(maf) == 1)) {
    warning("MAF of subpopulations not provided. Will use Balding-Nichols model\n")
    alpha <- maf * (1 - Fst) / Fst
    beta <- (1 - maf) * (1 - Fst) / Fst

    # Generate allele frequencies for the two populations
    q <- rbeta(numPop, shape1 = alpha, shape2 = beta)
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
  ## These will be 2*numpop vectors, first numpop is for E = 1, last numpop is for E = 0

  sampled_pop_E_counts <-
    rmultinom(
      n = 1,
      size = ntrios - nControl,
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
  ## Gives us a two column, 2 * numPop row matrix where col1 and col2 are the counts
  ## for number of cases and controls, respectively,
  ## where first numpop rows represent the subpopulations with E = 1, and last
  ## numpop rows represent subpops with E = 0.

  ## REMEMBER THAT THIS ONLY WORKS IF WE ASSUME D AND E ARE INDEPENDENT

  # Create all the datasets separately in each population
  for (i in 1:numPop) {
    ## NEW CODE
    ntrios.pop.E.case <- sampled_pop_E_counts[i, 1]
    ntrios.pop.notE.case <- sampled_pop_E_counts[i + numPop, 1]
    ntrios.pop.E.control <- sampled_pop_E_counts[i, 2]
    ntrios.pop.notE.control <- sampled_pop_E_counts[i + numPop, 2]

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
      full.tables <- list(condition = c(paste0("pop=", i), "case=1", "E=0"), table = caseE0$datFull)
    } else {
      caseE0.all <- mergeCounts(caseE0.all, caseE0)
      full.tables <- c(full.tables, list(
        condition = c(paste0("pop=", i), "case=1", "E=0"),
        table = caseE0$datFull
      ))
    }

    if (includeE == TRUE) {
      caseE1 <- simulateDataSubset(
        ntrios = ntrios.pop.E.case, maf = q[i], R = R, S = S, mtCoef = mtCoef,
        V = V, includeE = TRUE, Einteraction = Einteraction, Im = Im, If = If
      )
      if (i == 1) {
        caseE1.all <- caseE1
      } else {
        caseE1.all <- mergeCounts(caseE1.all, caseE1)
      }
      full.tables <- c(full.tables, list(
        condition = c(paste0("pop=", i), "case=1", "E=1"),
        table = caseE1$datFull
      ))
    }

    # Simulate control trios
    if (includeControl == TRUE) {
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
      full.tables <- c(full.tables, list(
        condition = c(paste0("pop=", i), "case=0", "E=0"),
        table = controlE0$datFull
      ))

      if (includeE == TRUE) {
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
        full.tables <- c(full.tables, list(
          condition = c(paste0("pop=", i), "case=0", "E=1"),
          table = controlE1$datFull
        ))
      }
    }
  }

  # Create final dataset by stacking all the individual tables (E=T/F, control=T/F)
  finaldat <- list(dat4R = caseE0.all$dat4R, dat4haplin = caseE0.all$dat4haplin)

  if (includeE == TRUE) {
    finaldat <- stackCounts(finaldat, caseE1.all)
  }
  if (includeControl == TRUE) {
    finaldat <- stackCounts(finaldat, controlE0.all)
    if (includeE == TRUE) {
      finaldat <- stackCounts(finaldat, controlE1.all)
    }
  }

  peddat <- createPed(finaldat$dat4R)

  return(list(
    dat4R = finaldat$dat4R, dat4haplin = finaldat$dat4haplin, dat4EMIM = peddat,
    datAll = full.tables
  ))
}
