#' Running a loglinear analysis of trio data.
#'
#' @param mtmodel Mating type model to use in the analysis, can be "`HWE`" for
#' Hardy-Weinberg Equilibrium, "`MS`" for Mating Symmetry, and "`MaS`" for Mating Asymmetry.
#' @param effects A vector listing the effects, as strings, to include
#' in the model.  Example effects include:
#' \describe{
#'  \item{"`C`"}{Child effects.}
#'  \item{"`M`"}{Maternal effects.}
#'  \item{"`Im`"}{Maternal imprinting effects.}
#'  \item{"`If`"}{Paternal imprinting effects.}
#' }
#' @param includeE A logical value indicating whether to include environment
#' interaction effects.  If set to "`FALSE`", any exposed counts in `dat` are
#' combined with the respective unexposed count.
#' @param Einteraction A string indicating what variable environmental effects
#' interact with.  Can be "`Im`", "`If`", "`C`", or "`M`".
#' @param Estrat A logical value indicating whether to use a stratified approach
#' for environmental interactions equivalent to that of EMIM and/or Haplin.
#' @param includeD A logical value indicating whether to use the hybrid
#' model with controls.  If set to "`FALSE`", any control trios will be removed
#' from the data set prior to analysis.
#' @param dat A data frame with triad data, with the formatting of
#' [example_dat4R].
#' @param PStest A logical value indicating whether to perform a population
#' stratification test on the data.
#' @param includeIm A logical value indicating whether to include maternal
#' imprinting in the model, equivalent to adding "`Im`" in the "`effects`" vector.
#' @param includeIf A logical value indicating whether to include paternal
#' imprinting in the model, equivalent to adding "`If`" in the "`effects`" vector.
#' @param Minit Initial proportion of maternal inheritence to split the triple
#' heterozygote cell by if the EM algorithm is necessary.
#' @param max.iter Maximum number of iterations for the EM algorithm.
#' @param EM.diag A logical value indicating whether to show diagnostic messages
#' for the EM algorithm.
#'
#' @return An object of class "`glm`".
#' @export
#'
#' @examples
#' res <- TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im"), dat = example_dat4R)
#' res %>% summary() %>% coef()
TriLLIEM <- function(mtmodel = "MS", effects = c("C", "M"), dat, PStest = FALSE,
                     includeE = FALSE, Einteraction = "M", Estrat = FALSE,
                     includeD = FALSE, includeIm = FALSE,
                     includeIf = FALSE, Minit = 0.5, max.iter = 30, EM.diag = FALSE) {
  ## Test for violation of HWE when pop strat (sim MS data and compare HWE vs MS)
  ## Only one of E:M or E:Im
  ## test out with other code
  ## p values should be smaller than haplin
  ## Test with D compared to EMIM and haplin
  if (includeIm){
    effects <- c(effects, "Im")
  }
  if (includeIf){
    effects <- c(effects, "If")
  }
  # If Einteraction is not in effect, give warning and add it to effect
  if (!(Einteraction %in% effects) && includeE == TRUE){
    warning("Einteraction is not in effects, including it manually.")
    effects <- c(effects, Einteraction)
  }

  if (all(c("C", "Im", "If") %in% effects)){
    stop("Cannot include maternal and paternal imprinting with child effects.")
  }

  if(length(unique(dat$E)) < 2 && includeE){
    stop("E column must have at least 2 distinct values.")
  }
  if(length(unique(dat$D)) < 2 && includeD){
    stop("D column must have at least 2 distinct values.")
  }

  # Portion of model equation depends on mating type model
  dat <- dat %>% dplyr::mutate(offset = dplyr::case_when(type == 9 ~ 2, .default = 1))
  if (mtmodel == "HWE") {
    dat$HWgeno <- dat$M + dat$F
    mteffect <- "HWgeno"
    modelformula <- "count~" # Must include intercept for HW model because of log(1-p) term
  } else if (mtmodel == "MS") {
    mteffect <- "as.factor(mt_MS)"
    modelformula <- "count~-1+" # I think I can remove the intercept for MS model
  } else if (mtmodel == "MaS") {
    if (length(unique(dat$D)) == 1) {
      stop("Only 1 phenotype in the phenotype column. Mating asymmetry models require
            both cases and controls, and setting includeD to TRUE\n")
    }
    mteffect <- "as.factor(mt_MaS)"
    modelformula <- "count~-1+" # I think I can remove intercept
  }

  # Environmental effects
  Eeffects <- c()
  if(includeE){
    ## If we want the same results as the stratified approach of emim and haplin, MUST include E:everything,
    ## otherwise, don't include by default to save on power
    if(Estrat){
      Eeffects <-
        c(mteffect, effects) %>%
        paste0(":E", sep = "") %>%
        {if(mtmodel == "HWE") append(., "E") else .}
    } else {
      Eeffects <-
        c(mteffect, Einteraction) %>%
        paste0(":E", sep = "") %>%
        {if(mtmodel == "HWE") append(., "E") else .}
    }
  } else {
    # If includeE is FALSE, treat all counts as unexposed
    dat <- dat %>%
      summarize(count = sum(count), .by = c(-E, -count)) %>%
      dplyr::mutate(E = 0)
  }

  # Hybrid model
  Deffects <- c()
  if(includeD){
    dat <-
      dat %>%
      dplyr::mutate(C = C * D,
                    M = M * D)
    Deffects <- c("D", if (includeE) "E:D")
  } else if (sum(dat$D == 0) != 0) {
    base::warning("Control trios detected but includeD set to FALSE.  Ignoring all control trios...\n")
    dat <- dat %>% dplyr::filter(D != 0)
  }

  # Portion of model equation and offset depends on mating type model
  origDat <- add_PoO_data(dat, Mprop = c(Minit, if (includeE) Minit), includeE = includeE)

  modeleffects <- c(mteffect, Eeffects, Deffects, effects)

  linpred <- paste(modeleffects,
                   collapse = "+")
  modelformula <- paste0(modelformula, linpred)

  # Setup results objects
  resVec <- vector(length = length(effects))
  names(resVec) <- effects
  pvalVec <- vector(length = length(effects))
  names(pvalVec) <- effects
  resVecPS <- NULL
  pvalVecPS <- NULL

  # Include test and results under population stratification
  ## try catch for PStest in case glm fails due to low counts
  if (PStest == TRUE) {
    resVecPS <- vector(length = length(effects))
    names(resVecPS) <- effects
    pvalVecPS <- vector(length = length(effects))
    names(pvalVecPS) <- effects

    if (!includeD) {
      stop("Can only test for population stratification if there are control trios\n")
    } else {
      PSeffect <- paste0(mteffect, ":D")
      modelformula.PS <- paste(modelformula, PSeffect, sep = "+")
    }
  }

  # Run model and save results

  # EM for Imprinting, using same stopping criteria as Haplin...
  if(any(c("Im", "If") %in% effects)){
    counter <- 0
    if(EM.diag){
      message(paste0("Initial proportion for maternal inheritance cell = ", Minit))
    }
    prev_coeffs <- -1
    prev_dev <- -1
    repeat{
      counter <- counter + 1
      res <- suppressWarnings(
        glm(as.formula(modelformula), data = origDat, family = poisson(), x = TRUE)
      )
      class(res) <- c("TriLLIEM", "glm", "lm")

      Imhat <- ifelse(is.na(exp(res$coefficients["Im"])),
                      1,
                      exp(res$coefficients["Im"])
                      )
      Ifhat <- ifelse(is.na(exp(res$coefficients["If"])),
                      1,
                      exp(res$coefficients["If"])
                      )
      if (includeE){
        ImEhat <- ifelse(is.na(exp(res$coefficients["E:Im"])),
                        1,
                        exp(res$coefficients["E:Im"])
        )
        IfEhat <- ifelse(is.na(exp(res$coefficients["E:If"])),
                        1,
                        exp(res$coefficients["E:If"])
        )
      }

      if(EM.diag){
        message(paste0("\nIteration ", counter, " of EM algorithm.
                        \nIm hat = ", Imhat,
                       "\nIf hat = ", Ifhat,
                       if (includeE) paste0("\nE:Im hat = ", ImEhat),
                       if (includeE) paste0("\nE:If hat = ", IfEhat),
                       "\nProportion for maternal inheritance cell = ", Imhat/(Ifhat + Imhat),
                       if (includeE) paste0("\nProportion for maternal inheritance cell = ", ImEhat*Imhat/(ImEhat*Imhat + IfEhat*Ifhat)),
                       "\nDifference in deviance: ", abs(deviance(res) - prev_dev),
                       "\nMax difference in coefficients: ", max(abs(coef(res) - prev_coeffs))))
      }

      origDat <- add_PoO_data(dat,
                              Mprop = c(Imhat/(Ifhat + Imhat),
                                        if (includeE) ImEhat*Imhat/(ImEhat*Imhat + IfEhat*Ifhat)),
                              includeE = includeE
                              )
      ## Use a proper deviance function for imprinting
      ## Check Haplin LogLik code
      ## Make conv criteria parameters
      if(abs(deviance(res) - prev_dev) < 2e-006 &&
         max(abs(coef(res) - prev_coeffs)) < 1e-006){
        break
      }
      else if(counter == max.iter){
        stop("Max iterations reached without convergence.")
        break
      }
      prev_coeffs <- coef(res)
      prev_dev <- deviance(res)
    }

    # filled_inds <- which(dat$type == 9)
    # res$aic <- -2 * (sum(dpois(x = res$y[-c(filled_inds, filled_inds+1)], lambda = res$fitted.values[-c(filled_inds, filled_inds+1)], log = TRUE),
    #                      dpois(x = dat$count[filled_inds], lambda = (Imhat + Ifhat) * exp(res$coefficients[["as.factor(mt_MS)4"]] + sum(res$coefficients[modeleffects %>% setdiff(c("Im", "If", mteffect))])), log = TRUE)
    #                      )
    #                  ) +
    #   2 * res$rank
    ## saturated model set lambda = counts (perfect)
    ## Instead of using subclass, return a Trilliem object only that has the glm returned as a part of it (with incorrect), and the relevant correct info, avoids the redundant functions that return nonsense
    ## Add EM coefficients
    ## EM works with 16-rows due to the "proof"
    ## Test out emax.glm
    ## Compare EM with emax.glm performance
    ## Data abstraction and problem solving with cpp Chpt 1
  } else {
    ## Only include offset in this case when we don't split the 1,1,1 cell
    res <- glm(as.formula(modelformula), data = dat, offset=log(dat$offset), family = poisson())
    class(res) <- c("TriLLIEM", "glm", "lm")
  }

  # R is not consistent about how interaction is specified. Even though it
  # is fit as E:M, sometimes R flips it to M:E in the output of results.
  # if (is.element("M:E", names(coef(res)))) {
  #   effects[effects == "E:M"] <- "M:E"
  # }
  #
  # for (j in 1:length(effects)) {
  #   resVec[j] <- exp(summary(res)$coef[effects[j], 1])
  #   pvalVec[j] <- summary(res)$coef[effects[j], 4]
  # }
  #
  # test.res <- NULL
  # if (PStest == TRUE) {
  #   res.PS <- glm(as.formula(modelformula.PS), data = origDat, family = poisson())
  #   test.res <- anova(res, res.PS, test = "LRT")
  #
  #   for (j in 1:length(effects)) {
  #     resVecPS[j] <- exp(summary(res.PS)$coef[effects[j], 1])
  #     pvalVecPS[j] <- summary(res.PS)$coef[effects[j], 4]
  #   }
  # }

  return(res)
}
