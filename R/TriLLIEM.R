#' Fit the log-linear model to trio data.
#'
#' @description
#' This function is used to fit the user-specified log-linear model to trio
#' count data.
#'
#'
#' @param mtmodel Mating type model to use in the analysis, can be "`HWE`" for
#' Hardy-Weinberg Equilibrium, "`MS`" for Mating Symmetry (default), and "`MaS`"
#' for Mating Asymmetry.
#' @param effects A vector listing the effects, as strings, to include
#' in the model. Effects can include:
#' \describe{
#'  \item{"`C`"}{Child effects.}
#'  \item{"`M`"}{Maternal effects.}
#'  \item{"`Im`"}{Maternal imprinting effects.}
#'  \item{"`If`"}{Paternal imprinting effects.}
#'  \item{"`E:C`"}{Child gene environment interactions.}
#'  \item{"`E:M`"}{Maternal gene environment interactions.}
#'  \item{"`E:Im`"}{Maternal imprinting by environment interactions.}
#'  \item{"`E:If`"}{Paternal imprinting by environment interactions.}
#' }
#' Default is `c("C", "M")`.
#' @param dat A data frame with triad data, with the formatting of
#' [example_dat4R].
#' @param includeE A logical value indicating whether to include environment
#' interaction effects. If set to "`FALSE`", any exposed counts in `dat` are
#' combined with the respective unexposed count (treating `dat` as if all
#' `E = 0`).  Default is `FALSE`.
#' @param Estrat A logical value indicating whether to use a stratified approach
#' for environmental interactions. Default is `FALSE`. See details for more
#'  information.
#' @param Eanova A logical value indicating if this is for the sake of running
#' anova to compare a model with environmental interactions to a model that does
#' not include environmental interactions.  Should be left as "`FALSE`" if not
#' for this purpose, as the degrees of freedom will be incorrect.  See details
#' for more information. Default is `FALSE`.
#' @param includeD A logical value indicating whether to use the hybrid
#' model with controls.  If set to "`FALSE`", any control trios will be removed
#' from the data set prior to analysis. Default is `FALSE`.
#' @param Minit Initial value for the proportion of triple heterozygote (`M=1`,
#' `F=1`, `C=1`) category where the '1' allele is passed from the mother to child.
#' This is used to initialize the EM algorithm. Default is 0.5.
#' @param max.iter Maximum number of iterations for the EM algorithm.  Default
#' is 30.
#' @param EM.diag A logical value indicating whether to show diagnostic messages
#' for the EM algorithm.  Default is `FALSE`.
#'
#' @details
#' Fits the specified log-linear model of \insertCite{Wein+98;textual}{TriLLIEM} to
#' `dat` using R's [glm] framework.  This includes \insertCite{WeinUmba05;textual}{TriLLIEM}'s
#' hybrid model allowing for control triads.  When imprinting effects ("`Im`" or
#' "`If`") are included, an EM algorithm is run to estimate the counts of
#' \eqn{(M, F, C) = (1, 1, 1)} triads which are maternally and paternally
#' inherited.  Normally, if this was done using the [glm] function, the computed
#' likelihoods would use these estimated counts, which is incorrect. This function,
#' alongside its methods like [summary.TriLLIEM] and [anova.TriLLIEM],
#' performs the EM algorithm while using the original observed counts to correctly
#' obtain likelihoods.
#'
#' The model formula supplied to [glm] is made up of the genetic effects supplied
#' in `effects`, alongside several nuisance parameters.  These include:
#'  - Mating type parameters, dependent on the specified `mtmodel`,
#'  - `D`, if `includeD = TRUE`,
#'  - `E:D`, if `includeD = TRUE` and `includeE = TRUE`,
#'  - Mating type parameter by `E` interactions, if `includeE = TRUE`; `E`, if
#'  `includeE = TRUE` and `mtmodel = "HWE"` (both necessary as shown in
#'  \insertCite{Shin+10;textual}{TriLLIEM}).
#'  These nuisance parameters are omitted when printing the function output,
#'  but may be viewed by using the [coef] function on the output.
#'
#' `Estrat = TRUE` forces every every listed effect in `effects` to have its
#' gene environment interaction included in the model, regardless of whether the
#' user has specified them explicitly or not. This essentially stratifies the
#' model by `E`, and is useful for replicating the stratified models necessary
#' for analyzing gene environment interactions in `Haplin`
#' \insertCite{GjesLie06}{TriLLIEM} and `EMIM` \insertCite{HoweCord12}{TriLLIEM}.
#'
#' `Eanova = TRUE` allows the model to be fit when \eqn{E \neq 0} rows are
#' present but `includeE = FALSE`, without modifying these exposed rows.
#' This is only useful for using `anova.TriLLIEM`
#' to determine if models with gene-environment interactions are statistically
#' different from models without gene-environment interactions.  Since this option
#' keeps these \eqn{E \neq 0} categories without using the `E` parameter though,
#' the fitted model will use incorrect degrees of freedom in its significance
#' tests, and hence should not be used for any inferences besides this very
#' specific case of model comparison.
#'
#' All `TriLLIEM` objects are of family `poisson_em`, a modified version of the
#' poisson family to account for additional data created by the EM algorithm when
#' computing residuals and AIC.
#'
#' @return An object of class "`TriLLIEM`", which inherits from class "`glm`"
#' and has the same components as the output of [glm], with the following
#' modifications:
#' \describe{
#'  \item{`df.residual`}{if imprinting effects are specified, subtracted by the
#'  number of additional rows introduced by the EM algorithm.}
#'  \item{`df.null`}{if imprinting effects are specified, subtracted by the
#'  number of additional rows introduced by the EM algorithm.}
#'  \item{`EM_iter`}{number of EM algorithm iterations before convergence.}
#'  \item{`y_initial`}{initial supplied data frame, `dat`, before addition of rows
#'  by EM.}
#'  \item{`grp`}{list of vectors, with each vector containing indices of grouped
#'  rows in the data after the EM algorithm is applied (`y`). Aggregating these grouped
#'  rows by sum will yield `y_initial`.}
#' }
#'
#' @export
#'
#' @examples
#' res1 <- TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im"), dat = example_dat4R)
#' summary(res1)
#'
#' dat <-
#'   simulateData(
#'     nControl = 1000,
#'     propE = c(0.1, 0.4),
#'     propE.control = c(0.2, 0.2),
#'     nPop = 2,
#'     maf = c(0.3, 0.4),
#'     prev.byPop = c(0.2, 0.3),
#'     prop.byPop = c(0.6, 0.4)
#'   )
#' ## Obtain the non-stratified and stratified models and compare them via anova
#' res2 <-
#'   TriLLIEM(
#'     mtmodel = "HWE",
#'     effects = c("C", "M", "Im", "E:Im"),
#'     dat = dat,
#'     includeE = TRUE,
#'     includeD = TRUE
#'   )
#' res3 <-
#'   TriLLIEM(
#'     mtmodel = "HWE",
#'     effects = c("C", "M", "Im", "E:Im"),
#'     dat = dat,
#'     includeE = TRUE,
#'     Estrat = TRUE,
#'     includeD = TRUE
#'   )
#' anova(res2, res3)
#'
#' ## Compare non-stratified model to a model without E by setting Eanova = TRUE
#' res4 <-
#'   TriLLIEM(
#'     mtmodel = "HWE",
#'     effects = c("C", "M", "Im"),
#'     dat = dat,
#'     Eanova = TRUE,
#'     includeD = TRUE
#'   )
#'
#' anova(res2, res4)
#'
#' @references{
#' \insertAllCited{}
#' }
#' @importFrom rlang .data
TriLLIEM <- function(mtmodel = "MS", effects = c("C", "M"), dat,
                     includeE = FALSE, Estrat = FALSE, Eanova = FALSE,
                     includeD = FALSE, Minit = 0.5, max.iter = 30, EM.diag = FALSE) {
  for (effect in effects){
    if (grepl("E:", effect)){
      noeeff <- sub("E:", "", effect)
      if (!(noeeff %in% effects)){
        warning(paste0(effect, " is given as an environmental interaction but ",
                       noeeff, " is not listed as an effect. Model inferences
                       may be invalid.\n"))
      }
    }
    else if (grepl(":E", effect)){
      noeeff <- sub(":E", "", effect)
      if (!(noeeff %in% effects)){
        warning(paste0(effect, " is given as an environmental interaction but ",
                       noeeff, " is not listed as an effect. Model inferences
                       may be invalid.\n"))
      }
    }
  }

  if (all(c("C", "Im", "If") %in% effects)){
    stop("Cannot include maternal and paternal imprinting with child effects.\n")
  }

  if(length(unique(dat$E)) < 2 && includeE){
    stop("E column must have at least 2 distinct values.\n")
  }
  if(length(unique(dat$D)) < 2 && includeD){
    stop("D column must have at least 2 distinct values.\n")
  }

  cal <- match.call()

  # Portion of model equation depends on mating type model
  dat <- dat |> dplyr::mutate(offset = dplyr::case_when(.data$type == 9 ~ 2, .default = 1))
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
  hasE <- FALSE
  if(includeE){
    ## If we want the same results as the stratified approach of emim and haplin, MUST include E:everything,
    ## otherwise, don't include by default to save on power
    if(Estrat){
      Eeffects <-
        c(mteffect, effects) |>
        paste0(":E", sep = "")
      if(mtmodel == "HWE")
        Eeffects <- append(Eeffects, "E")
    } else {
      Eeffects <-
        c(mteffect) |>
        paste0(":E", sep = "")
      if(mtmodel == "HWE")
        Eeffects <- append(Eeffects, "E")
    }
  } else if (!Eanova) {
    # If includeE is FALSE and not using Eanova, treat all counts as unexposed
    dat <- dat |>
      dplyr::summarize(count = base::sum(.data$count), .by = c(-"E", -"count")) |>
      dplyr::mutate(E = 0)
  } else {
    # If includeE is FALSE but we're still using Eanova
    if (any(dat$E != 0)) {
      hasE <- TRUE
    }
  }

  # Hybrid model
  Deffects <- c()
  if(includeD){
    dat <-
      dat |>
      dplyr::mutate(C = .data$C * .data$D,
                    M = .data$M * .data$D)
    Deffects <- c("D", if (includeE) "E:D")
  } else if (base::sum(dat$D == 0) != 0) {
    base::warning("Control trios detected but includeD set to FALSE.  Ignoring all control trios...\n")
    dat <- dat |> dplyr::filter(.data$D != 0)
  }

  # Portion of model equation and offset depends on mating type model
  origDat <- add_PoO_data(dat, Mprop = c(Minit, if (includeE || hasE) Minit), includeE = (includeE || hasE))

  modeleffects <- c(mteffect, Eeffects, Deffects, effects)

  linpred <- paste(modeleffects,
                   collapse = "+")
  modelformula <- paste0(modelformula, linpred)

  ## groupings of categories for em
  grp <- list()

  ## Groupings for EM
  grp <- c(grp, list(c(9,10)))
  if (includeE || hasE || includeD) {
    grp <- c(grp, list(c(25,26)))
  }
  if ((includeE || hasE) && includeD) {
    grp <- c(grp, list(c(41,42)))
    grp <- c(grp, list(c(57,58)))
  }

  # Run model and save results

  # EM for Imprinting, using same stopping criteria as Haplin...
  counter <- 0

  if(any(c("Im", "If") %in% effects)){
    if(EM.diag){
      message(paste0("Initial proportion for maternal inheritance cell = ", Minit))
    }
    prev_coeffs <- -1
    prev_dev <- -1
    repeat{
      counter <- counter + 1
      res <- suppressWarnings(
        stats::glm(stats::as.formula(modelformula), data = origDat, family = poisson_em(grp = grp), x = TRUE)
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
                       "\nDifference in deviance: ", abs(stats::deviance(res) - prev_dev),
                       "\nMax difference in coefficients: ", max(abs(stats::coef(res) - prev_coeffs))))
      }

      origDat <- add_PoO_data(dat,
                              Mprop = c(Imhat/(Ifhat + Imhat),
                                        {if (includeE) ImEhat*Imhat/(ImEhat*Imhat + IfEhat*Ifhat) else if (hasE) Imhat/(Ifhat + Imhat)}),
                              includeE = (includeE || hasE)
                              )
      ## Use a proper deviance function for imprinting
      ## Check Haplin LogLik code
      ## Make conv criteria parameters
      if(abs(stats::deviance(res) - prev_dev) < 2e-006 &&
         max(abs(stats::coef(res) - prev_coeffs)) < 1e-006){
        break
      }
      else if(counter == max.iter){
        stop("Max iterations reached without convergence.")
        break
      }
      prev_coeffs <- stats::coef(res)
      prev_dev <- stats::deviance(res)
    }

    # filled_inds <- which(dat$type == 9)
    # res$aic <- -2 * (sum(dpois(x = res$y[-c(filled_inds, filled_inds+1)], lambda = res$fitted.values[-c(filled_inds, filled_inds+1)], log = TRUE),
    #                      dpois(x = dat$count[filled_inds], lambda = (Imhat + Ifhat) * exp(res$coefficients[["as.factor(mt_MS)4"]] + sum(res$coefficients[modeleffects |> setdiff(c("Im", "If", mteffect))])), log = TRUE)
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
    ## Only include offset if we don't split the 1,1,1 cell
    ## I'm going to keep using the split 1,1,1 cell always so that anova LRT's work.
    ## EM doesn't happen in this case
    #res <- glm(as.formula(modelformula), data = dat, offset=log(dat$offset), family = poisson_em(grp = grp))
    res <- stats::glm(stats::as.formula(modelformula), data = origDat, family = poisson_em(grp = grp), x = TRUE)
    class(res) <- c("TriLLIEM", "glm", "lm")
  }

  res$EM_iter <- counter

  ## Accounting for the fact that some observations were not observed...
  res$y_initial <- dat |> dplyr::select("type", "mt_MS", "mt_MaS", "M", "F", "C", "D", "E", "count", "offset")
  res$grp <- grp
  res$df.null <- res$df.null - length(grp)
  res$df.residual <- res$df.residual - length(grp)

  res$call <- cal
  return(res)
}
