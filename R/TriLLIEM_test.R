## Run a loglinear analysis with the specified model
## mtmodel = "HWE", "MS", "MaS"
## effects = "C", "M", "E:M", ... could examine environment child effects
## For EM start with 50/50 M/F, could also have sensitivity analysis
## Minit from 0 to 1, initial "guess" for proportion of maternal origin in 1,1,1
## MatImp can be true/false based on if we assume maternal/paternal imprinting
## respectively (for now assuming can't be both)
## EM.diag = show EM diagnostics
TriLLIEM_test <- function(mtmodel = "MS", effects = c("C", "M"), dat, PStest = FALSE,
                     includeIm = FALSE, includeIf = FALSE, Minit = 0.5, max.iter = 12,
                     EM.diag = FALSE) {
  if(includeIm){
    effects <- c(effects, "Im")
  }
  if(includeIf){
    effects <- c(effects, "If")
  }

  if (all(c("C", "Im", "If") %in% effects)){
    stop("Cannot include maternal and paternal imprinting with child effects.")
  }

  # Portion of model equation and offset depends on mating type model
  heteroInds <- with(dat, which((M == 1) & (F == 1) & (C == 1)))
  origDat <- add_PoO_data_15(dat)

  # Portion of model equation depends on mating type model
  # No offset as we split the (1,1,1) case
  if (mtmodel == "HWE") {
    origDat$HWgeno <- origDat$M + origDat$F
    mteffect <- "HWgeno"
    modelformula <- "count~" # Must include intercept for HW model because of log(1-p) term
  } else if (mtmodel == "MS") {
    mteffect <- "as.factor(mt_MS)"
    modelformula <- "count~-1+" # I think I can remove the intercept for MS model
  } else if (mtmodel == "MaS") {
    if (length(unique(origDat$D)) == 1) {
      stop("Only 1 phenotype in the phenotype column. Mating asymmetry models require \n
            both cases and controls\n")
    }
    mteffect <- "as.factor(mt_MaS)"
    modelformula <- "count~-1+" # I think I can remove intercept
  }

  if (is.element("E:M", effects)) { # check gene-environment paper on adding other effects ex E:C
    if (sum(origDat$D == 0) > 0) { # There are controls; include environmental interaction
      origDat$C <- origDat$C * origDat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
      origDat$M <- origDat$M * origDat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
      # (Note that the E:M term in model will be from crossing this
      # variable with E=1, which is exactly what is needed.
      if (mtmodel == "HWE") { # Include main effect of E to give different intercept for HW+E case
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects, "E", "D", "E:D")
      } else {
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects, "D", "E:D")
      }
    } else { # No controls

      if (mtmodel == "HWE") { # Include main effect of E to give different intercept for HW+E case
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects, "E")
      } else {
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects)
      }
    }
  } else { # No E:M effect

    if (sum(origDat$D == 0) > 0) {
      origDat$C <- origDat$C * origDat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
      origDat$M <- origDat$M * origDat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 O
      modeleffects <- c(mteffect, effects, "D")
    } else {
      modeleffects <- c(mteffect, effects)
    }
  }

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
  if (PStest == TRUE) {
    resVecPS <- vector(length = length(effects))
    names(resVecPS) <- effects
    pvalVecPS <- vector(length = length(effects))
    names(pvalVecPS) <- effects

    if (sum(origDat$D == 0) == 0) {
      stop("Can only test for population stratification if there are control trios\n")
    } else {
      PSeffect <- paste0(mteffect, ":D")
      modelformula.PS <- paste(modelformula, PSeffect, sep = "+")
    }
  }


  # Run model and save results

  # EM for Imprinting
  if(any(c("Im", "If") %in% effects)){
    counter <- 0
    if(EM.diag){
      message(paste0("Initial proportion for maternal inheritance cell = ", Minit))
    }
    if(all(c("Im", "If") %in% effects)){
      prev_Imhat <- 10
      prev_Ifhat <- 10
      repeat{
        counter <- counter + 1
        res <- glm(as.formula(modelformula),
                   data = origDat,
                   family = poisson(),
                   offset = c(rep(0, 8),
                              log((prev_Imhat^2 + prev_Ifhat^2)/(prev_Imhat + prev_Ifhat)),
                              rep(0,6))
                   )

        Imhat <- exp(res$coefficients["Im"])
        Ifhat <- exp(res$coefficients["If"])

        if(EM.diag){
          message(paste0("\nIteration ", counter, " of EM algorithm.
                          \nIm hat = ", Imhat,
                         "\nIf hat = ", Ifhat,
                         "\nProportion for maternal inheritance cell = ", Imhat/(Ifhat + Imhat)))
        }

        origDat <- add_PoO_data_15(dat) %>%
          dplyr::mutate(HWgeno = M + F)
        if(Imhat == prev_Imhat && Ifhat == prev_Ifhat){
          message("EM convergence achieved.")
          break
        }
        else if(counter == max.iter){
          stop("Max iterations reached without convergence.")
          break
        }
        prev_Imhat <- Imhat
        prev_Ifhat <- Ifhat
      }
    }
    else{
      prev_Ihat <- -1
      repeat{
        counter <- counter + 1
        res <- glm(as.formula(modelformula), data = origDat, family = poisson())

        Ihat <- ifelse("Im" %in% effects, exp(res$coefficients["Im"]), exp(res$coefficients["Ip"]))

        if(EM.diag){
          message(paste0("\nIteration ", counter, " of EM algorithm.
                       \nI", ifelse("Im" %in% effects,"m","f"), " hat = ", Ihat,
                         "\nProportion for maternal inheritance cell = ", ifelse("Im" %in% effects,
                                                                                 Ihat/(1 + Ihat),
                                                                                 1/(1 + Ihat))))
        }

        origDat <- add_PoO_data_15(dat) %>%
          dplyr::mutate(HWgeno = M + F)
        if(Ihat == prev_Ihat){
          message("EM convergence achieved.")
          break
        }
        else if(counter == max.iter){
          warning("Max iterations reached without convergence.")
          break
        }
        prev_Ihat <- Ihat
      }
    }
  } else{
    res <- glm(as.formula(modelformula), data = origDat, family = poisson())
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
