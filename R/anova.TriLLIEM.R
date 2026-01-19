
#' ANOVA method for TriLLIEM objects
#'
#' This method is modelled after the `anova.glm` method.  It produces an
#' analysis of deviance table for multiple nested models.
#'
#' @param object An object of class `TriLLIEM`.
#' @param ... Additional `TriLLIEM` objects for comparison to `object`
#'
#' @details
#' Like `anova.glm`, each model's residual degrees of freedom and deviances are
#' given, alongside their respective differences between the models.
#' Models should be nested for these
#' results to be statistically interpretable. The last column shows the p-value
#' from chi-squared tests comparing the difference in deviance for each model.
#'
#' `TriLLIEM` objects modelling any sort of imprinting effect must use this
#' function, as the EM algorithm used in `TriLLIEM` causes the `anova.glm()`
#' function to treat the estimated values as having been truly observed,
#' modifying the degrees of freedom.
#'
#' @returns An object of class "`anova.TriLLIEM`" inheriting from class
#' "`anova`".
#'
#' @seealso [anova.glm()]
#' @export
#'
#' @examples
#' model1 <- TriLLIEM(dat = example_dat4R, effects = c("C"))
#' model2 <- TriLLIEM(dat = example_dat4R, effects = c("C", "M"))
#' anova(model1, model2)
anova.TriLLIEM <- function (object, ...) {
  ## All anova tests should use these parameters...
  dispersion <- NULL
  test <- "LRT"

  #stop()## Make it only work when multiple objects specified

  res <- stats:::anova.glm(object, ..., dispersion = dispersion, test = test)
  ## ANOVA will not call TriLLIEM(), instead it directly uses glm.fit, meaning
  ## df's will be incorrectly shifted for certain effects...
  df_shift <- length(object$grp)

  res[2, "Df"] <- res[2, "Df"] + df_shift
  res[2:(nrow(res) - 1), "Resid. Df"] <- res[2:(nrow(res) - 1), "Resid. Df"] - df_shift

  res[nrow(res), "Df"] <- res[nrow(res), "Df"] - df_shift

  ## Using code from anova.glm
  ## Only using LRT chisq, so limit to just that, remove F test
  fam <- object$family
  dispersion <- summary(object)$dispersion
  df.dispersion <- Inf

  ## df.dispersion = Inf for poisson family, and dispersion should always be 1
  anova_res <- stat.anova(table = res[, -ncol(res)], test = test, scale = 1, df.scale = Inf, n = TriLLIEM::nobs.TriLLIEM(object))

  ## Correcting res with the proper test
  res[, ncol(res)] <- anova_res[, ncol(anova_res)]

  coefs <- coef(object)
  regex_filter <-
    "^as\\.factor\\(mt\\_MS\\)[1-6](\\:E)?$|^\\(Intercept\\)$|^HWgeno(\\:E)?$|^as\\.factor\\(mt\\_MaS\\)[1-9](\\:E)?$|^E$|^D$|^E\\:D$"
  positions <- grep(regex_filter, names(coefs))
  pasted_factors <- paste(names(coefs)[-positions], collapse = ", ")

  class(res) <- c("anova.TriLLIEM", class(res))
  attr(res, "heading") <- paste0("Analysis of Deviance Table\n\nFactors included in model: ", pasted_factors, "\n\nResponse: count \n\nTerms added sequentially (first to last)\n\n")
  return(res)
}
