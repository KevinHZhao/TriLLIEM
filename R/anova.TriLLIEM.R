
#' Analysis of Deviance for TriLLIEM objects
#'
#' @param object
#' @param ...
#' @param dispersion
#' @param test
#'
#' @returns
#' @export
#'
#' @examples
anova.TriLLIEM <- function (object, ..., dispersion = NULL, test = "LRT") {
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
  if (is.null(test)) {
    test <- if (!is.null(dispersion))
      "Chisq"
    else if (!is.null(fam$dispersion))
      if (is.na(fam$dispersion))
          "F"
      else "Chisq"
    else FALSE
  }
  if (!isFALSE(test)) {
    if (is.null(dispersion)) {
      dispersion <- summary(object)$dispersion
    }
    df.dispersion <- if (is.null(fam$dispersion))
      if (isTRUE(dispersion == 1))
        Inf
      else object$df.residual
    else if (is.na(fam$dispersion))
      object$df.residual
    else Inf
    if (isTRUE(test == "F") && df.dispersion == 0) {
      test <- FALSE
    }
  }

  ## df.dispersion = Inf for poisson family, and dispersion should always be 1
  anova_res <- stat.anova(table = res[, -ncol(res)], test = test, scale = 1, df.scale = Inf, n = TriLLIEM::nobs.TriLLIEM(object))

  ## Correcting res with the proper test
  res[, ncol(res)] <- anova_res[, ncol(anova_res)]

  coefs <- coef(x)
  regex_filter <-
    "^as\\.factor\\(mt\\_MS\\)[1-6](\\:E)?$|^\\(Intercept\\)$|^HWgeno(\\:E)?$|^as\\.factor\\(mt\\_MaS\\)[1-9](\\:E)?$|^E$|^D$|^E\\:D$"
  positions <- grep(regex_filter, names(coefs))
  pasted_factors <- paste(names(coefs)[-positions], collapse = ", ")

  class(res) <- c("anova.TriLLIEM", class(res))
  attr(res, "heading") <- paste0("Analysis of Deviance Table\n\nFactors included in model: ", pasted_factors, "\n\nResponse: count \n\nTerms added sequentially (first to last)\n\n")
  return(res)
}
