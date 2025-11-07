#' Summary function for TriLLIEM functions
#'
#' @param res Resulting object of class "`glm` from the [TriLLIEM()] function.
#'
#' @return A list with the following components:
#' \describe{
#'  \item{effects}{The maximum likelihood estimates of the coefficients in the model.}
#'  \item{se}{The standard errors of the maximum likelihood estimates, based on the observed information matrix.}
#'  \item{pvals}{The p-values for each of the maximum likelihood estimates.}}
#' @export
#'
#' @examples
#' res <- TriLLIEM(dat = example_dat4R)
#' summ_trill(res)
summary.TriLLIEM <- function (object)
{
  df.r <- object$df.residual
  aliased <- is.na(coef(object))
  p <- object$rank
  if (p > 0) {
    p1 <- 1L:p
    Qr <- qr(object)
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat <- vcov(object)
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    pvalue <- 2 * pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
    df.f <- NCOL(covmat)
  }
  keep <- match(c("call", "terms", "family", "deviance", "aic",
                  "contrasts", "df.residual", "null.deviance", "df.null",
                  "iter", "EM_iter", "na.action"), names(object), 0L)
  ans <- c(
    object[keep],
    list(
      deviance.resid = TriLLIEM:::residuals.TriLLIEM(object, type = "deviance"),
      coefficients = coef.table,
      aliased = aliased,
      dispersion = 1,
      df = c(object$rank, df.r, df.f),
      cov.unscaled = covmat,
      cov.scaled = covmat
    )
  )
  class(ans) <- c("summary.TriLLIEM", "summary.glm")
  return(ans)
}

## summary for trill objects
## Add anova.TriLLIEM for LRT
## Document that summary.glm from base r was modified to get this...
