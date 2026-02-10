#' Summary function for TriLLIEM functions
#'
#' @param object an object of class "`TriLLIEM`", resulting from a call to
#' [TriLLIEM].
#' @param ... arguments passed to or from other methods.
#'
#' @details
#' Due to [TriLLIEM] using the EM algorithm when fitting imprinting effects,
#' calculation of residuals is modified compared to in [summary.glm], where
#' the original data (`object$y_initial`) is used instead of `object$y` to
#' ensure the correct residuals are computed.  See [TriLLIEM] for more details.
#'
#' @return A list with the same components as those returned by [summary.glm],
#' but with the addition of the following:
#' \describe{
#'  \item{terms}{the same component from `object`, always included for "`TriLLIEM`" objects.}
#'  \item{aic}{the same component from `object`, always included for "`TriLLIEM`" objects}
#'  \item{EM_iter}{the same component from `object`, always included for "`TriLLIEM`" objects}
#' }
#'
#' @export
#'
#' @examples
#' res <- TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im"), dat = example_dat4R)
#' summary(res)
#'
#' @importFrom stats vcov
summary.TriLLIEM <- function (object, ...)
{
  df.r <- object$df.residual
  aliased <- is.na(stats::coef(object))
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
    pvalue <- 2 * stats::pnorm(-abs(tvalue))
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
      deviance.resid = residuals.TriLLIEM(object, type = "deviance"),
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
