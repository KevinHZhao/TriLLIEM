residuals.TriLLIEM <- function(object, type = c("deviance", "pearson", "working",
                                                "response", "partial"), ...)  {
  ## Modifying the code for residuals.glm
  type <- match.arg(type)
  y <- TriLLIEM:::sum_grp(object$y, object$grp)
  ## Fix working residuals formula
  r <- TriLLIEM:::sum_grp(object$residuals * object$fitted.values, object$grp) / TriLLIEM:::sum_grp(object$fitted.values, object$grp)
  mu <- TriLLIEM:::sum_grp(object$fitted.values, object$grp)
  wts <- TriLLIEM:::mean_grp(object$prior.weights, object$grp)
  grp <- object$grp
  res <- switch(type, deviance = if (object$df.residual > 0) {
    d.res <- sqrt(pmax((object$family$dev.resids)(y, mu,
                                                  wts), 0))
    ifelse(y > mu, d.res, -d.res)
  } else rep.int(0, length(mu)), pearson = (y - mu) * sqrt(wts)/sqrt(object$family$variance(mu)),
  working = r, response = y - mu, partial = r)
  if (!is.null(object$na.action))
    res <- naresid(object$na.action, res)
  if (type == "partial")
    res <- res + predict(object, type = "terms")
  res
}
