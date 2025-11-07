#' Print method for `summary.TriLLIEM` objects
#'
#' See \code{\link[stats]{print.summary.glm}}.
#' @export
print.summary.TriLLIEM <- function (x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    symbolic.cor = x$symbolic.cor,
                                    signif.stars = getOption("show.signif.stars"),
                                    show.residuals = FALSE,
                                    ...
)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if (show.residuals) {
    cat("\nDeviance Residuals: \n")
    if (x$df.residual > 5) {
      x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                            na.rm = TRUE), c("Min", "1Q", "Median", "3Q",
                                                             "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  }

  regex_filter <-
    "^as\\.factor\\(mt\\_MS\\)[1-6](\\:E)?$|^\\(Intercept\\)$|^HWgeno(\\:E)?$|^as\\.factor\\(mt\\_MaS\\)[1-9](\\:E)?$|^E$|^D$|^E\\:D$"
  positions <- grep(regex_filter, rownames(coef(x)))

  if (length(x$aliased[-positions]) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df <- if ("df" %in% names(x))
      x[["df"]]
    else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients[-positions, , drop = FALSE]
    if (!is.null(aliased <- x$aliased[-positions]) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                               colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients[-positions,]
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }

  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
      "\n\n", "Number of Fisher Scoring iterations: ", x$iter,
      "\n", sep = "")
  if (x$EM_iter > 0L)
    cat("Number of EM iterations: ", x$EM_iter, "\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2L), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}
## Copy and pasted from getS3method("print", "summary.glm")
