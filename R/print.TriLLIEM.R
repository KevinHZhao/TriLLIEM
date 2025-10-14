#' Title
#'
#' @param x
#' @param digits
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
#'
print.TriLLIEM <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  coefs <- coef(x)
  regex_filter <-
    "^as\\.factor\\(mt\\_MS\\)[1-6](\\:E)?$|^\\(Intercept\\)$|^HWgeno(\\:E)?$|^as\\.factor\\(mt\\_MaS\\)[1-9](\\:E)?$|^E$|^D$|^E\\:D$"
  positions <- grep(regex_filter, names(coefs))
  coefs <- coefs[-positions]

  if (length(coefs)) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ", apply(cbind(names(co), co),
                                  1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(coefs, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
      x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Null Deviance:\t   ", format(signif(x$null.deviance,
                                           digits)), "\nResidual Deviance:", format(signif(x$deviance,
                                                                                           digits)), "\tAIC:", format(signif(x$aic, digits)))
  cat("\n")
  invisible(x)
}
# code copied from getS3method("print", "glm")
