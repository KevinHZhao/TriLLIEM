#' Print method for `TriLLIEM` objects
#'
#' @param x an object of class `TriLLIEM`, usually, a result of a call to
#' [TriLLIEM].
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods. (FORMATTING HERE
#' WAS WORD FOR WORD TAKEN DIRECTLY FROM THE DOCUMENTATION OF ?print.glm
#' , MAKING NOTE OF THIS IN CASE THIS NEEDS TO BE CITED SOMEHOW)
#'
#' @returns Prints details of the `TriLLIEM` model, with nuisance parameters
#' (e.g., mating type parameters) omitted.  To view all fitted parameters,
#' run `coef(x)`.
#'
#' @examples
#' res <- TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im"), dat = example_dat4R)
#' print(res)
#'
#' @export
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
# mention in details code adapted from print.glm
