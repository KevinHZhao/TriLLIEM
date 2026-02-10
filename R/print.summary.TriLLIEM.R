#' Print method for `summary.TriLLIEM` objects
#'
#' @param x an object of class "`summary.TriLLIEM`", usually, a result of a call
#' to [summary.TriLLIEM].
#' @param digits the number of significant digits to use when
#' printing.
#' @param signif.stars logical. If `TRUE`, 'significance stars' are printed for
#' each coefficient, with significance codes shown under the coefficients table.
#' @param ... additional arguments passed to [printCoefmat].
#'
#' @returns Prints summary of the `TriLLIEM` model, displaying the original call
#' to the function, the matrix of coefficients, the AIC, and the number of Fisher
#' Scoring and EM iterations.
#'
#' @seealso [print.summary.glm]
#'
#' @examples
#' res <- TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im"), dat = example_dat4R)
#' print(summary(res))
#'
#' @export
print.summary.TriLLIEM <- function (x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"),
                                    ...
)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  regex_filter <-
    "^as\\.factor\\(mt\\_MS\\)[1-6](\\:E)?$|^\\(Intercept\\)$|^HWgeno(\\:E)?$|^as\\.factor\\(mt\\_MaS\\)[1-9](\\:E)?$|^E$|^D$|^E\\:D$"
  positions <- grep(regex_filter, rownames(stats::coef(x)))

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
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }

  if (nzchar(mess <- stats::naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
      "\n\n", "Number of Fisher Scoring iterations: ", x$iter,
      "\n", sep = "")
  if (x$EM_iter > 0L)
    cat("Number of EM iterations: ", x$EM_iter, "\n", sep = "")
  correl <- x$correlation
  cat("\n")
  invisible(x)
}
## Copy and pasted from getS3method("print", "summary.glm")
