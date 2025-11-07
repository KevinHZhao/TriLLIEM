
#' Accessing number of observations in TriLLIEM objects
#'
#' @param object a `TriLLIEM` object
#'
#' @returns an integer representing the number of observations in the data used
#' to get `object`.
#' @export
#'
#' @examples
#' model <- TriLLIEM(dat = example_dat4R)
#' nobs(model)
nobs.TriLLIEM <- function (object) {
  nrow(object$y_initial)
}
