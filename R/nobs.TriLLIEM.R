
#' Accessing number of observations in TriLLIEM objects
#'
#' @param object
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
nobs.TriLLIEM <- function (object, ...) {
  nrow(object$y_initial)
}
