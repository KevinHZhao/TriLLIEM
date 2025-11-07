## Calculates the Fisher info matrix, accounting for EM

#' Calculate Variance-Covariance Matrix for a `TriLLIEM` Object
#'
#' @param res a `TriLLIEM` object.
#'
#' @return A matrix of covariances between parameters, accounting for the EM
#' algorithm.
#' @export
#' @keywords internal
#'
#' @examples
#' model <- TriLLIEM(dat = example_dat4R, effects = c("C", "M", "Im"))
#' vcov(model)
vcov.TriLLIEM <- function(res) {
  proven <- FALSE # need to fix to work with controls and E
  ## Calculating SE's
  ## IX is the 16-row info matrix
  ## IY is the 15-row info matrix

  n <- length(res$y)

  if(n %% 16 != 0) {
    return(NextMethod())
  }

  ## I_X <- vcov(res) is slightly less precise
  Z <- res$x
  I_X <- t(Z) %*% diag(res$fitted.values) %*% Z

  ## Matrix that transforms complete data (X) into observed data (Y)
  L <- diag(n)
  L[seq(9, n, by = 16),] <-
    L[seq(9, n, by = 16),] + L[seq(10, n, by = 16),]
  L[sapply(10:15, FUN = seq, to = n, by = 16) %>% t() %>% c(),] <-
    L[sapply(11:16, FUN = seq, to = n, by = 16) %>% t() %>% c(),]
  L <- L[-(seq(16, n, by = 16)),]

  # Try to clean up so don't need to calculate all the zeros
  # Looks like var_XgY is only non-zero in the 9 and 10 rows/cols, likely scales
  # up by 16 each time
  # MUST PROVE that covariance at 9,9 = 10,10 = -9,10 = -10,9
  if(proven){
    var_XgY <- vapply(X = 1:n,
                      FUN = function(x){
                        vapply(X = 1:n,
                               cov_trill,
                               a = x,
                               y = L %*% res$y,
                               L = L,
                               mu = fitted(res),
                               FUN.VALUE = double(1)
                        )
                      },
                      FUN.VALUE = double(n)
    )
  } else{
    var_XgY <- matrix(data = 0, nrow = n, ncol = n)
    for(i in seq(9, n, by = 16)){
      var <- cov_trill(y = L %*% res$y, L = L, a = 9, b = 9, mu = fitted(res))
      ## Can test cov_trill = var...
      var_XgY[i, i] <- var
      var_XgY[i+1, i+1] <- var
      var_XgY[i+1, i] <- -var
      var_XgY[i, i+1] <- -var
    }
  }
  I_XgY <- t(Z) %*% var_XgY %*% Z
  I_Y <- I_X - I_XgY
  return(solve(I_Y))
}
