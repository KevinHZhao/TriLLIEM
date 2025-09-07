poisson_em <- function(link = "log", grp) {
  ## grp = list(c(a,b,c), (d,e), ...), represents groupings
  ## of elements to be added together to get observed data
  fam <- poisson(link = link)
  fam$family = "poisson_em"
  fam$grp = grp
  fam$dev.resids <- function(y, mu, wt) {
    if (length(mu) == 1){
      ## I think in ANOVA, last step just feeds in null model (all mu's equal),
      ## but they just recycle one mu value instead of using a vector of proper length
      mu <- rep(mu, length(y)) # This should deal with the last step in ANOVA, since incorrect behaviour if length is not same as y...
    }
    y <- TriLLIEM:::sum_grp(y, fam$grp)
    mu <- TriLLIEM:::sum_grp(mu, fam$grp)
    wt <- TriLLIEM:::mean_grp(wt, fam$grp) # This probably shouldn't matter since wt always is 1...

    r <- mu * wt
    p <- which(y > 0)
    r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
    2 * r
  }
  fam$aic <- function(y, n, mu, wt, dev) {
    y <- TriLLIEM:::sum_grp(y, fam$grp)
    mu <- TriLLIEM:::sum_grp(mu, fam$grp)
    wt <- TriLLIEM:::mean_grp(wt, fam$grp) # This probably shouldn't matter since wt always is 1...

    -2 * sum(dpois(y, mu, log = TRUE) * wt)
  }
  return(fam)
}
