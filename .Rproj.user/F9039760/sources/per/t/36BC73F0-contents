getHWE <- function(dat) {
  mom <- tapply(dat$count, dat$M, sum)
  dad <- tapply(dat$count, dat$F, sum)
  both <- dad + mom
  res <- c(
    Mom.exact = HWExact(mom, verbose = F)$pval,
    Dad.exact = HWExact(dad, verbose = F)$pval,
    Both.exact = HWExact(both, verbose = F)$pval
  )
  return(res)
}
