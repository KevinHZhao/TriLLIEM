runHaplin <- function(effects = c("C", "M"), dat, includeD = FALSE,
                      includeE = FALSE, PoO = FALSE, verbose = FALSE) {
  ## Test with pooxe
  ## Set up temp wd so Haplin files don't show up
  withr::local_dir(new = withr::local_tempdir())

  # Number of columns that are not genotype. First column is E, second is Phenotype
  nvars <- 2

  # Write out the data for haplin, then read it in using haplin
  write.table(dat, "haplin_temp.dat", quote = F, row = F, col = F, sep = " ")
  write(c("chrom snp a", "1 rs1 0"), "temp.map", ncol = 1)
  dat.raw <- invisible(Haplin::genDataRead("haplin_temp.dat",
    format = "haplin", overwrite = TRUE, n.vars = nvars,
    map.file = "temp.map"
  ))

  if (includeE) {
    if (includeD) {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "cc.triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplinStrat(dat.processed,
        response = "mult", verbose = verbose,
        design = "cc.triad", ccvar = 2,
        printout = FALSE, strata = 1, reference = "ref.cat", maternal = TRUE,
        poo = PoO
      ))
    } else {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplinStrat(dat.processed,
        response = "mult", verbose = verbose,
        printout = FALSE, strata = 1, reference = "ref.cat", maternal = TRUE,
        poo = PoO
      ))
    }
  } else {
    if (is.element("M", effects)) {
      includeMaternal <- TRUE
    } else {
      includeMaternal <- FALSE
    }

    if (includeD) {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "cc.triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplin(dat.processed,
        response = "mult", verbose = verbose, design = "cc.triad",
        ccvar = 2, printout = FALSE, reference = "ref.cat", maternal = includeMaternal,
        poo = PoO
      ))
    } else {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplin(dat.processed,
        response = "mult", verbose = verbose,
        printout = FALSE, reference = "ref.cat", maternal = includeMaternal,
        poo = PoO
      ))
    }

    # res <- Haplin::haptable(res)[2, ]
  }
  # system("rm temp.map haplin_temp.dat")
  return(res)
  # return(list(effects = resVec, pvals = pvalVec))
}



summ_haplin <- function(res, PoO = FALSE, includeE = FALSE){
  resVec <- c()
  pvalVec <- c()
  if(includeE && PoO){
    GEtest <- Haplin::gxe(res)
    resVec["M[E=1]"]=Haplin::haptable(res)[6,"RRm.est."] #RR for E=1, M=1
    resVec["M[E=0]"]=Haplin::haptable(res)[4,"RRm.est."]  #RR for E=0, M=1
    resVec["E:M"]=resVec["M[E=1]"]/resVec["M[E=0]"]
    pvalVec["E:M"]=GEtest$gxe.test[3,"pval"] # pval for stratified test
    pvalVec["M"]=Haplin::haptable(res)[2,"RRm.p.value"] # pval for unstratified analysis
    resVec["Im"] <- Haplin::haptable(res)[2,"RRcm.est."]
    pvalVec["Im"] <- Haplin::haptable(res)[2,"RRcm.p.value"]
    resVec["If"] <- Haplin::haptable(res)[2,"RRcf.est."]
    pvalVec["If"] <- Haplin::haptable(res)[2,"RRcf.p.value"]
  }
  else if(includeE){
    GEtest <- Haplin::gxe(res)
    resVec["M[E=1]"]=Haplin::haptable(res)[6,"RRm.est."] #RR for E=1, M=1
    resVec["M[E=0]"]=Haplin::haptable(res)[4,"RRm.est."]  #RR for E=0, M=1
    resVec["E:M"]=resVec["M[E=1]"]/resVec["M[E=0]"]
    pvalVec["E:M"]=GEtest$gxe.test[3,"pval"] # pval for stratified test
    pvalVec["M"]=Haplin::haptable(res)[2,"RRm.p.value"] # pval for unstratified analysis
    resVec["C"]=Haplin::haptable(res)[2,"RR.est."]
    pvalVec["C"]=Haplin::haptable(res)[2,"RR.p.value"]
  }
  else if(PoO){
    res <- Haplin::haptable(res)[2,]
    resVec["M"] <- res$RRm.est
    pvalVec["M"] <- res$RRm.p.value
    resVec["Im"] <- res$RRcm.est
    pvalVec["Im"] <- res$RRcm.p.value
    resVec["If"] <- res$RRcf.est
    pvalVec["If"] <- res$RRcf.p.value
  }
  else{
    res <- Haplin::haptable(res)[2,]
    resVec["C"] <- res$RR.est.
    pvalVec["C"] <- res$RR.p.value
    resVec["M"] <- res$RRm.est.
    pvalVec["M"] <- res$RRm.p.value
  }
  return(list(effects = resVec, pvals = pvalVec))
}
