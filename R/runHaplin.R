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
