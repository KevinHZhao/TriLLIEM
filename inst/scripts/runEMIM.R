runEMIM <- function(mtmodel = "MS", effects = c("C", "M"), peddat,
                    emimpath = "~/EMIM/", Einteraction = "M",
                    includeIm = FALSE, includeIf = FALSE, weinberg = FALSE,
                    includeE = FALSE) {
  ## Set up temp wd so EMIM files don't show up
  withr::local_dir(new = withr::local_tempdir())

  if(includeIm){
    effects <- c(effects, "Im")
  }
  if(includeIf){
    effects <- c(effects, "If")
  }

  if (all(c("C", "Im", "If") %in% effects)){
    warning("Cannot include maternal and paternal imprinting with child effects.")
  }
  # If Einteraction is not in effect, give warning and add it to effect
  if (!(Einteraction %in% effects) && includeE == TRUE){
    warning("Einteraction is not in effects, including it manually.")
    effects <- c(effects, Einteraction)
  }

  ## This is for ensuring EMIM knows to use "2" as the risk allele (otherwise it
  ## will default to the least common allele)
  write(c("1", "A", "2"), "emim_rfile", ncol = 3)

  # Setup results objects
  neweffects <- c(
    setdiff(effects, c(paste0("E:", Einteraction), Einteraction)),
    paste0(Einteraction, "[E=0]"),
    paste0(Einteraction, "[E=1]"),
    paste0(Einteraction, "[E=1]/", Einteraction, "[E=0]")
  )


  # Set up options for the parameter file
  options <- " -a -so -rfile emim_rfile"

  if (is.element("C", effects)) { # Multiplicative allele model for C effect
    options <- c(options, "-ct")
  }
  if (is.element("Im", effects)) { # Multiplicative allele model for C effect
    options <- c(options, "-im")
  }
  if (is.element("If", effects)) { # Multiplicative allele model for C effect
    options <- c(options, "-ip")
  }

  if (is.element("M", effects)) { # Multiplicative allele model for M effect
    options <- c(options, "-mt")
  }
  options <- paste0(options, " ", collapse = " ")


  # Note that ped data includes an E column, which must
  # must be removed
  if (includeE) {
    ## Subset the data into cases where E=0 and once with E=1
    subset0 <- peddat |> dplyr::filter(E == 0)
    subset1 <- peddat |> dplyr::filter(E == 1)

    # All
    peddat <- peddat[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_all.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree_all.map", ncol = 4)

    # E=0 data
    peddat <- subset0[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_0.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree_0.map", ncol = 4)

    # E=1 data
    peddat <- subset1[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_1.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree_1.map", ncol = 4)

    ## Run PREMIM and EMIM
    for (i in c("all", 0, 1)) {
      # Run PREMIM
      command <- paste0(emimpath, "premim", options, "temp_pedigree_", i, ".ped temp_pedigree_", i, ".map")
      system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

      # Change parameter in options file for mating symmetry
      if (mtmodel == "MS") {
        params <- readLines("emimparams.dat")
        params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
        write(params, "emimparams.dat", ncol = 1)
      }
      if (mtmodel == "MaS") {
        params <- readLines("emimparams.dat")
        params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
        params[18] <- "1   << use CPG likelihood (estimate 9 mu parameters)"
        write(params, "emimparams.dat", ncol = 1)
      }

      # run EMIM and rename results
      command <- paste0(emimpath,
                        "emim",
                        ifelse(.Platform$OS.type == "unix",
                               "",
                               ".exe")
                        )
      system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
      system(paste0("mv emimsummary.out emimsummary_", i, ".out"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    }

    ## Read in results and parse
    resAll <- read.table("emimsummary_all.out", header = T)
    res0 <- read.table("emimsummary_0.out", header = T)
    res1 <- read.table("emimsummary_1.out", header = T)
    res <- mget(c("resAll", "res0", "res1"))
  } else {
    peddat <- peddat[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]

    write.table(peddat, "temp_pedigree.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree.map", ncol = 4)


    # Run PREMIM
    command <- paste0(emimpath,
                      "premim",
                      ifelse(.Platform$OS.type == "unix",
                             "",
                             ".exe"
                             ),
                      options,
                      "temp_pedigree.ped temp_pedigree.map")
    system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Change parameter for mating symmetry
    if (mtmodel == "MS") {
      params <- readLines("emimparams.dat")
      params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
      write(params, "emimparams.dat", ncol = 1)
    }
    if (mtmodel == "MaS") {
      params <- readLines("emimparams.dat")
      params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
      params[18] <- "1   << use CPG likelihood (estimate 9 mu parameters)"
      write(params, "emimparams.dat", ncol = 1)
    }

    # RUN EMIM
    command <- paste0(emimpath,
                      "emim",
                      ifelse(.Platform$OS.type == "unix",
                             "",
                             ".exe")
                      )
    system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Read in results
    res <- read.table("emimsummary.out", header = T)
  }

  attr(res, "includeE") <- includeE
  attr(res, "Einteraction") <- Einteraction

  return(res)
}

summ_emim <- function(res){
  resVec <- c()
  sdVec <- c()
  pvalVec <- c()
  if(attr(res, "includeE")) {
    list2env(res, envir = environment())
    resVec["C"] <- exp(resAll$lnR1)
    sdVec["C"] <- resAll$sd_lnR1
    resVec["M"] <- exp(resAll$lnS1)
    sdVec["M"] <- resAll$sd_lnS1
    pvalVec["C"] <- 2 * pnorm(abs(resAll$lnR1 / resAll$sd_lnR1), lower = F)
    pvalVec["M"] <- 2 * pnorm(abs(resAll$lnS1 / resAll$sd_lnS1), lower = F)
    resVec["Im"] <- exp(resAll$lnIm)
    resVec["If"] <- exp(resAll$lnIp)

    if(attr(res, "Einteraction") == "M"){
      resVec <- resVec[! names(resVec) == "M"]
      resVec["M[E=0]"] <- exp(res0$lnS1)
      resVec["M[E=1]"] <- exp(res1$lnS1)
      resVec["E:M"] <- resVec["M[E=1]"] / resVec["M[E=0]"]

      # Get a Wald-type GE test like Haplin
      z <- abs(res0$lnS1 - res1$lnS1) / sqrt(res0$sd_lnS1^2 + res1$sd_lnS1^2)
      pvalVec["E:M"] <- 2 * pnorm(z, lower = F)
    } else if(attr(res, "Einteraction") == "Im"){ ## CHeck over notes
      resVec <- resVec[! names(resVec) == "Im"]
      resVec["Im[E=0]"] <- exp(res0$lnIm)
      resVec["Im[E=1]"] <- exp(res1$lnIm)
      resVec["E:Im"] <- resVec["Im[E=1]"] / resVec["Im[E=0]"]

      # Get a Wald-type GE test like Haplin
      z <- abs(res0$lnIm - res1$lnIm) / sqrt(res0$sd_lnIm^2 + res1$sd_lnIm^2)
      pvalVec["E:Im"] <- 2 * pnorm(z, lower = F)
    } else if(attr(res, "Einteraction") == "If"){
      resVec <- resVec[! names(resVec) == "If"]
      resVec["If[E=0]"] <- exp(res0$lnIp)
      resVec["If[E=1]"] <- exp(res1$lnIp)
      resVec["E:If"] <- resVec["If[E=1]"] / resVec["If[E=0]"]

      # Get a Wald-type GE test like Haplin
      z <- abs(res0$lnIp - res1$lnIp) / sqrt(res0$sd_lnIp^2 + res1$sd_lnIp^2)
      pvalVec["E:If"] <- 2 * pnorm(z, lower = F)
    }

    sdVec["Im"] <- resAll$sd_lnIm
    pvalVec["Im"] <- 2 * pnorm(abs(resAll$lnIm / resAll$sd_lnIm), lower = F)
    sdVec["If"] <- resAll$sd_lnIp
    pvalVec["If"] <- 2 * pnorm(abs(resAll$lnIp / resAll$sd_lnIp), lower = F)
  } else {
    resVec["C"] <- exp(res$lnR1)
    sdVec["C"] <- res$sd_lnR1
    pvalVec["C"] <- 2 * pnorm(abs(res$lnR1 / res$sd_lnR1), lower = F)
    resVec["M"] <- exp(res$lnS1)
    sdVec["M"] <- res$sd_lnS1
    pvalVec["M"] <- 2 * pnorm(abs(res$lnS1 / res$sd_lnS1), lower = F)
    resVec["Im"] <- exp(res$lnIm)
    sdVec["Im"] <- res$sd_lnIm
    pvalVec["Im"] <- 2 * pnorm(abs(res$lnIm / res$sd_lnIm), lower = F)
    resVec["If"] <- exp(res$lnIp)
    sdVec["If"] <- res$sd_lnIp
    pvalVec["If"] <- 2 * pnorm(abs(res$lnIp / res$sd_lnIp), lower = F)
  }
  return(list(effects = resVec, se = sdVec, pvals = pvalVec))
}
