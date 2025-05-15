#' Title
#'
#' @param mtmodel
#' @param effects
#' @param peddat
#' @param emimpath
#' @param includeIm
#' @param includeIf
#' @param weinberg
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
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
    subset0 <- peddat %>% dplyr::filter(E == 0)
    subset1 <- peddat %>% dplyr::filter(E == 1)

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
