add_PoO_data <- function(dat, Mprop, includeE) {
  # Portion of model equation and offset depends on mating type model
  if (includeE && length(Mprop) != 2){
    stop("Environmental interactions included hence Mprop should be a vector of length two!")
  }
  heteroInds.E0 <- with(dat, which(type == 9 & E == 0))
  M.count.E0 <- dat |>
    dplyr::filter(.data$type == 9, .data$E == 0) |>
    dplyr::pull(.data$count) |>
    (\(x) x * Mprop[[1L]])()
  if (includeE){
    heteroInds.E1 <- with(dat, which(type == 9 & E == 1))
    M.count.E1 <- dat |>
      dplyr::filter(.data$type == 9, .data$E == 1) |>
      dplyr::pull(.data$count) |>
      (\(x) x * Mprop[[2L]])()

    PoO_dat <- dat |>
      dplyr::left_join(PoO_df |> dplyr::select("type", "matOrg", "patOrg"), by = "type") |>
      dplyr::mutate(count = base::replace(.data$count, is.na(.data$patOrg) & .data$E == 0, M.count.E0),
                    count = base::replace(.data$count, is.na(.data$patOrg) & .data$E == 1, M.count.E1),
                    matOrg = base::replace(.data$matOrg, is.na(.data$matOrg), 1),
                    patOrg = base::replace(.data$patOrg, is.na(.data$patOrg), 0)) |>
      dplyr::add_row(dat |>
                       dplyr::filter(dplyr::row_number() %in% heteroInds.E0) |>
                       dplyr::mutate(count = dat$count[heteroInds.E0] - M.count.E0,
                                     matOrg = 0,
                                     patOrg = 1)) |>
      dplyr::add_row(dat |>
                       dplyr::filter(dplyr::row_number() %in% heteroInds.E1) |>
                       dplyr::mutate(count = dat$count[heteroInds.E1] - M.count.E1,
                                     matOrg = 0,
                                     patOrg = 1)) |>
      dplyr::arrange(dplyr::desc(.data$D), dplyr::desc(.data$E), .data$type, dplyr::desc(.data$matOrg)) |>
      dplyr::mutate(typeOrig = base::rep(1:16, dplyr::n()/16),
                    Im = .data$matOrg * .data$D,
                    If = .data$patOrg * .data$D) |>
      dplyr::relocate("typeOrig")
  } else {
    PoO_dat <- dat |>
      dplyr::left_join(PoO_df |> dplyr::select("type", "matOrg", "patOrg"), by = "type") |>
      dplyr::mutate(count = base::replace(.data$count, is.na(.data$patOrg), M.count.E0),
                    matOrg = base::replace(.data$matOrg, is.na(.data$matOrg), 1),
                    patOrg = base::replace(.data$patOrg, is.na(.data$patOrg), 0)) |>
      dplyr::add_row(dat |>
                       dplyr::filter(dplyr::row_number() %in% heteroInds.E0) |>
                       dplyr::mutate(count = dat$count[heteroInds.E0] - M.count.E0,
                                     matOrg = 0,
                                     patOrg = 1)) |>
      dplyr::arrange(dplyr::desc(.data$D), dplyr::desc(.data$E), .data$type, dplyr::desc(.data$matOrg)) |>
      dplyr::mutate(typeOrig = base::rep(1:16, dplyr::n()/16),
                    Im = .data$matOrg * .data$D,
                    If = .data$patOrg * .data$D) |>
      dplyr::relocate("typeOrig")
  }
  return(PoO_dat)
}

add_PoO_data_15 <- function(dat) {
  # Portion of model equation and offset depends on mating type model
  heteroInds <- with(dat, which(type == 9))
  PoO_dat <- dat |>
    dplyr::left_join(PoO_df, by = c("M", "F", "C")) |>
    dplyr::mutate(matOrg = replace(.data$matOrg, is.na(.data$matOrg), 0),
                  patOrg = replace(.data$patOrg, is.na(.data$patOrg), 0)) |>
    dplyr::mutate(Im = .data$matOrg * .data$D,
                  If = .data$patOrg * .data$D) |>
    dplyr::arrange(dplyr::desc(.data$D), dplyr::desc(.data$E), .data$type, dplyr::desc(.data$matOrg))
  return(PoO_dat)
}

## elements of covariance matrix function, a = row index, b = col index
## mu = complete poisson means
## L must be a matrix of 0s and 1s with no more than one 1 per column
cov_trill <- function(y, L, a, b, mu){
  if (a == b) {
    if (1 %in% L[,a]) {
      r_a <- which(L[,a] == 1)
      return(y[[r_a]] *
               mu[[a]] / (t(L[r_a,]) %*% mu) *
               (1 - mu[[a]] / (t(L[r_a,]) %*% mu))
      )
    } else {
      return(mu[[a]])
    }
  } else {
    r_a <- which(L[,a] == 1)
    r_b <- which(L[,b] == 1)
    if(r_a == r_b){
      return(- y[[r_a]] *
               mu[[a]] / (t(L[r_a,]) %*% mu) *
               mu[[b]] / (t(L[r_b,]) %*% mu)
      )
    } else {
      return(0)
    }
  }
}

## Create df with rows for every possible trio
createGenoMat <- function() {
  M <- c(rep(0, 3), rep(1, 2), 0, 2, rep(1, 6), rep(2, 3))
  F <- c(0, rep(1, 2), rep(0, 2), 2, 0, rep(1, 4), rep(2, 2), rep(1, 2), 2)
  C <- c(rep(0, 2), 1, 0, rep(1, 3), 0, rep(1, 2), 2, 1, 2, 1, rep(2, 2))

  ## Indicator variables for parent-of-origin
  matOrg <- c(rep(0, 4), 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1)
  patOrg <- c(rep(0, 2), 1, rep(0, 2), 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1)

  return(data.frame(M, F, C, matOrg, patOrg))
}

createHaplinGeno <- function() {
  M <- c(rep("1;1", 3), rep("1;2", 2), "1;1", "2;2", rep("1;2", 6), rep("2;2", 3))
  F <- c("1;1", rep("1;2", 2), rep("1;1", 2), "2;2", "1;1", rep("1;2", 4), rep("2;2", 2), rep("1;2", 2), "2;2")
  C <- c(rep("1;1", 2), "1;2", "1;1", rep("1;2", 3), "1;1", rep("1;2", 2), "2;2", "1;2", "2;2", "1;2", rep("2;2", 2))
  return(data.frame(M, F, C))
}

mtmat <- function(maf = 0.4, C = c(1, 1, 1)) {
  genocat <- createGenoMat()
  mts <- unique(genocat[, 1:2])

  Mprobs <- stats::dbinom(mts[, 1], 2, prob = maf)
  Fprobs <- stats::dbinom(mts[, 2], 2, prob = maf)

  asymfactor <- c(1, 2 - C[1], C[1],
                  2 - C[2], C[2], 1,
                  2 - C[3], C[3], 1)
  mts <- cbind(mts, prMF = Mprobs * Fprobs * asymfactor)
  ## First three mating pairs: can make 0 or 1
  ## Middle three mating pairs: can make 0, 1, or 2
  ## Last three mating pairs: Can make 1 or 2

  mt_MS <- c(1, 2, 2, 3, 3, 4, 5, 5, 6)
  mt_MaS <- 1:9

  mts <- cbind(mt_MS, mt_MaS, mts)

  prCGivenMFOrg <-
    c(1,
      rep(1 / 2, 2),
      rep(1 / 2, 2),
      1,
      1,
      rep(1 / 4, 4),
      rep(1 / 2, 2),
      rep(1 / 2, 2),
      1)

  genocat <- cbind(genocat, prCGivenMFOrg)
  temp <- merge(genocat, mts)
  temp$prMFCOrg <- temp$prMF * temp$prCGivenMFOrg

  return(temp[
    order(temp$mt_MaS),
    c("mt_MS", "mt_MaS", "M", "F", "C", "matOrg", "patOrg", "prMF", "prCGivenMFOrg", "prMFCOrg")
  ])
}

# A simplified function for simulating a subset of the full data
# under particular conditions
# mtCoef = C1, C2, C4 in paper, used like this since it preserves HWE for tests
# Shoud work for imprinting & controls/environment (controls are not affected by imprinting)
# - Case Trios
# - Control Trios
# - With E or without
# - With MaS or not
# - With Imprinting or not (controls not required, for now only consider no controls or no MaS)
simulateDataSubset <- function(ntrios = 1000, maf = 0.3,
                               R = c(1, 1, 1), S = c(1, 1, 1),
                               mtCoef = c(1, 1, 1),
                               V = c(1, 1, 1), includeE = FALSE, Einteraction = "M",
                               includeControl = FALSE, Im = 1, If = 1) {
  genomat <- mtmat(maf, C = mtCoef)
  # Compute the values proportional to P(D|M,F,C,E) for the
  # selected model.

  # Add main genetic effects
  diseasefactor <- R[(genomat$C + 1)] * S[(genomat$M + 1)]

  # Add imprinting effects
  impFactorm <- ifelse(genomat$matOrg, Im, 1)
  impFactorf <- ifelse(genomat$patOrg, If, 1)
  diseasefactor <- diseasefactor * impFactorm * impFactorf

  # Add environmental effects for trios that will have E=1
  if (Einteraction == "M") {
    diseasefactor <- diseasefactor * V[(genomat$M + 1)]
  } else if (Einteraction == "C") {
    diseasefactor <- diseasefactor * V[(genomat$C + 1)]
  } else if (Einteraction == "Im") {
    diseasefactor <- diseasefactor * V[(genomat$matOrg + 1)]
  } else if (Einteraction == "If") {
    diseasefactor <- diseasefactor * V[(genomat$patOrg + 1)]
  }

  # Compute the values proportional to P(M,F,C,E|D) = P(D|M,F,C,E)P(M,F)P(C|M,F)
  # and re-scale so that probabilities sum to 1 (if only cases this happens)
  diseaseprob <- genomat$prMFC * diseasefactor
  diseaseprob <- diseaseprob / sum(diseaseprob)


  # Sample the trios
  triocounts <- sample(1:length(diseaseprob), ntrios, prob = diseaseprob, replace = TRUE)
  counts <- data.frame(table(triocounts))
  colnames(counts) <- c("typeOrig", "count")

  # Prepare data for R log-linear
  genomat <- data.frame(typeOrig = 1:length(diseaseprob), genomat, PrCMFEgivenD = diseaseprob)
  triodat <- merge(genomat, counts, all = TRUE, by = "typeOrig")
  triodat[is.na(triodat[, "count"]), "count"] <- 0

  if (includeE == TRUE) {
    E <- rep(1, nrow(triodat))
  } else {
    E <- rep(0, nrow(triodat))
  }

  if (includeControl == TRUE) {
    D <- rep(0, nrow(triodat))
  } else {
    D <- rep(1, nrow(triodat))
  }

  triodat <- data.frame(triodat[, 1:6], E = E, D = D, triodat[, 7:ncol(triodat)])
  subdat <- data.frame(triodat[, c(1:8)], count = triodat$count) |>
    dplyr::group_by(.data$mt_MS, .data$mt_MaS, .data$M, .data$F, .data$C, .data$E, .data$D) |>
    dplyr::summarise_at(dplyr::vars("count"), sum) |>
    dplyr::ungroup() |>
    dplyr::mutate(type = 1:15) |>
    dplyr::relocate("type") |>
    as.data.frame()

  # Prepare data for Haplin. Haplin has 1 row per trio, with the row
  # consisting of M,F,C genotype columns in format 1;1, 1;2 or 2;2
  # The first column will be the environmental covariate
  # The second column will be the disease status
  # The third, fourth and fifth columns are the genotype for M, F, and C
  haplingeno <- createHaplinGeno()
  # if (includeControl==TRUE){ # Only parents so make child genotype NA
  #  MF=unique(haplingeno[,1:2])
  #  haplingeno=data.frame(type=1:nrow(MF),MF,C=rep(NA,nrow(MF)))
  # } else{
  haplingeno <- data.frame(typeOrig = 1:16, haplingeno)
  # }
  haplindat <- merge(haplingeno, triodat, by = "typeOrig") |>
    dplyr::group_by(.data$M.x, .data$F.x, .data$C.x, .data$mt_MS, .data$mt_MaS, .data$E, .data$D, .data$prMF) |>
    dplyr::summarise_at(dplyr::vars("count"), sum) |>
    dplyr::ungroup() |>
    dplyr::mutate(type = 1:15) |>
    dplyr::relocate("type") |>
    as.data.frame()

  finaldat <- NULL
  for (j in 1:nrow(haplindat)) {
    if (haplindat[j, "count"] > 0) {
      genos <- matrix(unlist(rep(haplindat[j, 2:4], haplindat[j, "count"])), ncol = 3, byrow = TRUE)
      genos <- cbind(
        rep(haplindat[j, "E"], length = haplindat[j, "count"]),
        rep(haplindat[j, "D"], length = haplindat[j, "count"]),
        genos
      )
      finaldat <- rbind(finaldat, genos)
    }
  }
  finaldat <- as.data.frame(finaldat)

  # If D=0, remove the genotypes for the children since these
  # are not used for log-linear and EMIM
  # if (includeControl==TRUE){
  #  subdat$C=rep(-9,nrow(subdat))
  # }

  return(list(dat4R = subdat, dat4haplin = finaldat, datFull = triodat))
}




# To be used when there is hidden population substructure. This
# function merges the count column for two hidden populations.
mergeCounts <- function(dat1, dat2) {
  triodat1 <- dat1$dat4R
  triodat2 <- dat2$dat4R

  # Merge for log linear analysis
  alltriodat <- merge(triodat1, triodat2, by = "type", all = TRUE)
  counts <- alltriodat$count.x + alltriodat$count.y
  alltriodat <- data.frame(triodat1[, c(1:8)], count = counts)

  # Stack for haplin
  allhaplindat <- rbind(dat1$dat4haplin, dat2$dat4haplin)


  return(list(dat4R = alltriodat, dat4haplin = allhaplindat))
}

# To be used with environmental effects and control trios where
# the datasets with E=0 or 1 and D=0 or 1 are stacked.
stackCounts <- function(dat1, dat2) {
  triodat1 <- dat1$dat4R
  triodat2 <- dat2$dat4R

  # Stack for log linear analysis
  alltriodat <- rbind(triodat1, triodat2)

  # Stack for haplin
  allhaplindat <- rbind(dat1$dat4haplin, dat2$dat4haplin)

  return(list(dat4R = alltriodat, dat4haplin = allhaplindat))
}


# Create a ped file
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype
createPed <- function(dat) {
  # Create the fixed columns
  ntrios <- sum(dat[, "count"])

  famid <- as.vector(t(matrix(rep(1:ntrios, 3), ncol = 3, byrow = F)))
  indid <- seq(1:length(famid))
  pid <- rep(0, length(famid))
  pid[3 * (1:ntrios)] <- indid[3 * (1:ntrios) - 2]
  mid <- rep(0, length(famid))
  mid[3 * (1:ntrios)] <- indid[3 * (1:ntrios) - 1]
  sex <- as.vector(rbind(
    rep(1, ntrios),
    rep(2, ntrios),
    sample(c(1, 2), ntrios, replace = TRUE)
  ))


  # Create the columns with data that changes from row
  # to row
  finaldat <- NULL

  for (i in 1:nrow(dat)) {
    ntrios <- dat[i, "count"]

    if (ntrios > 0) {
      phenotype <- rep(c(-9, -9, (dat[i, "D"] + 1)), ntrios)
      E <- rep(dat[i, "E"], 3 * ntrios)
      genotype1 <- unlist(rep(dat[i, c("F", "M", "C")], ntrios))
      genotype1[genotype1 == 2] <- 1
      genotype2 <- unlist(rep(dat[i, c("F", "M", "C")], ntrios)) - 1
      genotype2[genotype2 == -1] <- 0
      genotype1 <- genotype1 + 1
      genotype2 <- genotype2 + 1

      # Deal with missing genotypes. They are coded -9, so will be changed
      # to -8 by above code. Use "0" for missing
      # genotype1[genotype1==-8]=0
      # genotype2[genotype2==-8]=0

      tempdat <- data.frame(E, D = phenotype, genotype1, genotype2)

      finaldat <- rbind(finaldat, tempdat)
    }
    tempdat <- NULL
  }

  finaldat <- data.frame(famid, indid, pid, mid, sex, finaldat)

  return(finaldat)
}

sum_grp <- function(x, grp){
  replace_inds <- unlist(lapply(grp, min))
  remove_inds <- setdiff(unlist(grp), replace_inds)
  sums <-
    unlist(
      lapply(
        grp,
        function(y) sum(x[y])
      )
    )
  x[replace_inds] <- sums
  if (length(remove_inds) > 0){
    x <- x[-remove_inds]
  }
  names(x) <- 1:length(x)
  return(x)
}

mean_grp <- function(x, grp){
  replace_inds <- unlist(lapply(grp, min))
  remove_inds <- setdiff(unlist(grp), replace_inds)
  means <-
    unlist(
      lapply(
        grp,
        function(y) mean(x[y])
      )
    )
  x[replace_inds] <- means
  if (length(remove_inds) > 0){
    x <- x[-remove_inds]
  }
  names(x) <- 1:length(x)
  return(x)
}
