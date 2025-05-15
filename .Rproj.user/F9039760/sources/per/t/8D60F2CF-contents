add_PoO_data <- function(dat, Mprop, includeE) {
  # Portion of model equation and offset depends on mating type model
  if (includeE && length(Mprop) != 2){
    stop("Environmental interactions included hence Mprop should be a vector of length two!")
  }
  heteroInds.E0 <- with(dat, which(type == 9 & E == 0))
  M.count.E0 <- dat %>%
    dplyr::filter(type == 9, E == 0) %>%
    dplyr::pull(count) %>%
    magrittr::multiply_by(Mprop[[1L]])
  if (includeE){
    heteroInds.E1 <- with(dat, which(type == 9 & E == 1))
    M.count.E1 <- dat %>%
      dplyr::filter(type == 9, E == 1) %>%
      dplyr::pull(count) %>%
      magrittr::multiply_by(Mprop[[2L]])

    PoO_dat <- dat %>%
      dplyr::left_join(PoO_df %>% dplyr::select(type, matOrg, patOrg), by = "type") %>%
      dplyr::mutate(count = base::replace(count, is.na(patOrg) & E == 0, M.count.E0),
                    count = base::replace(count, is.na(patOrg) & E == 1, M.count.E1),
                    matOrg = base::replace(matOrg, is.na(matOrg), 1),
                    patOrg = base::replace(patOrg, is.na(patOrg), 0)) %>%
      dplyr::add_row(dat %>%
                       dplyr::filter(dplyr::row_number() %in% heteroInds.E0) %>%
                       dplyr::mutate(count = dat$count[heteroInds.E0] - M.count.E0,
                                     matOrg = 0,
                                     patOrg = 1)) %>%
      dplyr::add_row(dat %>%
                       dplyr::filter(dplyr::row_number() %in% heteroInds.E1) %>%
                       dplyr::mutate(count = dat$count[heteroInds.E1] - M.count.E1,
                                     matOrg = 0,
                                     patOrg = 1)) %>%
      dplyr::arrange(dplyr::desc(D), dplyr::desc(E), type, dplyr::desc(matOrg)) %>%
      dplyr::mutate(typeOrig = base::rep(1:16, dplyr::n()/16),
                    Im = matOrg * D,
                    If = patOrg * D) %>%
      dplyr::relocate(typeOrig)
  } else {
    PoO_dat <- dat %>%
      dplyr::left_join(PoO_df %>% dplyr::select(type, matOrg, patOrg), by = "type") %>%
      dplyr::mutate(count = base::replace(count, is.na(patOrg), M.count.E0),
                    matOrg = base::replace(matOrg, is.na(matOrg), 1),
                    patOrg = base::replace(patOrg, is.na(patOrg), 0)) %>%
      dplyr::add_row(dat %>%
                       dplyr::filter(dplyr::row_number() %in% heteroInds.E0) %>%
                       dplyr::mutate(count = dat$count[heteroInds.E0] - M.count.E0,
                                     matOrg = 0,
                                     patOrg = 1)) %>%
      dplyr::arrange(dplyr::desc(D), dplyr::desc(E), type, dplyr::desc(matOrg)) %>%
      dplyr::mutate(typeOrig = base::rep(1:16, dplyr::n()/16),
                    Im = matOrg * D,
                    If = patOrg * D) %>%
      dplyr::relocate(typeOrig)
  }
  return(PoO_dat)
}

add_PoO_data_15 <- function(dat) {
  # Portion of model equation and offset depends on mating type model
  heteroInds <- with(dat, which(type == 9))
  PoO_dat <- dat %>%
    dplyr::left_join(PoO_df, by = c("M", "F", "C")) %>%
    dplyr::mutate(matOrg = replace(matOrg, is.na(matOrg), 0),
                  patOrg = replace(patOrg, is.na(patOrg), 0)) %>%
    dplyr::mutate(Im = matOrg * D,
                  If = patOrg * D) %>%
    dplyr::arrange(desc(D), desc(E), type, desc(matOrg))
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

  Mprobs <- dbinom(mts[, 1], 2, prob = maf)
  Fprobs <- dbinom(mts[, 2], 2, prob = maf)

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
  subdat <- data.frame(triodat[, c(1:8)], count = triodat$count) %>%
    dplyr::group_by(mt_MS, mt_MaS, M, F, C, E, D) %>%
    dplyr::summarise_at(dplyr::vars(count), sum) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(type = 1:15) %>%
    dplyr::relocate(type) %>%
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
  haplindat <- merge(haplingeno, triodat, by = "typeOrig") %>%
    dplyr::group_by(M.x, F.x, C.x, mt_MS, mt_MaS, E, D, prMF) %>%
    dplyr::summarise_at(dplyr::vars(count), sum) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(type = 1:15) %>%
    dplyr::relocate(type) %>%
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
