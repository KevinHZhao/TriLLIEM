## Like V2 but HW, PS and MaS tests added


MaSlrt=function(dat){

  subdat=dat[c(2:7,11:14),]
  subdat=aggregate(subdat$count,by=list(M=subdat$M,F=subdat$F),sum)
  subdat=subdat[c(1,3,2,5,4,6),]
  colnames(subdat)[3]="count"

  num=3*(subdat$count[5]+subdat$count[6])+2*(subdat$count[3]+subdat$count[4])+
    (subdat$count[1]+subdat$count[2])
  denom=4*sum(subdat$count)
  p=num/denom

  ishet.M=subdat$M==1
  ishet.F=subdat$F==1
  factor.M=rep(1,nrow(subdat))
  factor.M[ishet.M]=2
  factor.F=rep(1,nrow(subdat))
  factor.F[ishet.F]=2
  HWfreqs=factor.M*p^subdat$M*(1-p)^(2-subdat$M)*factor.F*p^subdat$F*(1-p)^(2-subdat$F)


  c4=2*subdat$count[5]/(subdat$count[5]+subdat$count[6])
  c2=2*subdat$count[3]/(subdat$count[3]+subdat$count[4])
  c1=2*subdat$count[1]/(subdat$count[1]+subdat$count[2])
  Cfactor=c(c1,2-c1,c2,2-c2,c4,2-c4)

  l1=sum(subdat$count*log(HWfreqs*Cfactor))
  l0=sum(subdat$count*log(HWfreqs))
  LRT=2*(l1-l0)
  return(pchisq(LRT,df=3,lower=F))
}


getHWE=function(dat){
  require("HardyWeinberg")
  mom=tapply(dat$count, dat$M, sum)
  dad=tapply(dat$count, dat$F, sum)
  both=dad+mom
  res=c(Mom.exact=HWExact(mom, verbose=F)$pval,
           Dad.exact=HWExact(dad, verbose=F)$pval,
           Both.exact=HWExact(both, verbose=F)$pval)
  return(res)
}



## Run a loglinear analysis with the specified model
runLoglin=function(mtmodel, effects, dat, PStest=FALSE)
{

  # Portion of model equation and offset depends on mating type model
  dat$offset=rep(NA,nrow(dat))
  if (mtmodel=="HW"){

    dat$HWgeno=dat[,"M"]+dat[,"F"]
    mteffect="HWgeno"

    if (sum(dat$E==1)>0){ # E variable in dataset

      # Assumes order type=1:15 and controls children are in dataset
      dat$offset[dat$E==0]=c(rep(1,8),2,rep(1,6))
      dat$offset[dat$E==1]=c(rep(1,8),2,rep(1,6))

    } else {
      # Assumes order type=1:15 and controls children are in dataset
      dat$offset=c(rep(1,8),2,rep(1,6))
    }

    modelformula="count~" # Must include intercept for HW model because of log(1-p) term

  } else if (mtmodel=="MS"){
    mteffect="as.factor(mt_MS)"

    if (sum(dat$E==1)>0 ) { # Includes E in dataset
      # Assumes order type=1:15 and controls children are in dataset
      dat$offset[dat$E==1]=c(rep(1,8),2,rep(1,6))
      dat$offset[dat$E==0]=c(rep(1,8),2,rep(1,6))
    } else {
      # Assumes order type=1:15 and controls children are in dataset
      dat$offset=c(rep(1,8),2,rep(1,6))
    }

    modelformula="count~-1+" # I think I can remove the intercept for MS model
  } else if (mtmodel=="MaS"){

    if (length(unique(dat$D))==1){
      stop("Only 1 phenotype in the phenotype column. Mating asymmetry models require \n
            both cases and controls\n")
    }
    mteffect="as.factor(mt_MaS)"

    if (sum(dat$E==1)>0){ # Includes environmental variable
      dat$offset[dat$E==0]=c(rep(1,8),2,rep(1,6))
      dat$offset[dat$E==1]=c(rep(1,8),2,rep(1,6))
    } else { # No environmental variable
      # Assumes order type=1:15 (D=1)
      dat$offset=c(rep(1,8),2,rep(1,6))
    }

    modelformula="count~-1+" # I think I can remove intercept
  }

  if (is.element("E:M",effects)){

    if (sum(dat$D==0)>0){ # There are controls; include environmental interaction
      dat$C=dat$C*dat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
      dat$M=dat$M*dat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
                        # (Note that the E:M term in model will be from crossing this
                        # variable with E=1, which is exactly what is needed.
      if (mtmodel=="HW"){ #Include main effect of E to give different intercept for HW+E case
        modeleffects=c(mteffect,paste0(mteffect,":E"),effects,"E","D","E:D")
      } else {
        modeleffects=c(mteffect,paste0(mteffect,":E"), effects, "D","E:D")
      }

    } else { # No controls

      if (mtmodel=="HW"){ #Include main effect of E to give different intercept for HW+E case
        modeleffects=c(mteffect,paste0(mteffect,":E"),effects,"E")
      } else {
        modeleffects=c(mteffect,paste0(mteffect,":E"),effects)
      }
    }

  } else { # No E:M effect

    if (sum(dat$D==0)>0){
      dat$C=dat$C*dat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
      dat$M=dat$M*dat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 O
      modeleffects=c(mteffect, effects, "D")



    } else {
      modeleffects=c(mteffect,effects)
    }
  }

  linpred=paste(modeleffects, collapse="+")
  modelformula=paste0(modelformula,linpred)

  # Setup results objects
  resVec=vector(length=length(effects))
  names(resVec)=effects
  pvalVec=vector(length=length(effects))
  names(pvalVec)=effects
  resVecPS=NULL
  pvalVecPS=NULL

  # Include test and results under population stratification
  if (PStest==TRUE){

    resVecPS=vector(length=length(effects))
    names(resVecPS)=effects
    pvalVecPS=vector(length=length(effects))
    names(pvalVecPS)=effects

    if (sum(dat$D==0)==0){
      stop("Can only test for population stratification if there are control trios\n")
    } else{
      PSeffect=paste0(mteffect,":D")
      modelformula.PS=paste(modelformula,PSeffect,sep="+")
    }
  }


  # Run model and save results
  res=glm(as.formula(modelformula), data=dat, offset=log(dat$offset), family=poisson())

  # R is not consistent about how interaction is specified. Even though it
  # is fit as E:M, sometimes R flips it to M:E in the output of results.
  if (is.element("M:E",names(coef(res)))){
    effects[effects=="E:M"]="M:E"
  }

  for (j in 1:length(effects)){
    resVec[j]=exp(summary(res)$coef[effects[j],1])
    pvalVec[j]=summary(res)$coef[effects[j],4]
  }

  test.res=NULL
  if (PStest==TRUE){
    res.PS=glm(as.formula(modelformula.PS), data=dat, offset=log(dat$offset), family=poisson())
    test.res=anova(res,res.PS, test="LRT")

    for (j in 1:length(effects)){
      resVecPS[j]=exp(summary(res.PS)$coef[effects[j],1])
      pvalVecPS[j]=summary(res.PS)$coef[effects[j],4]
    }
  }

  return(res)
  # return(list(effects=resVec, pvals=pvalVec, PS.test=test.res$`Pr(>Chi)`[2],
  #             effectsPS=resVecPS, pvalsPS=pvalVecPS))
}

## Run a haplin analysis with the specified model
runHaplin=function(effects, dat, haplinControls=FALSE)
{

  require("Haplin")
  #browser()

  # Setup results objects
  if (is.element("E:M",effects)){
    neweffects=c(setdiff(effects,c("E:M","M")),"M[E=0]", "M[E=1]","M[E=1]/M[E=0]")
    resVec=vector(length=length(neweffects))
    pvalVec=vector(length=length(effects))
    names(resVec)=neweffects
    names(pvalVec)=effects
  } else {
    resVec=vector(length=length(effects))
    names(resVec)=effects
    pvalVec=vector(length=length(effects))
    names(pvalVec)=effects
  }

  # Number of columns that are not genotype. First column is E, second is Phenotype
  nvars=2

  # Write out the data for haplin, then read it in using haplin
  write.table(dat,"haplin_temp.dat",quote=F,row=F,col=F,sep=" ")
  write(c("chrom snp a","1 rs1 0"), "temp.map",ncol=1)
  dat.raw=invisible(genDataRead("haplin_temp.dat",
                                format="haplin",overwrite=TRUE, n.vars=nvars,
                                map.file="temp.map"))

  if (is.element("E:M",effects)){

    if (haplinControls){
      dat.processed=invisible(genDataPreprocess( data.in =dat.raw, design = "cc.triad",
                                                 overwrite=TRUE, map.file="temp.map"))
      res=invisible(haplinStrat(dat.processed, response="mult", verbose=FALSE,
                                design = "cc.triad", ccvar=2,
                                printout=FALSE, strata=1,reference = "ref.cat", maternal=TRUE))
    } else {
      dat.processed=invisible(genDataPreprocess( data.in =dat.raw, design = "triad",
                                                 overwrite=TRUE, map.file="temp.map"))
      res=invisible(haplinStrat(dat.processed, response="mult", verbose=FALSE,
                                printout=FALSE, strata=1,reference = "ref.cat", maternal=TRUE))
    }
  } else {

    if (is.element("M",effects)){
      includeMaternal=TRUE
    } else {includeMaternal=FALSE}

    if (haplinControls==TRUE){
      dat.processed=invisible(genDataPreprocess(data.in =dat.raw, design = "cc.triad",
                                                overwrite=TRUE, map.file="temp.map"))
      res=invisible(haplin(dat.processed, response="mult", verbose=FALSE, design="cc.triad",
                           ccvar=2, printout=FALSE, reference = "ref.cat", maternal=includeMaternal))
    } else {
      dat.processed=invisible(genDataPreprocess(data.in =dat.raw, design = "triad",
                                                overwrite=TRUE, map.file="temp.map"))
      res=invisible(haplin(dat.processed, response="mult", verbose=FALSE,
                           printout=FALSE, reference = "ref.cat", maternal=includeMaternal))
    }

    res=haptable(res)[2,]
    if (is.element("M",effects)){
      resVec["M"]=res$RRm.est.
      pvalVec["M"]=res$RRm.p.value
    }
    if (is.element("C",effects)){
      resVec["C"]=res$RR.est.
      pvalVec["C"]=res$RR.p.value
    }
  }

  system("rm temp.map haplin_temp.dat")
  return(list(effects=resVec, pvals=pvalVec))
}


## Run a haplin analysis with the specified model
runEMIM=function(mtmodel, effects, peddat, emimpath="~/StatGenTools/emim-v3.22-code/")
{

  # Setup results objects
  if (is.element("E:M",effects)){
    neweffects=c(setdiff(effects,c("E:M","M")),"M[E=0]", "M[E=1]","M[E=1]/M[E=0]")
    resVec=vector(length=length(neweffects))
    pvalVec=vector(length=length(effects))
    names(resVec)=neweffects
    names(pvalVec)=effects
  } else {
    resVec=vector(length=length(effects))
    names(resVec)=effects
    pvalVec=vector(length=length(effects))
    names(pvalVec)=effects
  }


  # Set up options for the parameter file
  options=" -a -so"

  if (is.element("C",effects)){ # Multiplicative allele model for C effect
    options=c(options,"-ct")
  }

  if (is.element("M",effects)||is.element("E:M",effects)){ # Multiplicative allele model for M effect
    options=c(options,"-mt")
  }
  options=paste0(options," ",collapse=" ")


  # Note that ped data includes an E column, which must
  # must be removed
  if (is.element("E:M",effects)){

    ## Subset the data into cases where E=0 and once with E=1
    subset0=subset(peddat,peddat$E==0)
    subset1=subset(peddat,peddat$E==1)

    # All
    peddat=peddat[,c("famid","indid","pid","mid","sex","D","genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_all.ped",col.names=F,row.names=F,quote=F)
    write(c("1","A","0","0"), "temp_pedigree_all.map",ncol=4)

    # E=0 data
    peddat=subset0[,c("famid","indid","pid","mid","sex","D","genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_0.ped",col.names=F,row.names=F,quote=F)
    write(c("1","A","0","0"), "temp_pedigree_0.map",ncol=4)

    # E=1 data
    peddat=subset1[,c("famid","indid","pid","mid","sex","D","genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_1.ped",col.names=F,row.names=F,quote=F)
    write(c("1","A","0","0"), "temp_pedigree_1.map",ncol=4)

    ## Run PREMIM and EMIM
    for (i in c("all",0,1)){

      # Run PREMIM
      command=paste0(emimpath,"premim",options,"temp_pedigree_",i,".ped temp_pedigree_",i,".map")
      system(command, intern=TRUE)

      # Change parameter in options file for mating symmetry
      if (mtmodel=="MS"){
        params=readLines("emimparams.dat")
        params[16]="0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
        write(params,"emimparams.dat",ncol=1)
      }
      if (mtmodel=="MaS"){
        params=readLines("emimparams.dat")
        params[16]="0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
        params[18]="1   << use CPG likelihood (estimate 9 mu parameters)"
        write(params,"emimparams.dat",ncol=1)
      }

      # run EMIM and rename results
      command=paste0(emimpath,"emim")
      system(command, intern=TRUE)
      system(paste0("mv emimsummary.out emimsummary_",i,".out"))
    }

    ## Read in results and parse
    resAll=read.table("emimsummary_all.out", header=T)
    res0=read.table("emimsummary_0.out", header=T)
    res1=read.table("emimsummary_1.out", header=T)

    if (is.element("C",effects)){
      resVec["C"]=exp(resAll$lnR1)
      pvalVec["C"]=2*pnorm(abs(resAll$lnR1/resAll$sd_lnR1),lower=F)
    }
    pvalVec["M"]=2*pnorm(abs(resAll$lnS1/resAll$sd_lnS1),lower=F)

    resVec["M[E=0]"]=exp(res0$lnS1)
    resVec["M[E=1]"]=exp(res1$lnS1)
    resVec["M[E=1]/M[E=0]"]=resVec["M[E=1]"]/resVec["M[E=0]"]

    # Get a Wald-type GE test like Haplin
    z=abs(res0$lnS1-res1$lnS1)/sqrt(res0$sd_lnS1^2+res1$sd_lnS1^2)
    pvalVec["E:M"]=2*pnorm(z,lower=F)

  } else {

    peddat=peddat[,c("famid","indid","pid","mid","sex","D","genotype1", "genotype2")]

    write.table(peddat, "temp_pedigree.ped",col.names=F,row.names=F,quote=F)
    write(c("1","A","0","0"), "temp_pedigree.map",ncol=4)


    # Run PREMIM
    command=paste0(emimpath,"premim",options,"temp_pedigree.ped temp_pedigree.map")
    system(command, intern=TRUE)

    # Change parameter for mating symmetry
    if (mtmodel=="MS"){
      params=readLines("emimparams.dat")
      params[16]="0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
      write(params,"emimparams.dat",ncol=1)
    }
    if (mtmodel=="MaS"){
      params=readLines("emimparams.dat")
      params[16]="0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
      params[18]="1   << use CPG likelihood (estimate 9 mu parameters)"
      write(params,"emimparams.dat",ncol=1)
    }

    # RUN EMIM
    command=paste0(emimpath,"emim")
    system(command, intern=TRUE)

    # Read in results
    res=read.table("emimsummary.out", header=T)

    if (is.element("M",effects)){
      resVec["M"]=exp(res$lnS1)
      pvalVec["M"]=2*pnorm(abs(res$lnS1/res$sd_lnS1),lower=F)
    }
    if (is.element("C",effects)){
      resVec["C"]=exp(res$lnR1)
      pvalVec["C"]=2*pnorm(abs(res$lnR1/res$sd_lnR1),lower=F)
    }
  }

  system("rm temp_pedigree*")
  return(list(effects=resVec, pvals=pvalVec))
}



# Run a simulation with the specified parameters.
# eligible effects are "M","C", "M:E", "D" (terms included in model)
# eligible MT models are "HW", "MS", "MaS"
# includeE, PopStrat, and Control means that these conditions are added to
#   the dataset. However, to include an environmental effect to the analysis
#   must have V not be all 0 and to include controls in analysis than a "D"
#   term will need to be added to the effects.
runTrioReps=function(nreps=100, ntrios=1000, maf=0.3,
                     R=c(1,1,1), S=c(1,1,1),  V=c(1,1,1),
                     mtCoef=c(1,1,1), effects=c("M","C"),
                     includeE=FALSE, envint="Mother", prE=0,
                     includePopStrat=FALSE, numPop=1, Fst=0.005, prCase.byPop=NULL, prControl.byPop=NULL,
                     includeControl=FALSE,prControl=0, prE.control=prE,
                     Loglin=TRUE, loglinMT=c("HW","MS"),
                     Haplin=TRUE, haplinControls=FALSE,
                     EMIM=TRUE, EMIMMT=c("HW","MS"),
                     PStest=FALSE, HWtest=FALSE, MStest=FALSE){


  # Set up results tables. There will be two tables for each method that
  # will be run. There may be objects for HW, mating symmetry and population
  # stratification tests.
  pvalPS=NULL
  resLoglinPS=NULL
  pvalLoglinPS=NULL
  resLoglin=NULL
  pvalLoglin=NULL
  resHaplin=NULL
  pvalHaplin=NULL
  resEMIM=NULL
  pvalEMIM=NULL
  HWres=NULL
  MSres=NULL

  if (Loglin){
    ncolLoglin=length(effects)*length(loglinMT)
    temp=expand.grid(effects, loglinMT)[,2:1]
    namesLoglin=paste(temp[,1],temp[,2], sep=".")
    resLoglin=matrix(nrow=nreps,ncol=ncolLoglin)
    colnames(resLoglin)=namesLoglin
    pvalLoglin=matrix(nrow=nreps,ncol=ncolLoglin)
    colnames(pvalLoglin)=namesLoglin

    if (PStest==TRUE){
      if (prControl==0){
        stop("PStest only possible when there are controls\n")
      } else {
        pvalPS=matrix(nrow=nreps,ncol=length(loglinMT))
        colnames(pvalPS)=loglinMT
        resLoglinPS=matrix(nrow=nreps,ncol=ncolLoglin)
        colnames(resLoglinPS)=namesLoglin
        pvalLoglinPS=matrix(nrow=nreps,ncol=ncolLoglin)
        colnames(pvalLoglinPS)=namesLoglin
      }
    }
  }

  if (Haplin){

    if (is.element("E:M",effects)){
      neweffects=c(setdiff(effects,c("E:M","M")),"M[E=0]", "M[E=1]","M[E=1]/M[E=0]")
      ncolHaplin=length(neweffects)
      resHaplin=matrix(nrow=nreps,ncol=ncolHaplin)
      colnames(resHaplin)=neweffects
      pvalHaplin=matrix(nrow=nreps,ncol=length(effects))
      colnames(pvalHaplin)=effects
    } else {
      ncolHaplin=length(effects)
      namesHaplin=effects
      resHaplin=matrix(nrow=nreps,ncol=ncolHaplin)
      colnames(resHaplin)=namesHaplin
      pvalHaplin=matrix(nrow=nreps,ncol=ncolHaplin)
      colnames(pvalHaplin)=namesHaplin
    }

  }

  if (EMIM){

    if (is.element("E:M",effects)){
      neweffects=c(setdiff(effects,c("E:M","M")),"M[E=0]", "M[E=1]","M[E=1]/M[E=0]")
      ncolEMIM=length(EMIMMT)*length(neweffects)
      temp.pval=expand.grid(effects, EMIMMT)[,2:1]
      temp.res=expand.grid(neweffects, EMIMMT)[,2:1]
      resEMIM=matrix(nrow=nreps,ncol=ncolEMIM)
      colnames(resEMIM)=paste(temp.res[,1],temp.res[,2], sep=".")
      pvalEMIM=matrix(nrow=nreps,ncol=nrow(temp.pval))
      colnames(pvalEMIM)=paste(temp.pval[,1],temp.pval[,2], sep=".")
    } else {
      ncolEMIM=length(EMIMMT)*length(effects)
      temp=expand.grid(effects, EMIMMT)[,2:1]
      namesEMIM=paste(temp[,1],temp[,2], sep=".")
      resEMIM=matrix(nrow=nreps,ncol=ncolEMIM)
      colnames(resEMIM)=namesEMIM
      pvalEMIM=matrix(nrow=nreps,ncol=ncolEMIM)
      colnames(pvalEMIM)=namesEMIM
    }

  }

  if (HWtest==TRUE){
    if (includeControl==TRUE){
      HWres=matrix(nrow=nreps, ncol=6)
      colnames(HWres)=c("Mothers.Case","Fathers.Case","Both.Case",
                        "Mothers.Control","Fathers.Control","Both.Control")
    } else {
      #browser()
      HWres=matrix(nrow=nreps, ncol=3)
      colnames(HWres)=c("Mothers","Fathers","Both")
    }

  }

  if (MStest==TRUE){
    MSres=vector(length=nreps)
  }


  # Run the simulation repetitions
  for (i in 1:nreps){
    dat=simulateData(ntrios=ntrios, maf=maf,
                     R=R, S=S,  V=V,
                     mtCoef=mtCoef, mtmodel=mtmodel,
                     includeE=includeE, envint=envint, prE=prE,
                     includePopStrat=includePopStrat, numPop=numPop, Fst=Fst,
                     includeControl=includeControl, prControl=prControl, prE.control=prE.control,
                     prCase.byPop=prCase.byPop, prControl.byPop=prControl.byPop)

    if (HWtest==TRUE){

      if (includeControl==TRUE){
        HWres[i,1:3]=getHWE(subset(dat$dat4R,dat$dat4R$D==1))
        HWres[i,4:6]=getHWE(subset(dat$dat4R,dat$dat4R$D==0))
      } else {
        HWres[i,]=getHWE(dat$dat4R)
      }
    }

    if (MStest==TRUE){
      MSres[i]=MaSlrt(dat$dat4R)
    }


    # For the ith dataset, run the log linear approach if requested
    if (Loglin){
      for (j in 1:length(loglinMT)){ # Loop over mating models selected
        res=runLoglin(loglinMT[j], effects, dat$dat4R, PStest=PStest)
        colvals=paste(loglinMT[j],names(res$effects),sep=".")
        names(res$effects)=colvals
        names(res$pvals)=colvals
        resLoglin[i,colvals]=res$effects[colvals]
        pvalLoglin[i,colvals]=res$pvals[colvals]

        if (PStest==TRUE){
          #browser()
          names(res$effectsPS)=colvals
          names(res$pvalsPS)=colvals
          resLoglinPS[i,colvals]=res$effectsPS[colvals]
          pvalLoglinPS[i,colvals]=res$pvalsPS[colvals]
          pvalPS[i,loglinMT[j]]=res$PS.test
        }
      }
    }

    # For the ith dataset, run haplin if requested
    if (Haplin){
      res=runHaplin(effects,dat$dat4haplin, haplinControls=haplinControls)
      colvals.effects=names(res$effects)
      colvals.pvals=names(res$pvals)
      resHaplin[i,colvals.effects]=res$effects[colvals.effects]
      pvalHaplin[i,colvals.pvals]=res$pvals[colvals.pvals]
    }

    # For the ith dataset, run EMIM if requested
    if (EMIM){
      for (j in 1:length(EMIMMT)){ # Loop over mating models selected
        res=runEMIM(EMIMMT[j], effects, dat$dat4EMIM)
        colvals.res=paste(EMIMMT[j],names(res$effects),sep=".")
        colvals.pvals=paste(EMIMMT[j],names(res$pvals),sep=".")
        names(res$effects)=colvals.res
        names(res$pvals)=colvals.pvals
        resEMIM[i,colvals.res]=res$effects[colvals.res]
        pvalEMIM[i,colvals.pvals]=res$pvals[colvals.pvals]
      }
    }
  }

  return(list(resLoglin=resLoglin, pvalLoglin=pvalLoglin,
              resLoglinPS=resLoglinPS, pvalLoglinPS=pvalLoglinPS,
              resHaplin=resHaplin, pvalHaplin=pvalHaplin,
              resEMIM=resEMIM, pvalEMIM=pvalEMIM,
              pvalHW=HWres, pvalMaS=MSres, pvalPS=pvalPS))
}


## Get ratio of cases/controls by population. Used only under
## population stratification
## Can I remove this??
getRatio=function(alpha, maf, omega, R=c(1,1,1), S=c(1,1,1),V=c(1,1,1),
                  mtcoef=c(1,1,1)){

  numPop=length(alpha)
  K=vector(length=numPop)

  if ((length(maf)!=numPop)||(length(omega)!=numPop)){
    stop("Number of elements of maf, alpha and omega must equal the
         number of subpopulations\n")
  }
  if (!is.list(mtcoef)){
    mtcoef=list(mtcoef)
    mtcoef=rep(mtcoef,numPop)
  }

  browser()
  for (i in 1:numPop){
    genomat=mtmat(maf[i],C=mtcoef[[i]])
    diseasefactor=R[(genomat$C+1)]*S[(genomat$M+1)]
    diseaseprob=alpha[i]*genomat$prMFC*diseasefactor
    K[i]=sum(diseaseprob)
  }

  prPopGivenCase=K*omega/sum(K*omega)

  prPopGivenControl=(1-K)*omega/sum((1-K)*omega)


  return(list(prevalence=K, prPopGivenCase=prPopGivenCase,
              prPopGivenControl=prPopGivenControl))

}
