

# Examples using simulate_data_functions_V3.R and the functions
# in this file 
source("simulate_data_functions_V3.R")


## Example 1 - no E and no control trios
dat=simulateData()

# Returns the p-value of a MaS test
MaSlrt(dat$dat4R)

# Run the log-linear analysis
runLoglin(mtmodel="HW", effects=c("C","M"), dat=dat$dat4R, PStest=FALSE)
runLoglin(mtmodel="MS", effects=c("C","M"), dat=dat$dat4R, PStest=FALSE)

## Example 2 - E and no control trios
dat=simulateData(includeE=TRUE,prE=0.3)

# Returns the p-value of a MaS test
MaSlrt(dat$dat4R)

# Run the log-linear analysis
runLoglin(mtmodel="HW", effects=c("C","M","E:M"), dat=dat$dat4R, PStest=FALSE)
runLoglin(mtmodel="MS", effects=c("C","M","E:M"), dat=dat$dat4R, PStest=FALSE)


## Example 3 - E and control trios
dat=simulateData(ntrios=2000,includeE=TRUE, prE=0.3, includeControl = TRUE, prControl = 0.5)


# Run the log-linear analysis
runLoglin(mtmodel="HW", effects=c("C","M","E:M"), dat=dat$dat4R, PStest=TRUE)
runLoglin(mtmodel="MS", effects=c("C","M","E:M"), dat=dat$dat4R, PStest=TRUE)
runLoglin(mtmodel="MaS", effects=c("C","M","E:M"), dat=dat$dat4R, PStest=TRUE)


## ---------------------------------------------------------------- ##
## Function to include in R package
## - copied from run_analysis_functions_V3.R
## ---------------------------------------------------------------- ##


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
  
  return(list(effects=resVec, pvals=pvalVec, PS.test=test.res$`Pr(>Chi)`[2], 
              effectsPS=resVecPS, pvalsPS=pvalVecPS))
}

