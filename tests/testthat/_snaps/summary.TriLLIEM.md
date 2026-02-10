# summary of default (C, M) model works

    Code
      summary(res)
    Output
      
      Call:
      TriLLIEM(dat = example_dat4R)
      
      Coefficients:
        Estimate Std. Error z value Pr(>|z|)
      C  0.05341    0.06970   0.766    0.443
      M -0.10573    0.06794  -1.556    0.120
      AIC: 108.58
      
      Number of Fisher Scoring iterations: 4
      

# summary of control and environment data model works

    Code
      summary(res)
    Output
      
      Call:
      TriLLIEM(mtmodel = "MaS", effects = c("C", "M", "Im", "E:M"), 
          dat = dat, includeE = TRUE, includeD = TRUE)
      
      Coefficients:
           Estimate Std. Error z value Pr(>|z|)
      C    0.065815   0.072799   0.904    0.366
      M    0.001332   0.077984   0.017    0.986
      Im  -0.109649   0.121629  -0.902    0.367
      E:M  0.156938   0.117689   1.334    0.182
      AIC: 408.36
      
      Number of Fisher Scoring iterations: 4
      Number of EM iterations: 10
      

---

    Code
      summary(resStrat)
    Output
      
      Call:
      TriLLIEM(mtmodel = "MaS", effects = c("C", "M", "Im", "E:M"), 
          dat = dat, includeE = TRUE, Estrat = TRUE, includeD = TRUE)
      
      Coefficients:
            Estimate Std. Error z value Pr(>|z|)
      C     0.068886   0.085289   0.808    0.419
      M     0.001365   0.082564   0.017    0.987
      Im   -0.112869   0.140278  -0.805    0.421
      E:C  -0.010714   0.164325  -0.065    0.948
      E:M   0.156646   0.155338   1.008    0.313
      E:Im  0.011527   0.281630   0.041    0.967
      AIC: 412.36
      
      Number of Fisher Scoring iterations: 4
      Number of EM iterations: 11
      

# population stratification model works

    Code
      summary(resHWE)
    Output
      
      Call:
      TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im", "E:Im"), 
          dat = dat, includeE = TRUE, includeD = TRUE)
      
      Coefficients:
           Estimate Std. Error z value Pr(>|z|)
      C    -0.09664    0.08260  -1.170    0.242
      M    -0.05193    0.08147  -0.637    0.524
      Im    0.16388    0.19056   0.860    0.390
      E:Im  0.18043    0.15735   1.147    0.252
      AIC: 394.4
      
      Number of Fisher Scoring iterations: 4
      Number of EM iterations: 11
      

---

    Code
      summary(resMS)
    Output
      
      Call:
      TriLLIEM(mtmodel = "MS", effects = c("C", "M", "Im", "E:Im"), 
          dat = dat, includeE = TRUE, includeD = TRUE)
      
      Coefficients:
           Estimate Std. Error z value Pr(>|z|)
      C    -0.09798    0.08283  -1.183    0.237
      M    -0.05463    0.08161  -0.669    0.503
      Im    0.17466    0.19471   0.897    0.370
      E:Im  0.17721    0.15853   1.118    0.264
      AIC: 379.55
      
      Number of Fisher Scoring iterations: 4
      Number of EM iterations: 11
      

---

    Code
      summary(resMaS)
    Output
      
      Call:
      TriLLIEM(mtmodel = "MaS", effects = c("C", "M", "Im", "E:Im"), 
          dat = dat, includeE = TRUE, includeD = TRUE)
      
      Coefficients:
           Estimate Std. Error z value Pr(>|z|)
      C    -0.09369    0.08729  -1.073    0.283
      M    -0.05567    0.08447  -0.659    0.510
      Im    0.15117    0.20115   0.752    0.452
      E:Im  0.20170    0.16764   1.203    0.229
      AIC: 387.62
      
      Number of Fisher Scoring iterations: 4
      Number of EM iterations: 11
      

