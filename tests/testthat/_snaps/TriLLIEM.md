# default (C, M) model works

    Code
      res
    Output
      
      Call:  TriLLIEM(dat = example_dat4R)
      
      Coefficients:
             C         M  
       0.05341  -0.10573  
      
      Degrees of Freedom: 15 Total (i.e. Null);  7 Residual
      Null Deviance:	    6940 
      Residual Deviance: 6.985 	AIC: 108.6

# control and environment data model works

    Code
      res
    Output
      
      Call:  TriLLIEM::TriLLIEM(mtmodel = "MaS", effects = c("C", "M", "Im", 
          "E:M"), dat = dat, includeE = TRUE, includeD = TRUE)
      
      Coefficients:
              C          M         Im        E:M  
       0.065815   0.001332  -0.109649   0.156938  
      
      Degrees of Freedom: 60 Total (i.e. Null);  36 Residual
      Null Deviance:	    25160 
      Residual Deviance: 35.4 	AIC: 408.4

---

    Code
      resStrat
    Output
      
      Call:  TriLLIEM(mtmodel = "MaS", effects = c("C", "M", "Im", "E:M"), 
          dat = dat, includeE = TRUE, Estrat = TRUE, includeD = TRUE)
      
      Coefficients:
              C          M         Im        E:C        E:M       E:Im  
       0.068886   0.001365  -0.112869  -0.010714   0.156646   0.011527  
      
      Degrees of Freedom: 60 Total (i.e. Null);  34 Residual
      Null Deviance:	    25160 
      Residual Deviance: 35.4 	AIC: 412.4

