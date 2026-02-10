# simple example_dat4R anova works

    Code
      anova(model_1, model_2)
    Output
      Analysis of Deviance Table
      
      Factors included in model: C, M
      
      Response: count 
      
      Terms added sequentially (first to last)
      
      
        Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
      1         6     6.9853                       
      2         5     1.0541  1   5.9313  0.01487 *
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Eanova works

    Code
      anova(res1, res2)
    Output
      Analysis of Deviance Table
      
      Factors included in model: C, M, Im, E:Im
      
      Response: count 
      
      Terms added sequentially (first to last)
      
      
        Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
      1        46      85.98                          
      2        50     653.97 -4  -567.99 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

