
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2024c/fna53.R")

#============================================================================
# Data Generation 
X = c(0.219, 0.386, 0.463, 0.489, 0.495, 0.526, 0.528, 0.528, 0.552, 0.554, 
      0.561, 0.585, 0.592, 0.672, 0.698, 0.796, 0.824, 0.945, 0.950, 0.971, 
      0.974, 1.042, 1.055, 1.107, 1.124, 1.234, 1.317, 1.361, 1.397, 1.504, 
      1.569, 1.673, 1.686, 1.761, 1.762, 1.798, 2.056, 2.093, 2.097, 2.144,
      2.168, 2.200, 2.207, 2.282, 2.471, 2.502, 2.918, 3,     3,     3)

d = list(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
         1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 
         1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 3L, c(2L,3L),  
         0L, 0L, 0L)
#============================================================================

##-----------------------------------------------------------------------
## Weibull Model
##-----------------------------------------------------------------------
  cat("\n==================\n")
  cat(" Weibull Model")
  cat("\n==================\n")
  weibull.cm.EM(X,d, eps=1.0E-5)

##-----------------------------------------------------------------------
## Inverse Weibull Model 
##-----------------------------------------------------------------------
  cat("\n==================\n")
  cat(" Inverse Weibull Model")
  cat("\n==================\n")
  invweibull.cr.EM(X,d, beta0=1:3, theta0=1:3, eps=1.0E-7)

