#======================================================================
# Example-D: Censored Weibull data (interval-censored)
#======================================================================

#----------------------------------------------------------------------
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2018/Rs2.R")
#----------------------------------------------------------------------

## Data from Nelson (1982). "Applied Life Data Analysis"
X1 = c( 0, 6.12, 19.92, 29.64, 35.40, 39.72, 45.24, 52.32, 63.48 )
X2 = c(    6.12, 19.92, 29.64, 35.40, 39.72, 45.24, 52.32, 63.48, Inf )
f  = c( 5, 16, 12, 18, 18, 2, 6, 17, 73 )

#-----------------------------------------------
Xinterval = cbind( rep(X1,f), rep(X2,f) )
para.exp  = ic.exp.EM(Xinterval, start=1, eps=1E-5)
para.wei  = ic.weibull.QEM(Xinterval, beta0=1, lam0=1, eps=1E-5)

cat("\n\n***** Exponential Model *****\n")
print(para.exp)


cat("\n\n***** Weibull Model ***** \n")
print(para.wei)


