#======================================================================
# Example-E: Censored Weibull data (interval-censored): QEM and MCEM comparison
#======================================================================

#----------------------------------------------------------------------
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2018/Rs2.R")
#----------------------------------------------------------------------

## Data from Nelson (1982). "Applied Life Data Analysis"
X1 = c( 0,    6.12, 19.92, 29.64, 35.40, 39.72, 45.24, 52.32, 63.48 )
X2 = c(6.12, 19.92, 29.64, 35.40, 39.72, 45.24, 52.32, 63.48, Inf )
f  = c(   5,    16,    12,    18,    18,     2,     6,    17,  73 )

#-----------------------------------------------
# Data
Xinterval = cbind( rep(X1,f), rep(X2,f) )

# QEM 
para.wei.QEM    = ic.weibull.QEM (Xinterval, beta0=1, lam0=1, K=10)

# MCEM 
para.wei.MCEMa  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=10)
para.wei.MCEMb  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=100)
para.wei.MCEMc  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=1000)

# QEM 
print(para.wei.QEM)    

# MCEM 
print(para.wei.MCEMa)
print(para.wei.MCEMb)
print(para.wei.MCEMc)

####################################





