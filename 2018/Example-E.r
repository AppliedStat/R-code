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
para.wei.QEM    = ic.weibull.QEM (Xinterval, beta0=1, lam0=1, K=100)

# MCEM 
para.wei.MCEMa  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=100)
para.wei.MCEMb  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=200)
para.wei.MCEMc  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=300)
para.wei.MCEMd  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=400)
para.wei.MCEMe  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=500)
para.wei.MCEMf  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=600)
para.wei.MCEMg  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=700)
para.wei.MCEMh  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=800)
para.wei.MCEMi  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=900)
para.wei.MCEMj  = ic.weibull.MCEM(Xinterval, beta0=1, lam0=1, K=1000)

# QEM 
print(para.wei.QEM)    

# MCEM 
print(para.wei.MCEMa)
print(para.wei.MCEMb)
print(para.wei.MCEMc)

#....

print(para.wei.MCEMj)


###########################################################################################
KK = seq(100, 1000, by=100)

BETA.MCEM = c(para.wei.MCEMa$beta, para.wei.MCEMb$beta, para.wei.MCEMc$beta, para.wei.MCEMd$beta,
              para.wei.MCEMe$beta, para.wei.MCEMf$beta, para.wei.MCEMg$beta, para.wei.MCEMh$beta,
              para.wei.MCEMi$beta, para.wei.MCEMj$beta)
plot(KK, BETA.MCEM, type="b")
abline(h=para.wei.QEM$beta, col="red")


LAM.MCEM = c(para.wei.MCEMa$lam, para.wei.MCEMb$lam, para.wei.MCEMc$lam, para.wei.MCEMd$lam,
              para.wei.MCEMe$lam, para.wei.MCEMf$lam, para.wei.MCEMg$lam, para.wei.MCEMh$lam,
              para.wei.MCEMi$lam, para.wei.MCEMj$lam)
plot(KK, LAM.MCEM, type="b")
abline(h=para.wei.QEM$lam, col="red")










