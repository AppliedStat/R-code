#============================================================================
source ("https://raw.githubusercontent.com/AppliedStat/R-code/master/2024b/fnmw41.R")
# If the above does not work, then just cut and past the R codes from 
#     https://github.com/AppliedStat/R-code/blob/master/2024b/fnmw41.R


#============================================================================
# Simulated Data
X = c(0.9189, 0.0610, 0.0783, 0.0686, 1.0000, 0.3857, 0.2394, 0.3201, 0.0911, 0.5572, 
      0.4759, 0.2178, 0.0684, 0.1477, 0.5894, 0.1543, 0.3881, 0.7628, 0.1326, 0.0824,
      0.0950, 0.0726, 0.4851, 0.3980, 1.0000, 0.3474, 0.0281, 1.0000, 0.1073, 0.1691, 
      1.0000, 0.2348, 0.0916, 0.1094, 0.0400, 1.0000, 0.8260, 0.1839, 0.9372, 0.0261,
      1.0000, 0.5459, 0.3801, 0.0573, 0.1087, 0.1768, 0.5211, 0.9155, 0.1643, 0.4421)
d = list(
1L, c(3,1), 3L, 3L, 0L, 3L, 2L, 2L, 1L, 1L, 
1L, 3L, 2L, 2L, 2L, 2L, 2L, 1L, 3L, 2L, 
2L, 3L, 3L, 2L, 0L, 1L, 2L, 0L, 3L, 1L, 
0L, 2L, 3L, 3L, 3L, 0L, 2L, 2L, 1L, 2L, 
0L, 1L, 2L, 3L, 2L, 1L, 1L, 2L, 3L, c(3,2,1))
Strength = X; Mode=d
#============================================================================

#============================================================================
library(survival)

dataALL = sort( jitter(Strength, factor=1E-1) )

idx     = order( Strength )
ModeALL = as.list( rep(NA, length(idx)) )
is.censored = numeric( length(idx) )

for ( i in 1:length(idx) ) {
    ModeALL[[i]] = Mode[[idx[i]]]
    is.censored [i] = !( all(ModeALL[[i]]==0) )   ## 0=censored, 1=not censored
}

fit = survfit ( Surv(dataALL, (is.censored) )~1 )  ## ~1 is needed for new version
Fall = 1 - fit$surv

##-----------------------------------------------------------------------
## BS Model 
##-----------------------------------------------------------------------
 SBS <- function(q, alpha, beta)  {   # Survival function of BS
    y1 = sqrt(q/beta[1])
    y2 = sqrt(q/beta[2])
    y3 = sqrt(q/beta[3])
    pnorm(y1-1/y1,sd=alpha[1],lower.tail=FALSE) * 
    pnorm(y2-1/y2,sd=alpha[2],lower.tail=FALSE) * 
    pnorm(y3-1/y3,sd=alpha[3],lower.tail=FALSE) 
 }
 para.BS = BS.cm.QEM(X,d, eps=1.0E-3, maxits=1000, K=1000)
 FBS = 1-SBS(sort(X), alpha=para.BS$alpha, beta=para.BS$beta)
 MSE= mean( (FBS-Fall)^2 )
##=========================================================================



##=========================================================================
cat("\n\n =====================\n")
cat(" MSE =", MSE)

cat("\n\n =====================\n")
cat(" BS parameter\n")
print(para.BS)

