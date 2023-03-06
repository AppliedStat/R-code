#============================================================================
source ("https://raw.githubusercontent.com/AppliedStat/R-code/master/2024b/fnmw41.R")
# If the above does not work, then just cut and past the R codes from 
#     https://github.com/AppliedStat/R-code/blob/master/2024b/fnmw41.R


#============================================================================
# Simulated Data
Strength = 
c(0.754, 0.279, 0.261, 0.902, 0.323, 1.16,  1.264, 0.632, 2,     0.174,
  0.135, 0.182, 0.205, 0.019, 0.101, 0.419, 0.106, 0.22,  1.109, 0.312 )
Mode = 
list(1L, 2L, 2L, 1L, 2L, 3L, 2L, 3L, 0L, c(3,1), 
     2L, 2L, 2L, 2L, 2L, 2L, 3L, 2L, 2L, c(2,1,3))

#============================================================================
X = Strength; d = Mode
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

logT0 = log(fit$time[ is.censored==0 ])
logT1 = log(fit$time[ is.censored> 0 ])

loglog0 = log( -log(fit$surv[ is.censored==0 ]) )
loglog1 = log( -log(fit$surv[ is.censored >0 ]) )

 xlim  = range( logT0, logT1, na.rm=TRUE )
 ylim1 = min( loglog0, loglog1, na.rm=TRUE )
 ylim2 = max( loglog0, loglog1, na.rm=TRUE )

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
 FBS = 1-SBS(X, alpha=para.BS$alpha, beta=para.BS$beta)
 MSE = mean( (FBS-Fall)^2 )
##=========================================================================



##=========================================================================
cat("\n\n ====================================\n")
cat(" MSE =", MSE)

cat("\n\n ====================================\n")
cat(" BS parameter\n")
print(para.BS)

