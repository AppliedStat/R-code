#============================================================================
source ("https://raw.githubusercontent.com/AppliedStat/R-code/master/2024b/fnmw41.R")
# If the above does not work, then just cut and past the R codes from 
#     https://github.com/AppliedStat/R-code/blob/master/2024b/fnmw41.R


#============================================================================
# Simulated Data
X = c(
0.4728, 0.3510, 0.2111, 0.0737, 0.6720, 0.1386, 1.2220, 0.0369, 0.1132, 0.6970,
1.5903, 0.1306, 1.0736, 0.7670, 0.3262, 0.0289, 0.0671, 0.1144, 2.0000, 0.4998,
1.6393, 0.0576, 0.2278, 0.3829, 0.4344, 0.1592, 1.2487, 0.5002, 0.1835, 0.9386,
0.4129, 0.0262, 1.3644, 1.1828, 0.0309, 0.0402, 0.0913, 0.2573, 0.0807, 1.0978,
0.4548, 0.1692, 0.0167, 0.1927, 0.0423, 0.6103, 0.2839, 0.2174, 0.0928, 0.0427)
d = list(
1L, 2L, 2L, 3L, 2L, 3L, 3L, 3L, 3L, c(1,2),
3L, 2L, 1L, 1L, 2L, 3L, 3L, 1L, 0L, c(1,2,3),
2L, 2L, 2L, 2L, 3L, 3L, 1L, 1L, 2L, 2L,
1L, 3L, 2L, 1L, 3L, 3L, 3L, 2L, 2L, 1L,
1L, 1L, 3L, 3L, 3L, 2L, 2L, 3L, 2L, c(3,2) )

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

