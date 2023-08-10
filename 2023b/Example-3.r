
## source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023b/Rcc22.R")
install.packages("weibullness")

library(weibullness)


# Full observations
#========================================================================
# Data from Example 8.1 on Page 222 of Leemis (2009). Reliability 2nd ed.
Bearings = c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 51.96,
             54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 
            105.12,105.84,127.92,128.04,173.40)
XX = cbind(Bearings,Bearings)
#========================================================================


# The EM method can find the global maximizer.
weibull.ic (XX, start=c(1,1))

weibull.ic (XX, start=c(5,5))

# =================================================================
# Log-likelihood function for Weibull distribution
# =================================================================
loglike.weibull = function(shape, scale, X){
   TINY = .Machine$double.eps
   ij = dim(X)
   n  = ij[1]
   if ( ij[2] > 2 ) stop(" Warning: The data X should be n x 2 matrix.");
   if ( any(X[,1] > X[,2]) )  stop(" Warning: The data should satisfy a <= b.");

   logL = 0 ;
   for ( i in 1:n ) {
       a = X[i,1]; b = X[i,2];
       if ( abs(a-b) < TINY ) {
         logL = logL + dweibull(a, shape, scale, log=TRUE)
       } else {
         logL = logL + log(pweibull(b,shape,scale)-pweibull(a,shape,scale))
       }
   }
   return(logL)
}
# ----------------------------------------------------------------------
# negative log-likelihood function (needed for optimization)
neg.loglike.weibull = function(para,X) {
   TINY = .Machine$double.eps
   HUGE = .Machine$double.xmax^0.5
   if( para[1] < TINY ) return( HUGE*(1-para[1]) )
   if( para[2] < TINY ) return( HUGE*(1-para[2]) )
   -loglike.weibull(para[1],para[2], X)
}
# ======================================================================


# The Newton-type method fails . 
# Also, the estimates depend on a starting value.
nlm(neg.loglike.weibull, p=c(1,1), X=XX)

nlm(neg.loglike.weibull, p=c(5,5), X=XX)


