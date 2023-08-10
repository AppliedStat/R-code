
## source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023b/Rcc22.R")
install.packages("weibullness")

library(weibullness)

#========================================================================
# Data from Table 3.10 (Radiotherapy and Chemotherapy) of Lawless (2003). 
# Statistical Models and Methods for Lifetime Data. Wiley, New York.
# Lower ends in [a,b]
X1 = c(
 8, 0, 24, 17, 17, 24, 16, 13, 11, 16,
18, 17, 32, 23, 44, 10, 0, 5, 12, 11,
33, 31, 13, 19, 34, 13, 16, 35, 15, 11,
22, 48, 30, 13, 10, 8, 4, 11, 14, 4,
34, 30, 18, 16, 35, 21, 11 )
# Upper ends in [a,b]
X2 = c(
 12, 22, 31, 27, 23, 30, 24, Inf, 13, 20,
 25, 26, Inf, Inf, 48, 35, 5, 8, 20, Inf,
 40, Inf, 39, 32, Inf, Inf, 24, Inf, 22, 17,
 32, Inf, 34, Inf, 17, 21, 9, Inf, 19, 8,
 Inf, 36, 24, 60, 39, Inf, 20)
XX = cbind(X1,X2)
#========================================================================

# The EM method can find the global maximizer.
weibull.ic (XX, start=c(1,1))

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


# The Newton-type method fails. 
nlm(neg.loglike.weibull, p=c(1,1), X=XX)

# The Newton-type method is successful with a certain starting value.
nlm(neg.loglike.weibull, p=c(2,30), X=XX)


