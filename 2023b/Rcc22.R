# =================================================================
# MLE of Weibull (with full observations)
# kappa=shape parameter  theta=scale parameter
#------------------------------------------------------------------
weibull.MLE <- function(x, interval) {
  if ( any (x <= 0) ) stop("The data should be positive")
  if (missing(interval)) {
     meanlog = mean(log(x))
     lower = 1 / ( log(max(x)) - meanlog )
     upper = sum( (x^lower)*log(x) ) / sum( x^lower ) - meanlog
     interval = c(lower,1/upper)
  }
  EE = function(kappa,x) {
     xkappa = x^kappa
     sum(log(x)*(xkappa)) / sum(xkappa) - 1/kappa - mean(log(x))
  }
  tmp = uniroot(EE, interval=interval, x=x)
  kappa = tmp$root ;  theta = (mean(x^kappa))^(1/kappa);
  names(kappa)="shape"; names(theta)="scale"
  list ( shape=kappa, scale=theta )
}
# =================================================================

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


# =================================================================
# Incomplete Upper Gamma Function 
#------------------------------------------------------------------
gamma.inc.old = function(a,b) { pgamma(b,a,lower.tail=FALSE)*gamma(a) }
lgamma.inc = function(a,b) {
    if( a < 0 ) stop("a must be non-negative")
    if( a == 0 ) a = .Machine$double.xmin
    return( pgamma(b,a,lower.tail=FALSE,log.p=TRUE) + lgamma(a) )
}
gamma.inc = function(a,b) { exp(lgamma.inc(a,b)) }
# =================================================================

# ===================================================================
# Version   : 1.0,   Dec. 25, 2016
#             1.1,   Feb. 14, 2020
#             1.2,   Jan.  4, 2023
# NOTE: Euler-Mascheroni constant gam = -digamma(1)
# ===================================================================
U.i.s = function(kappa.s, theta.s, ai, bi) { 
   TINY = .Machine$double.neg.eps

   if ( ai > bi ) return(0.0)
   if ( abs(bi-ai) < TINY ) return(log(ai))

   t.ai = (ai/theta.s)^kappa.s 
   exp.t.ai = exp(-t.ai) 

   t.bi = (bi/theta.s)^kappa.s 
   exp.t.bi = exp(-t.bi) 

   D.i.s = exp.t.ai - exp.t.bi

   if (ai>0) {
      tmpai = log(ai)*exp.t.ai + gamma.inc(0,t.ai)/kappa.s
   } else {
      tmpai = log(theta.s) + digamma(1)/kappa.s
   }

   BIG = .Machine$double.xmax^0.2
   if (bi >= BIG) {   ## .Machine$double.xmax
      tmpbi = 0
   } else {
      tmpbi = log(bi)*exp.t.bi + gamma.inc(0,t.bi)/kappa.s
   }

   return( (tmpai-tmpbi)/D.i.s )
}
#
U.i.s.0.inf = function(kappa.s, theta.s) { (log(theta.s) + digamma(1)/kappa.s) }
#
# U.i.s (2, 3, 1, Inf)
# U.i.s (2, 3, 0, Inf);  U.i.s.0.inf(2, 3)
# U.i.s (2, 3, 3, 3);   U.i.s (2, 3, 3, 3.0001)
#
# aa = seq(0,1, length=101) ; U  = numeric( length(aa) )
# for ( i in seq_along(aa) ) U[i] = U.i.s(1, 2, aa[i], 1)
# plot (aa, U, type="l")
# 
# bb = seq(2, 6 , l=101); U  = numeric( length(bb) )
# for ( i in seq_along(bb) ) U[i] = U.i.s(1, 2, 3, bb[i])
# plot (bb, U, type="l"); abline(h=log(3), col="green")
# 

# ===================================================================
# Version   : 1.0,   Dec. 25, 2016
#             1.1,   Feb. 14, 2020
#             1.2,   Jan.  4, 2023
# NOTE: Euler-Mascheroni constant gam = -digamma(1)
# ===================================================================
V.i.s = function(kappa, kappa.s, theta.s, ai, bi) { 
   TINY = .Machine$double.neg.eps
   if ( ai > bi ) return(0.0)
   if ( abs(bi-ai) < TINY ) return(ai^kappa)

   t.ai = (ai/theta.s)^kappa.s
   t.bi = (bi/theta.s)^kappa.s

   D.i.s = exp(-t.ai) - exp(-t.bi)

   tmpai = gamma.inc( (kappa+kappa.s)/kappa.s, (ai/theta.s)^kappa.s )
   tmpbi = gamma.inc( (kappa+kappa.s)/kappa.s, (bi/theta.s)^kappa.s )

   theta.s^kappa * (tmpai-tmpbi)/D.i.s
}
# V.i.s (1, 2, 3, 1, Inf)
# V.i.s (1, 2, 3, 9, 9);   V.i.s (1, 2, 3, 8.9999, 9)
# V.i.s (1, 2, 3, 1, 1);   V.i.s (1, 2, 3, 0.9999, 1)

##################################################################
## V.i.s = function(kappa, kappa.s, theta.s, ai, bi) 
## U.i.s = function(kappa.s, theta.s, ai, bi) 
weibull.ic.EM <-
function(X, start=c(1,1), maxits=10000, eps=1E-5){
   kappa = start[1]
   theta = start[2]
   ij = dim(X)
   n  = ij[1]
   if ( ij[2] > 2 ) stop(" Warning: The data X should be n x 2 matrix.")
   if ( any(X[,1] > X[,2]) )  stop(" Warning: The data should satisfy a <= b.");
   iter = 0
   converged = FALSE

   # Start the EM
   colnames(X) = NULL; rownames(X) = NULL

   ai = X[,1];  bi = X[,2]

   fn = function(newkappa) {
        sumU.i.s = 0
        sumV.i.s = 0
        for ( i in seq_len(n) ) sumU.i.s = sumU.i.s + U.i.s(kappa, theta, ai[i], bi[i])
        for ( i in seq_len(n) ) sumV.i.s = sumV.i.s + V.i.s(newkappa, kappa, theta, ai[i], bi[i])
        -1 * ( n*log(newkappa) + newkappa*sumU.i.s - n*log(sumV.i.s) )
   }

   while( (iter<maxits)&(!converged) ) {
      ## OPT = nlm(fn, kappa); newkappa = OPT$estimate
      OPT = nlminb(kappa,fn)
      newkappa = OPT$par
    
      meanV.i.s = 0
      for ( i in seq_len(n) ) meanV.i.s = meanV.i.s + V.i.s(newkappa, kappa, theta, ai[i], bi[i])/n
      newtheta = meanV.i.s^(1/newkappa)

      # assess convergence
      converged = (abs(newkappa-kappa)<eps*abs(newkappa)) & (abs(newtheta -theta)<eps*abs(newtheta))
      iter = iter+1
      kappa   = newkappa
      theta   = newtheta
   }
   list( shape=newkappa, scale=newtheta, iter=iter, conv=converged )
}
#------------------------------------------------------------------


