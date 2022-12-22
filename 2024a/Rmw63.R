#------------------------------
# Generalized MSE 
#------------------------------
MSE.gen = function(x,y, mux, muy) {
   N = length(x)
   a11 = sum( (x-mux)^2 )
   a22 = sum( (y-muy)^2 )
   a12 = sum( (x-mux)*(y-muy) )
   S = 1/N * matrix( c(a11,a12,a12,a22), nrow=2)
   det(S)
}

#===================================================================
# Weighted median
#===================================================================
#--------------
# new version on Oct. 31. 2019
# NOTE: a = 10L; b = exp(log(10)) ;  c(a,b); a==b
#  tiny=.Machine$double.eps; abs(a-b)<tiny  # good idea, but does not work
#  tiny = sqrt(.Machine$double.eps); abs(a-b)<tiny # It works!
#--------------
weighted.median <- function(x, w, interpolation=0.5) { 
  # Preparation
  if (missing(w)) w = rep(1L,length(x))
  if (length(w) != length(x)) stop("'x' and 'w' must have the same length")

  x = as.double(as.vector(x))
  w = as.double(as.vector(w))
  ok= complete.cases(x,w); x=x[ok]; w=w[ok]

  stopifnot(all(w >= 0))
  if(all(w <= 0)) stop("All weights are zero", call.=FALSE)

  orderx = order(x)
  x = x[orderx]
  w = w[orderx] / sum(w)
  Fn = cumsum(w)
  tiny = sqrt(.Machine$double.eps)

  # Main part
  if ( all( abs(Fn-0.5)>tiny ) ) {  # any values of Fn is not 1/2.
      k = sum( Fn < 0.5 )
      return( x[k+1] )
  } else {
    k = which.min ( signif(abs(Fn-0.5),digits=12) ) # Find k with Fn=0.5
    if (w[k+1] < tiny) {   # check if w[k+1] == 0 
        return( x[k+1] )
    } else {
      return( (1-interpolation)*x[k] + interpolation*x[k+1] )
    }
  }
}
#------------------------------
# This method seems simple and best!
Repeated.W.Median <- function (x,y, weight=FALSE, power=1) {
  n = length(x)
  dx = outer(x,x,"-"); dy = outer(y,y,"-");
  slope1 = dy/dx 
  weight1 = abs(dx)^power 
  slope2  = matrix( t(slope1)[lower.tri(slope1)  |upper.tri(slope1) ], byrow=TRUE, nrow=n)
  weight2 = matrix( t(weight1)[lower.tri(weight1)|upper.tri(weight1)], byrow=TRUE, nrow=n)
  tmp = numeric(n)
  if (weight) { 
     for ( j in 1L:n ) tmp[j] = weighted.median(slope2[j,], w=weight2[j,] )
  } else {
     for ( j in 1L:n ) tmp[j] = median(slope2[j,])
  }
  slope = median(tmp)
  intercept = median(y-slope*x)
  return(c(intercept,slope))
}
#------------------------------
# Different method for estimating intercept.
# More complex, and worse 
Repeated.W.Median1 <- function (x,y, weight=FALSE, power=1) {
  n = length(x)
  dx = outer(x,x,"-"); dy = outer(y,y,"-");
  slope1 = dy/dx
  xyslope1 = (outer(x,y)-outer(y,x))/dx
  weight1 = abs(dx)^power
  slope2  = matrix(   t(slope1)[lower.tri(slope1)  |upper.tri(slope1) ],  byrow=TRUE, nrow=n)
  xyslope2= matrix( t(xyslope1)[lower.tri(xyslope1)|upper.tri(xyslope1)], byrow=TRUE, nrow=n)
  weight2 = matrix(  t(weight1)[lower.tri(weight1) |upper.tri(weight1)],  byrow=TRUE, nrow=n)
  tmp0 = tmp1 = numeric(n)
  if (weight) {
     for ( j in 1L:n ) tmp0[j] = weighted.median(xyslope2[j,], w=weight2[j,] )
     for ( j in 1L:n ) tmp1[j] = weighted.median(slope2[j,],   w=weight2[j,] )
     intercept = median(tmp0)
  } else {
     for ( j in 1L:n ) tmp1[j] = median(slope2[j,])
     intercept = median(y-slope*x)
  }
  slope = median(tmp1)
  return(c(intercept,slope))
}
#------------------------------
# Breakdown point
RM.breakdown.point <- function(x, power=1) {
   n = length(x)
   dx = outer(x,x,"-");
   weight1 = abs(dx)^power 
   weight2 = matrix( t(weight1)[lower.tri(weight1)|upper.tri(weight1)], byrow=TRUE, nrow=n)
   WU = t( apply(weight2, 1, sort) ) 
   WL = t( apply(weight2, 1, sort, decreasing=TRUE) )
   BU = numeric(n)  # Breakdown (Upper)
   BL = numeric(n)  # Breakdown (Lower)
   for ( i in seq_len(n) ) {   
       BU[i] = sum( cumsum(WU[i,])/sum(WU[i,]) < 0.5 )
       BL[i] = sum( cumsum(WL[i,])/sum(WL[i,]) < 0.5 )
   }
   return(list(lower=min(BL)/(n-1), upper=max(BU)/(n-1)))
}
#------------------------------

# MLE of Weibull
weibull.MLE <- function (x, tol=.Machine$double.eps^0.25, maxiter=1000, trace=0) {
  # Setup for interval 
    meanlog = mean(log(x))
    lower = 1/(log(max(x)) - meanlog)
    upper = sum((x^lower) * log(x))/sum(x^lower) - meanlog
    interval = c(lower, 1/upper)
  # EE equation 
    EEweibull = function(alpha, x) {
        xalpha = x^alpha
        sum(log(x) * (xalpha))/sum(xalpha) - 1/alpha - mean(log(x))
    }
    tmp = uniroot(EEweibull, interval = interval, x = x, tol = tol,
        maxiter = maxiter, trace = trace)
    alpha = tmp$root
    beta = mean(x^alpha)^(1/alpha)
    return( c(alpha,beta) )
} # END of MLE of Weibull

#=================================================================
# Biranbaum-Saunders (See research/Lio/L01/R/fnL01.R)
# MLE of Biranbaum-Saunders
BS.MLE <- function(x) {
   r = 1/ mean(1/x)
   s =  mean(x)
   # lexical scoping of R (so, this won't work in S/Splus). 
   EE = function(beta) { 
      k = 1 / mean( 1/(beta+x) )
      return (beta^2 - beta*(2*r+k) + r*(s+k))
   }
   tmp = uniroot(EE, interval=c(r,s))
   beta = tmp$root
   alpha = sqrt ( s/beta + beta/r - 2 )
   ##list(alpha=alpha, beta=beta)
   return( c(alpha,beta) )
}
# CDF of Birnbaum-Saunders distribution 
pBS <- function(q, alpha=1, beta=1)  {
   y = sqrt(q/beta)
   pnorm(y-1/y, sd=alpha)
}
# pdf of Birnbaum-Saunders distribution 
dBS <- function(x, alpha=1, beta=1)  {
   y = sqrt(x/beta)
   (0.5/x) * (y+1/y) * dnorm(y-1/y, sd=alpha)
}
# random number of Birnbaum-Saunders distribution 
rBS <- function(n, alpha=1, beta=1)  {
   x = rnorm(n, sd=alpha/2)
   x2 = x*x
   beta*(1+2*x2+2*x*sqrt(1+x2))
}
# quantile of Birnbaum-Saunders distribution 
qBS <- function(p, alpha=1, beta=1)  {
   z = qnorm(p, sd=alpha/2)
   az = alpha*z
   .25*beta*(az + sqrt(az*az+4))^2
} ## NOTE qBS is wrong. The below is correct. (I found 2012-09-09)
qBS2 <-
function(p, alpha=1, beta=1)  {
   az = alpha*qnorm(p)
   .25*beta*(az + sqrt(az*az+4))^2
}

#=================================================================
rLL <- function(x, alpha=1, beta=1) {
    alpha * (1/(1/runif(n)-1))^(1/beta)
}
