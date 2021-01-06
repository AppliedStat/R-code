# RIGHT-CENSORING 
#=====================================================================
#  EM / MCEM / QEM algorithms for right censored normal sample
#  Written by C. Park
#---------------------------------------------------------------------
#  Required input arguments:
#       y = complete univariate data
#       R = right-censored failure time
#   start = starting values for parameters.
#           This should be a vector: c(mu,sigma)
#
#  Optional arguments
#   maxits = maximum number of iterations
#      eps = convergence criterion
#
#  Output value:
#  This function returns a list with the following components
#     para = final parameter estimates
#     iter = how many iterations were performed
#     converged = logical value indicateing whether it converged or not
#---------------------------------------------------------------------
# Example: 
#           y <- c(1.613, 1.644, 1.663, 1.732, 1.740, 1.763, 1.778)
#           R <- c(1.778, 1.778, 1.778)
#           rc.norm.em(y,R, start=c(mean(y),sqrt(var(y))) )
#           rc.norm.em(y,R, start=c(0,1), eps=0.005 )
#=======================================================================
#
#
#
#=======================================================================
# Normal Distribution 
#=======================================================================
# EM 
rc.norm.EM <-
function(y,R,start,maxits=500,eps=.0001){
   mu   <- start[1]
   sigma<- start[2]
   m    <- length(y)
   n    <- m + length(R)
   T1   <- sum(y)
   T2   <- sum(y^2)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out <- NULL
   while((iter<maxits)&(!converged)){
      r   <- (R-mu)/sigma
      sdr <- sigma*dnorm(r)
      pr  <- 1-pnorm(r)
      S1  <- (n-m)*mu + sum( sdr/pr )
      S2  <- (n-m)*(mu^2+sigma^2) + sum( (mu+R)*sdr/pr )
      newmu   <- (T1+S1)/n
      sigma2  <- (T2+S2)/n - newmu^2;
      newsigma<- sqrt(sigma2)
      # assess convergence
      converged <- (abs(newmu-mu)<eps*abs(mu)) &
                   (abs(newsigma-sigma)<eps*abs(sigma))
      iter  <- iter+1
      mu    <- newmu
      sigma <- newsigma
      ## cat("iter =", iter, ": ", tmp, "\n")
      out   <- rbind(out, c(mu,sigma))
   }
   cat("Done.\n")
   theta  <- c(newmu,newsigma)
   list(para=theta,iter=iter,conv=converged,out=out)
}
#---------------------------------------------------------------------
# MCEM 
rc.norm.MCEM <-
function(y,R,start,maxits=500,eps=.0001,K=10000) {
   mu   <- start[1]
   sigma<- start[2]
   m    <- length(y)
   n    <- m + length(R)
   T1   <- sum(y)
   T2   <- sum(y^2)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out <- NULL
   while((iter<maxits)&(!converged)){
      z   <- rtnorm( K*(n-m), mean=mu, sd=sigma, R=R)
      U1  <- sum(z)
      U2  <- sum(z^2)

      newmu   <- (T1 + U1/K) / n
      sigma2  <- (T2+U2/K)/n - newmu^2;
      newsigma<- sqrt(sigma2)
      # assess convergence
      converged <- (abs(newmu-mu)<eps*abs(mu)) &
                   (abs(newsigma-sigma)<eps*abs(sigma))
      iter  <- iter+1
      mu    <- newmu
      sigma <- newsigma
      ## out   <- formatC(c(mu,sigma),wid=8,dig=4,format="f")
      ## cat("iter =", iter, ": ", out, "\n")
      out   <- rbind(out, c(mu,sigma))
   }
   cat("Done.\n")
   theta  <- c(newmu,newsigma)
   list(para=theta,iter=iter,conv=converged,out=out)
}
#--------------------------------------------------------------------
# QEM 
rc.norm.QEM <- 
function(y,R,start,maxits=500,eps=.0001,K=10000) {
   mu   <- start[1]
   sigma<- start[2]
   m    <- length(y)
   n    <- m + length(R)
   T1   <- sum(y)
   T2   <- sum(y^2)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out = NULL;
   while((iter<maxits)&(!converged)){
      ## z   <- stnorm( K*(n-m), mean=mu, sd=sigma, R=R)
      z   <- stnorm( (n-m), mean=mu, sd=sigma, R=R,K=K)

      U1  <- sum(z)
      U2  <- sum(z^2)

      newmu   <- (T1 + U1/K) / n 
      sigma2  <- (T2 + U2/K) / n - newmu^2;
      ## newsigma<- sqrt(abs(sigma2))
      newsigma<- sqrt(sigma2)
      # assess convergence
      converged <- (abs(newmu-mu)<eps*abs(mu)) & 
                   (abs(newsigma-sigma)<eps*abs(sigma))
      iter  <- iter+1
      mu    <- newmu
      sigma <- newsigma
      ## cat("iter =", iter, ": ", round(c(mu,sigma),3), "\n")
      out   <- rbind(out, c(mu,sigma))
   }
   cat("Done.\n")
   theta  <- c(newmu,newsigma)
   list(para=theta,iter=iter,conv=converged, out=out)
}
#---------------


#=======================================================================
# Double Exponential Distribution
#=======================================================================
# MCEM 
rc.dexp.MCEM <-
function(y,R,start,maxits=500,eps=.0001,K=10000){
   mu   <- start[1]
   sigma<- start[2]
   m    <- length(y)
   n    <- m + length(R)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out <- NULL
   while((iter<maxits)&(!converged)){
      z   <- rtdexp( K*(n-m), loc=mu, scale=sigma, R=R)

      newmu   <- median( c(rep(y,K), z) )
      y.sum   <- sum( abs(y-newmu) )
      z.sum   <- sum( abs(z-newmu) )
      newsigma<- (y.sum+z.sum/K) / n
      # assess convergence
      converged <- (abs(newmu-mu)<eps*abs(mu)) &
                   (abs(newsigma-sigma)<eps*abs(sigma))
      iter  <- iter+1
      mu    <- newmu
      sigma <- newsigma
      # cat("iter =", iter, ": ", round(c(mu,sigma),3), "\n")
      out   <- rbind(out, c(mu,sigma))
   }
   cat("Done.\n")
   theta  <- c(newmu,newsigma)
   list(para=theta,iter=iter,conv=converged,out=out)
}
#---------------------------------------------------------------------
# QEM 
rc.dexp.QEM <-
function(y,R,start,maxits=500,eps=.0001,K=1000){
   mu   <- start[1]
   sigma<- start[2]
   m    <- length(y)
   n    <- m + length(R)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out = NULL;
   while((iter<maxits)&(!converged)){
      z   <- stdexp( (n-m), loc=mu, scale=sigma, R=R, K=K)

      newmu   <- median( c(rep(y,K), z) )
      y.sum   <- sum( abs(y-newmu) )
      z.sum   <- sum( abs(z-newmu) )
      newsigma<- (y.sum+z.sum/K) / n
      # assess convergence
      converged <- (abs(newmu-mu)<eps*abs(mu)) &
                   (abs(newsigma-sigma)<eps*abs(sigma))
      iter  <- iter+1
      mu    <- newmu
      sigma <- newsigma
      # cat("iter =", iter, ": ", round(c(mu,sigma),3), "\n")
      out   <- rbind(out, c(mu,sigma))
   }
   cat("Done.\n")
   theta  <- c(newmu,newsigma)
   list(para=theta,iter=iter,conv=converged, out=out)
}

# 
#=======================================================================
# Rayleigh Distribution 
#=======================================================================
# MCEM 
rc.rayleigh.MCEM <-
function(y,R,start,maxits=500,eps=.0001,K=10000) {
   para<- start
   m    <- length(y)
   n    <- m + length(R)
   T2   <- sum(y^2)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out <- NULL
   while((iter<maxits)&(!converged)){
      z   <- rtrayleigh( K*(n-m), scale=para, R=R)
      V2  <- sum(z^2)

      newpara<- sqrt( (T2+V2/K)/(2*n) )
      # assess convergence
      converged <- (abs(newpara-para)<eps*abs(para))
      iter  <- iter+1
      para  <- newpara
      ## cat("iter =", iter, ": ", round(para,3), "\n")
      out   <- c(out, para)
   }
   cat("Done.\n")
   theta  <- newpara
   list(para=theta,iter=iter,conv=converged,out=out)
}
#-----------------------------------------------------------------------
# QEM 
rc.rayleigh.QEM <-
function(y,R,start,maxits=500,eps=.0001,K=10000) {
   para<- start
   m    <- length(y)
   n    <- m + length(R)
   T2   <- sum(y^2)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out = NULL
   while((iter<maxits)&(!converged)){
      z   <- strayleigh( (n-m), scale=para, R=R, K=K)
      V2  <- sum(z^2)

      newpara<- sqrt( (T2+V2/K)/(2*n) )
      # assess convergence
      converged <- (abs(newpara-para)<eps*abs(para))
      iter  <- iter+1
      para  <- newpara
      ## out   <- formatC(para,wid=8,dig=5,format="f")
      ## cat("iter =", iter, ": ", formatC(para,wid=8,dig=3,format="f"), "\n")
      out   <- rbind(out, para)
   }
   cat("Done.\n")
   theta  <- newpara
   list(para=theta,iter=iter,conv=converged,out=out)
}
#
#=====================================================================
# Exponential Distribution
#=====================================================================
# EM 
rc.exp.EM <-
function(y,R,start,maxits=500,eps=.0001) {
   para<- start
   m    <- length(y)
   n    <- m + length(R)
   T1   <- sum(y)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out  = NULL
   while((iter<maxits)&(!converged)){
      newpara <- sum(y)+sum(R)+(n-m)*para
      newpara <- newpara/n

      # assess convergence
      converged <- (abs(newpara-para)<eps*abs(para))
      iter  <- iter+1
      para  <- newpara
      cat("iter =", iter, ": ", round(para,3), "\n")
      out = c(out,para)
   }
   cat("..... Done .....\n\n")
   theta  <- newpara
   list(para=theta,iter=iter,conv=converged,out=out)
}
#
#-----------------------------------------------------------------------
# MCEM 
rc.exp.MCEM <-
function(y,R,start,maxits=500,eps=.0001,K=1000) {
   para<- start
   m    <- length(y)
   n    <- m + length(R)
   T1   <- sum(y)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out  = NULL
   while((iter<maxits)&(!converged)){
      newpara <- sum(y)+sum(R)-para*(n-m)*mean(log(runif(K)))
      newpara <- newpara/n

      # assess convergence
      converged <- (abs(newpara-para)<eps*abs(para))
      iter  <- iter+1
      para  <- newpara
      cat("iter =", iter, ": ", round(para,3), "\n")
      out = c(out,para)
   }
   cat("..... Done .....\n\n")
   theta  <- newpara
   list(para=theta,iter=iter,conv=converged,out=out)
}
#-----------------------------------------------------------------------
# QEM 
rc.exp.QEM <-
function(y,R,start,maxits=500,eps=.0001,K=1000) {
   para<- start
   m    <- length(y)
   n    <- m + length(R)
   T1   <- sum(y)
   iter <- 0
   converged <- FALSE
   ## cat("Performing iterations of EM...\n")
   out  = NULL
   while((iter<maxits)&(!converged)){
      newpara <- sum(y)+sum(R)-para*(n-m)*mean(log(1-ppoints(K,a=1/2)))
      newpara <- newpara/n

      # assess convergence
      converged <- (abs(newpara-para)<eps*abs(para))
      iter  <- iter+1
      para  <- newpara
      cat("iter =", iter, ": ", round(para,3), "\n")
      out = c(out,para)
   }
   cat("..... Done .....\n\n")
   theta  <- newpara
   list(para=theta,iter=iter,conv=converged,out=out)
}
#---------------------------------------------------------------------

# INTERVAL-CENSORING 
#=====================================================================
#  EM / MCEM / QEM algorithms for interval-censored
#                  from exponential, normal and Weibull
#---------------------------------------------------------------------
#
#  Required arguments:
#       X = interval observation (n x 2 matrix)
#   start = starting values
#
#  Optional arguments
#   maxits = maximum number of iterations
#      eps = convergence criterion
#
#  Value:
#      iter = how many iterations were performed
# converged = logical value indicateing whether it converged or not
#
#=====================================================================
#=====================================================================
# Exponential Distribution
#=====================================================================
ic.exp.EM <-
function (X, start=1, maxits=500, eps=1E-5) {
   lam = start
   ij  = dim(X)
   n    = ij[1]
   if ( ij[2] > 2 ) stop (" The data X should be n x 2 matrix");
   iter = 0
   converged = FALSE

   # Start the EM
   TINY = .Machine$double.eps
   HUGE = .Machine$double.xmax^0.5
   while ( (iter<maxits)&(!converged) ) {
      sumA = 0
      for ( i in 1:n ) {
          a = X[i,1]; b = X[i,2];
          if ( a <= 0 )   a = 0 ;

          if ( abs(a-b) < TINY ) {
            sumA = sumA + a ;
          } else if ( b > HUGE ) {
            sumA = sumA + a + 1/lam ;
          } else {
            sumA = sumA + 1/lam +
                (a*exp(-lam*a)-b*exp(-lam*b))/(exp(-lam*a)-exp(-lam*b))
          }
      }
      newlam    = n/sumA
      converged = (abs(newlam-lam)<eps*abs(lam))
      iter  = iter + 1
      lam   = newlam
   }
   list(lam=newlam, iter=iter, conv=converged )
}
#---------------------------------------------------------------------
# para = c(beta, lam)  #beta=shape, lam=rate
loglikelihood.Weibull <- function(para, X) {
   n = dim(X)[1]
   if ( dim(X)[2] > 2 ) stop (" The data X should be n x 2 matrix");
   TINY = .Machine$double.eps
   beta = para[1]
   lam  = para[2]
   theta = lam^(-1/beta)
   loglike = 0
   for ( i in seq_len(n) ) {
       a = X[i,1]; b = X[i,2];
       tmp = ifelse( abs(a-b) < TINY, dweibull(a, shape=beta, scale=theta, log=TRUE),
               log( pweibull(b,shape=beta,scale=theta)-pweibull(a,shape=beta,scale=theta) ) )
       loglike = loglike + tmp
   }
   return(loglike)
}
#---
negative.loglikelihood.Weibull <-  function(para, X) { -loglikelihood.Weibull(para,X) } 
#---
likelihood.Weibull <- function(para, X) { exp(loglikelihood.Weibull(para,X)) }
#---
repara <- function(para) { 
  if( is.list(para) ) {
      shape = para$beta
      scale = para$lam^(-1/shape)
  } else {
      shape = para[1]
      scale = para[2]^(-1/shape)
  }
  newpara = c(shape, scale)
  names(newpara) = c("shape", "scale") 
  return(newpara)
}
#---------------------------------------------------------------------

#=====================================================================
# Normal Distribution
#=====================================================================
ic.norm.EM <-
function(X, start=c(0,1), maxits=1000, eps=1E-5){
   mu = start[1]
   sd = start[2]
   ij = dim(X)
   n  = ij[1]
   if ( ij[2] > 2 ) stop (" The data X should be n x 2 matrix");
   iter = 0
   converged = FALSE

   # Start the EM
   TINY = .Machine$double.eps
   HUGE = .Machine$double.xmax^0.5
   while( (iter<maxits)&(!converged) ) {
      sumA = 0 ; sumB = 0
      for ( i in 1:n ) {
          a = X[i,1]; b = X[i,2];
          if ( abs(a-b) < TINY ) {
            sumA = sumA + a^2;  sumB = sumB + a   ;
          } else if ( a < -HUGE && b > HUGE ) {
            sumA = sumA + mu^2 + sd^2;  sumB = sumB + mu ;
          } else if ( b > HUGE  ) {
            z = (a-mu)/sd
            psi = dnorm(z);  PSI = pnorm(z);
            sumA = sumA + mu^2 + sd^2 + (mu+a)*sd*psi/(1-PSI) ;
            sumB = sumB + mu + sd *psi/(1-PSI) ;
          } else if ( a < -HUGE  ) {
            z = (b-mu)/sd
            psi = dnorm(z);  PSI = pnorm(z);
            sumA = sumA + mu^2 + sd^2 - (mu+b)*sd*psi/PSI ;
            sumB = sumB + mu - sd *psi/PSI ;
          } else {
            za   = (a-mu)/sd ;  zb   = (b-mu)/sd ;
            psia = dnorm(za) ;  psib = dnorm(zb) ;
            PSIa = pnorm(za) ;  PSIb = pnorm(zb) ;
            sumA = sumA + mu^2 +sd^2 -sd*((mu+b)*psib-(mu+a)*psia)/(PSIb-PSIa);
            sumB = sumB + mu - sd*(psib-psia)/(PSIb-PSIa) ;
          }
      }
      newmu = sumB/n ;
      newsd  = sqrt ( sumA/n - newmu^2 )
      # assess convergence
      converged = (abs(newmu-mu)<eps*abs(mu)) & (abs(newsd -sd)<eps*abs(sd))
      iter = iter+1
      mu   = newmu
      sd   = newsd
   }
   list( mu=newmu, sd=newsd, iter=iter, conv=converged )
}
#------------------------------------------------------------------

#=====================================================================
# Weibull Distribution
#=====================================================================
# QEM 
ic.weibull.QEM <-
function(X, beta0=1, lam0=1, maxits=1000, eps=1E-3, K=100){
   beta = beta0
   lam   = lam0
   ij = dim(X)
   n  = ij[1]
   if ( ij[2] > 2 ) stop (" The data X should be n x 2 matrix");
   iter = 0
   converged = FALSE

   # Start the EM
   TINY = .Machine$double.eps
   HUGE = .Machine$double.xmax^0.5
   xi   = (1:K-.5)/K
   EEbeta = function(beta, qik) {
          qbeta = qik^beta
          1/beta + mean(log(qik)) - sum(qbeta*log(qik))/sum(qbeta)
   }
   while( (iter<maxits)&(!converged) ) {
      qik  = array(0, dim=c(n,K) )
      for ( i in 1:n ) {
          a = X[i,1]; b = X[i,2];
          if ( abs(a-b) < TINY ) {
             qik[i,] = rep(a, K)
          } else if ( b > HUGE  ) {
             qik[i,] = (-1/lam * (log(1-xi)-lam*a^beta) )^(1/beta)
          } else if ( a <= 0  ) {
             qik[i,] = (-1/lam*log((1-xi) + xi*exp(-lam*b^beta)) )^(1/beta)
          } else {
             qik[i,] = (-1/lam*log((1-xi)*exp(-lam*a^beta)+xi*exp(-lam*b^beta)) )^(1/beta)
          }
      }
      qmax = max(qik)
      lower = 1 / mean( log(qmax) - log(qik) )
      upper =  1 / (1/lower - EEbeta(lower, qik))
      tmp = uniroot( EEbeta, interval=c(lower,upper), qik=qik)

      newbeta = tmp$root
      newlam  = 1 / mean( qik^newbeta)
      # assess convergence
      converged = (abs(newbeta-beta)<eps*abs(beta)) & (abs(newlam-lam)<eps*abs(lam))
      iter = iter+1
      beta = newbeta
      lam  = newlam
      cat(".")
   }
   cat("\n * Done * \n\n")
   list( beta=newbeta, lam=newlam, iter=iter, conv=converged )
}
#------------------------------------------------------------------


#=====================================================================
# Misc. Functions 
#=====================================================================
rtexp <-
function (n, scale = 1, R) {
   if(length(n) > 1) n <- length(n)

   scale <- rep(scale,    length = n)
       R <- rep(R,        length = n)

   R - scale*log(runif(n))
}
#--------------------------------------------------
stexp <-
function (n, scale = 1, R) {
   if(length(n) > 1) n <- length(n)

   scale <- rep(scale,    length = n)
       R <- rep(R,        length = n)

   R - scale*log(ppoints(n,a=1/2))
}
#--------------------------------------------------------------------
stnorm <-
function(n, mean, sd, R, K) {
   if(length(n) > 1) n <- length(n)

   mean <- rep(mean, length = n)
     sd <- rep(sd,   length = n)
      R <- rep(R,    length = n)

   R1    = rep(R,   rep(K,n))
   mean1 = rep(mean,rep(K,n))
   sd1   = rep(sd,  rep(K,n))
   r1    = (R1-mean1)/sd1
   U1    = rep(ppoints(K,a=1/2), n)
   x     = mean1 + sd1 * qnorm( (1-pnorm(r1))*U1 + pnorm(r1) )
   return(x)

}
#--------------------------------------------------------------------
rtnorm <-
function(n, mean, sd, R) {
   lmean<- length(mean)
   lsd  <- length(sd)
   lR   <- length(R)

   if (lmean< n)  mean <- rep(mean, length = n)
   if (lsd  < n)    sd <- rep(sd,   length = n)
   if (lR   < n)     R <- rep(R,    length = n)

   U   <-  runif(n)
   r   <- (R-mean)/sd

   mean + sd * qnorm( (1-pnorm(r))*U + pnorm(r) )
}
#=====================================================================

#--------------------------------------------------------------------
stdexp <- 
function (n, location = 0, scale = 1, R, K) { 
   if(length(n) > 1) n <- length(n)

   location <- rep(location, length = n)
      scale <- rep(scale,    length = n)
          R <- rep(R,        length = n)

   U        = rep(ppoints(K,a=1/2), n)
   R        = rep(R,       rep(K,n))
   location = rep(location,rep(K,n))
   scale    = rep(scale,   rep(K,n))
   r        = (R-location)/scale

   id1 <-  R >= location
   id2 <-  sign( U- (1-exp(r))/(2-exp(r)) )  ##- sign(U-H0)

   X1 <- R - scale*log(U)
   X2 <- location - id2*scale*log(1+id2*(1-2*U)-id2*(1-U)*exp(r))
   id1*X1  + (1-id1)*X2
}
#---------------------------------------------------------------------
rtdexp <-
function (n, location = 0, scale = 1, R) {
   lloc <- length(location)
   lsc  <- length(scale)
   lR   <- length(R)

   if (lloc < n) location <- rep(location, length = n)
   if (lsc  < n)    scale <- rep(scale,    length = n)
   if (lR   < n)        R <- rep(R,        length = n)

   U   <-  runif(n)
   r   <- (R-location)/scale
   id1 <-  R >= location
   id2 <-  sign( U- (1-exp(r))/(2-exp(r)) )  ##- sign(U-H0)

   X1 <- R - scale*log(U)
   X2 <- location - id2*scale*log(1+id2*(1-2*U)-id2*(1-U)*exp(r))
   id1*X1  + (1-id1)*X2
}
#---------------------------------------------------------------------
strayleigh <-
function (n, scale = 1, R, K) {               
   lsc  <- length(scale)
   lR   <- length(R)

   scale <- rep(scale,    length = n)
       R <- rep(R,        length = n)
   U        = rep(ppoints(K,a=1/2), n)
   R        = rep(R,       rep(K,n))
   scale    = rep(scale,   rep(K,n))
   sqrt ( R^2 - 2 * scale^2 * log(U) )
}
#---------------------------------------------------------------------
rtrayleigh <-
function (n, scale = 1, R) {
   lsc  <- length(scale)
   lR   <- length(R)

   if (lsc  < n)    scale <- rep(scale,    length = n)
   if (lR   < n)        R <- rep(R,        length = n)

   U   <-  runif(n)
   sqrt ( R^2 - 2 * scale^2 * log(U) )
}
#---------------------------------------------------------------------

