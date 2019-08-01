##========================================================================
LoadShare.expo <-
function(X) {
   J = dim(X)[2] 
   n = dim(X)[1]
   lam = numeric(J)
   names(lam) = 0:(J-1)
   for ( j in 1:J ) lam[j] =  1 / ( (J-j+1)*mean(X[,j]) )
   return( list(mle=lam, bue=(n-1)/n*lam) )
}
## LoadShare.expo (XX)

##========================================================================
LoadShare.weibull <-
function(Y) {
   J = dim(Y)[2] 
   n = dim(Y)[1]
   lam <- alpha <- numeric(J)
   names(lam) = 0:(J-1)
   names(alpha) = 0:(J-1)
  
   #------------
   EE <- function(alpha,y) {
         y.a = y^alpha; logy = log(y);
         1/alpha + mean(logy) - sum(logy*y.a) / sum(y.a) 
   }
   #------------
   for ( j in 1:J ) {
       y = Y[,j] 
       meanlog = mean(log(y))
       lower = 1 / ( log(max(y)) - meanlog )
       upper = sum( (y^lower)*log(y) ) / sum( y^lower ) - meanlog
       interval = c(lower,1/upper)
       tmp = uniroot(EE, interval=interval, y=y)
       alpha[j] = tmp$root
       lam[j] =  1 / ( (J-j+1)*mean(y^alpha[j]) )
   }
   return( list(alpha=alpha, mle=lam,  bue=(n-1)/n*lam) )
  ## return( list(alpha=alpha, lam=lam) )

}
## LoadShare.weibull (XX)

##------------------------------------------------------------------------
VAR.weibull <-  function(alpha,lam)  {
   J = min ( length(alpha), length(lam) )
   V =  as.list( rep(NA,J) )
   names(V) = 0:(J-1)
   gam = -digamma(1)
   c0 = 6/(pi^2)
   for ( j in 1:J ) { 
       tmp = log((J-j+1)*lam[j]) + gam - 1
       V11 = c0 * alpha[j]^2
       V12 = c0 * alpha[j] * lam[j] * tmp 
       V22 = c0 * lam[j]^2 * tmp^2 + lam[j]^2 
       V[[j]] = matrix( c(V11,V12,V12,V22), nrow=2 )
   }
  return(V)
}
   

