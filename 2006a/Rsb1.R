# ==============================================================================
# File name : fnsb1.R 
# Authors   : Chanseok Park
# Version   : 1.0,   June 15, 2004 #  based on R Version 1.9.0 (2004-04-12)
#           : 1.1,   November 25, 2016
# ==============================================================================

#-----------------------------------------------------------------------
norm.ic.EM <-
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
   HUGE = .Machine$double.xmax^0.2 
   colnames(X) = NULL
   while( (iter<maxits)&(!converged) ) {
      ST1 = 0 ; ST2 = 0
      for ( i in 1:n ) {
          a = X[i,1]; b = X[i,2];
          if ( abs(a-b) < TINY ) {
            ST1 = ST1 + a^2;  ST2 = ST2 + a   ;
          } else if ( a < -HUGE && b > HUGE ) {
            ST1 = ST1 + mu^2 + sd^2;  ST2 = ST2 + mu ;
          } else if ( b > HUGE  ) {
            z = (a-mu)/sd
            logpsi  = dnorm( z, log=TRUE);
            logPSI1 = pnorm(-z, log.p=TRUE);
            tmp     = exp(logpsi-logPSI1)
            ST1 = ST1 + mu^2 + sd^2 + (mu+a)*sd*tmp;
            ST2 = ST2 + mu + sd*tmp;
          } else if ( a < -HUGE  ) {
            z = (b-mu)/sd
            psi = dnorm(z);  PSI = pnorm(z);
            logpsi = dnorm(z, log=TRUE);
            logPSI = pnorm(z, log.p=TRUE);
            tmp    = exp(logpsi-logPSI)
            ST1 = ST1 + mu^2 + sd^2 - (mu+b)*sd*tmp ;
            ST2 = ST2 + mu - sd*tmp ;
          } else {
            za   = (a-mu)/sd ;  zb   = (b-mu)/sd ;
            psia = dnorm(za) ;  psib = dnorm(zb) ;
            PSIa = pnorm(za) ;  PSIb = pnorm(zb) ;
            ST1 = ST1 + mu^2 + sd^2 -sd*((mu+b)*psib-(mu+a)*psia)/(PSIb-PSIa);
            ST2 = ST2 + mu - sd*(psib-psia)/(PSIb-PSIa) ;
          }
      }
      newmu = ST2/n ;
      newsd  = sqrt ( ST1/n - newmu^2 )

      # assess convergence
      converged = (abs(newmu-mu)<eps*abs(mu)) & (abs(newsd -sd)<eps*abs(sd))
      iter = iter+1
      mu   = newmu
      sd   = newsd 
   }
   list( mu=newmu, sd=newsd, iter=iter, conv=converged )
}
#------------------------------------------------------------------


#---------------------------------------------------------------------
loglike <-
function(X, mu, sd ){
   TINY = .Machine$double.eps
   ij = dim(X)
   n  = ij[1]
   if ( ij[2] > 2 ) stop (" The data X should be n x 2 matrix");

   # Start the EM
   ST1 = 0 ; 
   for ( i in 1:n ) {
       a = X[i,1]; b = X[i,2];
       if ( abs(a-b) < TINY ) {
         ST1 = ST1 + dnorm(a, mean=mu, sd=sd, log=TRUE)
       } else {
         ST1 = ST1 + log(pnorm(b, mean=mu, sd=sd)-pnorm(a, mean=mu, sd=sd)) 
       }
   } 
   return(ST1) 
}
#------------------------------------------------------------------


