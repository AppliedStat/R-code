# #########################################################################
# EM: Birnbaum-Saunders Distribution Model
# #########################################################################
BS.cm.QEM <-
function(X, M, alpha0, beta0, maxits=1000, K=1000, eps=1.0E-3)  {
   nk = length(X)
   J  = max(unlist(M))
   idx = unique( unlist(M) )
   jj  = idx[ idx>0 ]
   if (is.vector(M)) M <- as.list(M)
   for ( i in 1:nk )  if ( any(M[[i]]<0) )  M[[i]] = jj
 # 
 # Setting the inital values
   if ( !missing(alpha0) && length(alpha0) != J ) alpha0 = rep(alpha0, l=J)
   if ( !missing(beta0)  && length(beta0)  != J )  beta0 = rep(beta0,  l=J)

   if ( missing(alpha0) ) alpha0 = rep(1, l=J)
   if ( missing(beta0)  ) beta0  = rep(1, l=J)
 # end of initial value setting

   alpha.new <- beta.new <- rep(NA, l=J)
   alpha   = alpha0 ;  beta = beta0
 # Start the EM algorithm 
   iter <- 0
   converged <- FALSE
   Q = array(dim=c(nk,K))
   TINY = .Machine$double.eps
   BIG  = .Machine$double.xmax^0.5

   EE = function(beta.new, X, U, Q,qbar,qbar.star, nk, K) {
       if (beta.new<TINY) return( (1-beta.new)*BIG )
       ALPHA2 = mean(U*X+ (1-U)*qbar)/beta.new + beta.new*mean(U/X+(1-U)/qbar.star)-2
       OUT = -0.5*nk/beta.new + sum( U/(beta.new+X) ) +
             sum((1-U)*apply(1/(beta.new+Q), 1, mean)) -
             0.5*sum(U*(1/X - X/(beta.new^2)))/ALPHA2 -
             0.5*sum((1-U)*(1/qbar.star - qbar/(beta.new^2)))/ALPHA2
       return(OUT)
   }

   while ((iter<maxits)&(!converged)){
      for ( j in jj ) {
          U = numeric(nk)
          for ( i in 1:nk ) {
              if (any(M[[i]]==j)) {
                 tmp  = hBS(x=X[i], alpha[M[[i]]], beta[M[[i]]])
                 U[i] = hBS(x=X[i], alpha[j], beta[j]) / sum(tmp)
              }
          }
          for ( i in seq_len(nk) ) {
              Q[i, ] = stBS(n=K, alpha[j], beta[j], R=X[i])
          }
          qbar = apply(Q,1,mean)
          qbar.star = 1/apply(1/Q,1,mean)

          lower = min(TINY,1/mean(1/X)); upper = mean(X)
          tmp = uniroot(EE, interval=c(lower,upper), extendInt ="yes",
                        X=X, U=U, Q=Q, qbar=qbar, qbar.star=qbar.star, nk=nk,K=K)
          beta.new[j] = tmp$root

          alpha.new[j]= sqrt(mean(U*X+ (1-U)*qbar)/beta.new[j] +beta.new[j]*mean(U/X+(1-U)/qbar.star)-2)
      }
      iter = iter + 1
      conv1 = all ( abs(alpha.new[jj]-alpha[jj]) < eps*abs(alpha.new[jj]) )
      conv2 = all ( abs(beta.new[jj] - beta[jj]) < eps*abs(beta.new[jj]) )
      converged = conv1 && conv2
      alpha = alpha.new; beta = beta.new ;
      cat(".")
    }
    cat("\n * Done (BS) *\n\n")
    list ( alpha=alpha.new, beta=beta.new, iter=iter, conv=converged )
}

#====================================================================
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
   az = alpha*qnorm(p)
   .25*beta*(az + sqrt(az*az+4))^2
}
hBS <- function(x, alpha=1, beta=1)  {
    y = sqrt(x/beta)
    dBS(x=x,alpha=alpha,beta=beta) / pnorm(y-1/y,sd=alpha,lower.tail=FALSE)
}
#--------------------------------------------------------------------
qtBS <- function(p, alpha=1, beta=1, R=0)  {   # truncated at R
   gamp = alpha*qnorm(p+(1-p)*pnorm((sqrt(R/beta)-sqrt(beta/R))/alpha))/sqrt(2)
   beta * (1+gamp^2 + gamp*sqrt(gamp^2+2))
}
stBS <- function(n, alpha=1, beta=1, R=0)  {   # truncated at R
   p = ppoints(n)
   gamp = alpha*qnorm(p+(1-p)*pnorm((sqrt(R/beta)-sqrt(beta/R))/alpha))/sqrt(2)
   return(beta*(1+gamp^2 + gamp*sqrt(gamp^2+2)))
}

#--------------------------------------------------------------------
SBS <- function(x, alpha, beta)  {   # Survival function of BS
   J = length(alpha) 
   S = rep(1, length(x)) 
   for (j in 1:J) {
       y = sqrt(x/beta[j])
       S = S * pnorm(y-1/y,sd=alpha[j],lower.tail=FALSE)
   }
   return(S)
}


