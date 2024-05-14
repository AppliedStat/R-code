# ------------------------------
# EM: inverse Weibull Distribution Model
# -------------------------------
invweibull.cr.EM <-
function(X, M, beta0, theta0, maxits=100, eps=1.0E-3)  {
   nk = length(X)
   J  = max(unlist(M))
   idx = unique( unlist(M) )
   jj  = idx[ idx>0 ]
   if (is.vector(M)) M <- as.list(M)
   for ( i in 1:nk )  if ( any(M[[i]]<0) )  M[[i]] = jj
 #
 # Setting the inital values
   if ( !missing(beta0)  && length(beta0)  !=J )  beta0 = rep(beta0, l=J)
   if ( !missing(theta0) && length(theta0) !=J ) theta0 = rep(theta0,l=J)
   X0 = as.list(NULL); length(X0) = J
   n1 = numeric(J)
   if ( missing(beta0) || missing(theta0) ) {
      for ( i in 1:nk ) {
          idx = M[[i]]
          if (length(idx) == 1) {
              if ( idx > 0 ) {
                      X0[[idx]] = c(X0[[idx]], X[i])
                       n1[idx]  = n1[idx] + 1
              } else if (idx == 0) for (j in jj) X0[[j]] = c(X0[[j]], X[i])
          }
      }
    }
    if  (  missing(beta0) ||  missing(theta0) ) {
        beta00 = theta00 = numeric(J)
        for ( j in jj ) {
           if ( is.null( X0[[j]] ) ) {
               tmp = invweibull.MLE ( X )
           } else {
               tmp = invweibull.MLE ( X0[[j]] )
           }
           beta00[j]  = tmp$shape
           theta00[j] = tmp$scale
        }
   }
   if ( missing(beta0) )   beta0 = beta00
   if ( missing(theta0) ) theta0 = theta00
   # end of initial value setting

   newtheta <- newbeta <- rep(NA, l=J)
   theta = theta0; beta = beta0
 #

 # Start the EM algorithm 
   iter <- 0
   converged <- FALSE

   EE2 = function(p,x,U) {  # p = parameters
       TINY = .Machine$double.eps
       HUGE = .Machine$double.xmax^0.25
       if( p[1] < TINY ) return( HUGE*(1-p[1]) )
       if( p[2] < TINY ) return( HUGE*(1-p[2]) )

       beta = p[1]
       theta= p[2]

       logS = pinvweibull(x, shape=beta, scale=theta, lower.tail=FALSE, log.p=TRUE )
     # tmp1 = U*(log(beta)-log(theta)+(beta+1)*(log(theta)-log(x)) -(theta/x)^beta) -
     #        log(1-exp(-(theta/x)^beta))
       tmp1 = U*(dinvweibull(x, shape=beta, scale=theta, log=TRUE ) - logS)

     # tmp2 = log(1-exp(-(theta/x)^beta))
       tmp2 = logS
       return( -sum(tmp1+tmp2) )
   }

   while ( (iter<maxits)&(!converged) ){
      for ( j in jj ) {
          U = rep(0,nk)
          for ( i in 1:nk ) {
              if (any(M[[i]]==j)) {
                 tmp1 = dinvweibull( X[i], shape=beta[M[[i]]], scale=theta[M[[i]]], log=TRUE )
                 tmp2 = pinvweibull( X[i], shape=beta[M[[i]]], scale=theta[M[[i]]], lower.tail=FALSE, log.p=TRUE )
                 deno= sum( exp(tmp1-tmp2) )

                 tmp1 = dinvweibull( X[i], shape=beta[j], scale=theta[j], log=TRUE )
                 tmp2 = pinvweibull( X[i], shape=beta[j], scale=theta[j], lower.tail=FALSE, log.p=TRUE )
                 num  = exp( tmp1-tmp2 )
                 U[i] = num / deno
              }
          }

          tmp = nlm(EE2, c(beta[j],theta[j]), x=X, U=U)
          newbeta[j]  = tmp$estimate[1]
          newtheta[j] = tmp$estimate[2]
      }
      iter = iter + 1
      conv1 = all ( abs(newbeta[jj]-beta[jj]) < eps*abs(newbeta[jj]) )
      conv2 = all ( abs(newtheta[jj]  -  theta[jj]) < eps*abs(newtheta[jj]) )
      converged = conv1 && conv2
      theta = newtheta ;  beta = newbeta
   }
   list ( beta=newbeta, theta=newtheta, iter=iter, conv=converged )
}



#=================================================================================
invweibull.MLE = function (x, interval, tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0) {
    y = 1/x
    if (missing(interval)) {
        meanlog = mean(log(y))
        lower = 1/(log(max(y)) - meanlog)
        upper.rev = sum((y^lower) * log(y))/sum(y^lower) - meanlog
        upper = 1/upper.rev
        if (is.nan(upper))
            upper = .Machine$double.xmax^0.2
        interval = c(lower, upper)
    }
    EEweibull = function(alpha, y) {
        TINY = .Machine$double.neg.eps
        yalpha = y^alpha
        if (sum(yalpha) < TINY)
            return(mean(log(y)) - 1/alpha - mean(log(y)))
        sum(log(y) * (yalpha))/sum(yalpha) - 1/alpha - mean(log(y))
    }
    tmp = uniroot(EEweibull, interval = interval, y = y, tol = tol,
        extendInt = "upX", maxiter = maxiter, trace = trace)
    alpha = tmp$root
    beta = 1/mean(y^alpha)^(1/alpha)
    list(shape = alpha, scale = beta)
}
#---------------------------------------------------------------------------------

pinvweibull = function (q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    k <- max(lq <- length(q), lshape <- length(shape), lscale <- length(scale))
    if (lq < k)
        q <- rep(q, length = k)
    if (lshape < k)
        shape <- rep(shape, length = k)
    if (lscale < k)
        scale <- rep(scale, length = k)
    id = which(q > 0)
    p = numeric(length(q))
    p[id] = pweibull(1/q[id], shape = shape[id], scale = 1/scale[id],
        lower.tail = !lower.tail, log.p = log.p)
    if (!is.null(Names <- names(q)))
        names(p) = rep(Names, length = k)
    return(p)
}
#---------------------------------------------------------------------------------
dinvweibull = function (x, shape, scale = 1, log = FALSE) {
    k <- max(lx <- length(x), lshape <- length(shape), lscale <- length(scale))
    if (lx < k)
        x <- rep(x, length = k)
    if (lshape < k)
        shape <- rep(shape, length = k)
    if (lscale < k)
        scale <- rep(scale, length = k)
    logd = numeric(k)
    id = (x > 0)
    logd[id] = log(shape[id]) - log(scale[id]) - (shape[id] +
        1) * log(x[id]/scale[id]) - 1/(x[id]/scale[id])^shape[id]
    logd[!id] = -Inf
    if (!is.null(Names <- names(x)))
        names(logd) = rep(Names, length = k)
    if (log) {
        return(logd)
    }
    else {
        return(exp(logd))
    }
}
#---------------------------------------------------------------------------------

#
# MLE of  Weibull 
#
weibull.MLE <- function(x, interval) {
  if ( any (x <= 0) ) stop("The data should be positive")
  if (missing(interval)) {
     meanlog = mean(log(x))
     lower = 1 / ( log(max(x)) - meanlog )
     upper = sum( (x^lower)*log(x) ) / sum( x^lower ) - meanlog
     interval = c(lower,1/upper)
  }
  EE = function(alpha,x) {
     xalpha = x^alpha
     sum(log(x)*(xalpha)) / sum(xalpha) - 1/alpha - mean(log(x))
  }
  tmp = uniroot(EE, interval=interval, x=x)
  alpha = tmp$root
  list ( alpha=alpha, lam=1/mean(x^alpha) )
}

# -------------------------------
# EM: Weibull Distribution Model
# -------------------------------
weibull.cm.EM <-
function(X, M, alpha0, lam0, maxits=100, eps=1.0E-3)  {
   nk = length(X)
   J  = max(unlist(M))
   idx = unique( unlist(M) )
   jj  = idx[ idx>0 ]
   if (is.vector(M)) M <- as.list(M)
   for ( i in 1:nk )  if ( any(M[[i]]<0) )  M[[i]] = jj 
 #
 # Setting the inital values
   if ( !missing(alpha0)&&length(alpha0)!=J ) alpha0 = rep(alpha0,l=J)
   if ( !missing(lam0)  && length(lam0) !=J )   lam0 = rep(lam0,  l=J)
   X0 = as.list(NULL); length(X0) = J
   n1 = numeric(J)
   if ( missing(alpha0) || missing(lam0) ) {
      for ( i in 1:nk ) {
          idx = M[[i]]
          if (length(idx) == 1) {
              if ( idx > 0 ) {
                      X0[[idx]] = c(X0[[idx]], X[i])
                       n1[idx]  = n1[idx] + 1 
              } else if (idx == 0) for (j in jj) X0[[j]] = c(X0[[j]], X[i])
          }
      }
    }
    if  (  missing(alpha0) ||  missing(lam0) ) {
        alpha00 = lam00 = numeric(J)
        for ( j in jj ) {
           if ( is.null( X0[[j]] ) ) {
                    tmp = weibull.MLE ( X )
           } else { tmp = weibull.MLE ( X0[[j]] ) }
            alpha00[j] = tmp$alpha 
              lam00[j] = tmp$lam 
        }
   }
   if (  missing(alpha0) ) alpha0 = alpha00
   if (  missing(lam0) )     lam0 = lam00
   # end of initial value setting

   newlam <- newalpha <- rep(NA, l=J)
   lam   = lam0 ;  alpha = alpha0
 #
 # Start the EM algorithm 
   iter <- 0
   converged <- FALSE
   sumx  = sum( X )
   EE2 = function(alpha,x,U) {
       xalpha = x^alpha
       sumU = sum(U)
       sumU*sum(xalpha*log(x))/sum(xalpha)-sumU/alpha-sum(U*log(x))
   }
   while ((iter<maxits)&(!converged)){
      for ( j in jj ) {
          U = rep(0,nk)
          for ( i in 1:nk ) {
              if (any(M[[i]]==j)) {
                 tmp  = alpha[M[[i]]]*lam[M[[i]]]*X[i]^(alpha[M[[i]]]-1) 
                 U[i] = alpha[j]*lam[j]*X[i]^(alpha[j]-1) / sum(tmp)
              }
          }
           sumU = sum(U)
          lower = sumU / sum( U*(log(max(X))-log(X)) ) 
            tmp = sumU*sum(X^lower * log(X))/sum(X^lower) - sum(U*log(X))
          upper = max( abs(sumU/tmp), lower )
          while ( EE2(lower,X,U)*EE2(upper,X,U) > 0 ) {
             upper = 2 * upper
          }
          tmp = uniroot(EE2, interval=c(lower,upper),x=X, U=U)
          newalpha[j] = tmp$root 
          newlam[j]   = sumU / sum(X^newalpha[j])
      }
      iter = iter + 1
      conv1 = all ( abs(newalpha[jj]-alpha[jj]) < eps*abs(newalpha[jj]) )
      conv2 = all ( abs(newlam[jj]  -  lam[jj]) < eps*abs(newlam[jj]) )
      converged = conv1 && conv2 
      lam = newlam ;  alpha = newalpha 
   }
   list ( lam=newlam, alpha=newalpha, iter=iter, conv=converged )
}
