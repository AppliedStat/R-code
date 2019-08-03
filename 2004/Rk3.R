#=====================================================================
#
#  Function   : expo.cr 
#               Parameter estimation of Exponential Dist. 
#                         with competing risks with masking
#
#    <input> :
#              X =  lifetime (group, 1,2...K)
#              d =  missing(-1), censor (0) and type of cause (1,2...m)
#              X and d are list ( [[1]], ... , [[K]] )
#   <output> :
#
#---------------------------------------------------------------------
#
#  Programmer : Park, Chanseok
#  Date       : Oct. 1, 2002
#
#  Usage      :
#             # data coding 
#               X1 <- c(1, 2, 2, 4, 3)    
#               X2 <- c(1, 2, 2, 8, 2)
#               d1 <- c(1, 2, 1, 3, 0)
#               d2 <- c(0, 2, 1, 3,-1)
#               X  <- list(GR1=X1,GR2=X2)
#               d  <- list(d1,d2)
#
#             # under H1:
#               expo.cr(X,d)
#
#             # under H(00):
#               expo.cr00(X,d)
#
#             # under H(10):
#               expo.cr10(X,d)
#
#             # under H(01):
#               expo.cr01(X,d)
#               NB: apply(expo.cr(X,d),2, mean)
#
#=====================================================================
expo.cr <-
function(X,d)  {
   if (!is.list(X))  X = list(X)
   if (!is.list(d))  d = list(d)
   X.names = names(X)
   K   = length(X); 
   J   = max(unlist(d))
   lam = array ( dim=c(J,K) )
   if (is.null(X.names)) X.names = paste("G",1:length(X),sep="") 
   d.names = paste("c", 1:J, sep="")
   dimnames(lam) = list(d.names, X.names)
   for ( k in 1:K ) {
       nk.miss = sum( d[[k]] < 0 )
       nk.s    = sum( d[[k]] > 0 )
       sumx    = sum( X[[k]] )
       J = unique( unlist(d[[k]]) )
       for ( j in J[J>0] ) {
           nk.j    = sum( d[[k]] ==  j )
           lam[j,k]= (1+nk.miss/nk.s)*nk.j/sumx
       }
   }
   return (lam)
}
#--------------------------------------------------------------
expo.cr00 <-
function(X,d)  {
   d   = unlist(d)
   x   = unlist(X);
   J   = max(d)
   n.. = sum(d > 0)
   n.miss = sum(d < 0)
   return( (n..+n.miss)/sum(x)/J )
}
#--------------------------------------------------------------
expo.cr10 <-
function(X,d)  {
   x   = unlist(X)
   d   = unlist(d)
   J   = max(d)
   n.miss = sum(d < 0)
   n.. = sum(d > 0)
   n.j = numeric(J)
   names(n.j) = paste("c", 1:J, sep="")
   for ( j in 1:J ) n.j[j] <- sum( d==j )
   return( n.j*(1+n.miss/n..)/sum(x) )
}
#--------------------------------------------------------------
expo.cr01 <-
function(X,d)  {
   if (!is.list(X))  X = list(X)
   if (!is.list(d))  d = list(d)
   X.names = names(X)
   if (is.null(X.names)) X.names = paste("G",1:length(X),sep="") 
   K   = length(X); 
   J   = max( unlist(d) )
   lam = numeric(K)
   names(lam) = X.names
   for ( k in 1:K ) {
        nk.miss = sum( d[[k]] < 0 )
        nk.s    = sum( d[[k]] > 0 )
        lam[k]= (nk.miss + nk.s) / (J*sum(X[[k]]))
   }
   return(lam)
}
#--------------------------------------------------------------
weibull.cr <-
function(X,d, interval=c(.Machine$single.eps,.Machine$single.xmax))  {
   fn <- function(alpha,X,d) {
         if (!is.list(X))  X = list(X)
         if (!is.list(d))  d = list(d)
         K   = length(X);
         J   = max(unlist(d))
         tmp = 0
         for ( k in 1:K ) {
             nk.miss = sum( d[[k]] < 0 )
             nk.s    = sum( d[[k]] > 0 )
             Xa      = X[[k]]^alpha
             logX    = log( X[[k]] )
             idx     = ( d[[k]] != 0 )
             tmp = tmp + (nk.s + nk.miss)/alpha + sum( logX[idx] ) -
                   (nk.s + nk.miss) * sum( Xa*logX) / sum( Xa )
         }
         return (tmp)
    }
    OUT = uniroot(fn, interval=interval, X=X, d=d )
    return (OUT$root)
}
##  weibull.cr(X,d, interval=c(0.1,10))

#--------------------------------------------------------------
#--------------------------------------------------------------
loglike  <-
function(lam,X,d)  {
   K   = length(X); 
   J   = max (unlist(d))
   like= 0
   nk  = numeric(J)
   for (k in 1:K) {
       for (j in 1:J) nk[j] = sum(d[[k]]==j) ;
       nk.miss = sum( d[[k]] < 0 )
       lam.k = sum(lam[,k])
       sumx  = sum(X[[k]])
       like = like + sum(nk*log(lam[,k])) - 
              lam.k*sumx+nk.miss*log(lam.k)
   }
   return(like)
} 
#--------
loglike00  <- function(lam,X,d)  {
   K = length(X) 
   J = max( unlist(d) )
   loglike(matrix(rep(lam,K*J),ncol=K), X,d)
}
loglike10  <- function(lam,X,d)  {
   K = length(X) 
   loglike(matrix(rep(lam,K),ncol=K), X,d)
}
loglike01  <- function(lam,X,d)  {
   J = max( unlist(d) )
   loglike(matrix(rep(lam,rep(J,length(lam))),nrow=J), X,d)
}
#
#--------------------------------------------------------------
rmvexp <-
function (n,rate) {
 matrix( rexp( n*length(rate),  rep(rate,n) ), nrow=n, byrow=TRUE)
}



##===================================================================
# 02/14/2003
# This is test for Multivariate Exponential Dist
# Ref.: JASA (1967) Vol. 62 pg. 30 - 44
##-------------------------------------------------------------------
##
## rmvexp2 <-
## function (n, b1=1,b2=1,b3=1, b12=1,b13=1,b23=1, b123=1) {
##    t1 = rexp(n,b1);  t2 = rexp(n,b2);  t3 = rexp(n,b3);
##    t12= rexp(n,b12); t13= rexp(n,b13); t23 = rexp(n,b23);
##    t123 = rexp(n,b123) ;
##    T1 = pmin(t123, t12, t13, t1)
##    T2 = pmin(t123, t12, t23, t2)
##    T3 = pmin(t123, t13, t23, t3)
##    ## cbind(T1, T2, T3)
##    cbind(jitter(T1), jitter(T2), jitter(T3) )
## }
#--------------------------------------------------------------------
