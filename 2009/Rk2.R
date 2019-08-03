#=====================================================================
#
#  Function   : KBw.test (With Weight)
#
#    <input> : X =  list of lifetime (group, 1,2...K)
#              d =  list of censor (0) and type of cause (1,2...m)
#
#   <output> : T1 : test statistics under H01
#              T2 : test statistics under H02
#              T3 : test statistics under H03
#
#=====================================================================
KBweight.test <- 
function(X,d, m, rho=0, t=Inf)  {
   TINY<- sqrt(.Machine$double.eps)
   K   <- length(X)
   if ( missing(m) ) m <- max( unlist(d) )
   x   <- unlist(X)  
   n   <- length(x)
   ## s   <- x[x < t+TINY]    ### t 
   s   <- unique ( x[x < t+TINY] )   ### t 
   s   <- c(0,sort(s))
   ls  <- length(s)
   N   <- array( dim=c(m,K,ls) )
   Y   <- array( dim=c(K,ls-1) )
   for ( k in 1:K ) {
      x <- X[[k]]
      for (si in 2:ls) Y[k,si-1] <- sum( x > s[si]-TINY )
      for ( j in 1:m ) {
          for (si in 1:ls) N[j,k,si] <- sum( (x<s[si]+TINY)*(d[[k]]==j) )
      }
   }  
   Y.  <- apply(Y,2,sum)
   Nj. <- apply(N, c(1,3), sum)
   N.k <- apply(N, c(2,3), sum)
   N.. <- apply(N, 3, sum)
   #-------------------------
   Z1  <- array( dim=c(m,K) )
   Z2  <- array( dim=c(m,K) )
   Z3  <- array( dim=c(m,K) )
   for ( j in 1:m ) {
      for ( k in 1:K ) {
          Wjk =  Y[k,] * (Y./n)^rho
          ## Wjk =  Y[k,] 
          Yk = Y[k,]; Yk[ Yk==0 ] = 1;
          Z1[j,k] <- sum(  Wjk*(diff(N[j,k,])/Yk-diff(N..)/m/Y.) )
          Z2[j,k] <- sum(  Wjk*(diff(N[j,k,])/Yk-diff(Nj.[j,])/Y.) )
          Z3[j,k] <- sum(  Wjk*(diff(N[j,k,])/Yk-diff(N.k[k,])/m/Yk) )
      }
   }

   #----- Find S1
   mK  <- m * K
   S1  <- array( dim=c(mK,mK) )
   for ( jk in  1:mK ) {
       j <- (jk+K-1) %/% K 
       k <- (jk+K-1) %%  K + 1
       for ( j1k1 in  1:mK ) {
           j1 <- (j1k1+K-1) %/% K 
           k1 <- (j1k1+K-1) %%  K + 1
           Wjk   =  Y[k,]  * (Y./n)^rho
           Wj1k1 =  Y[k1,] * (Y./n)^rho
           ## Wjk   =  Y[k,]  
           ## Wj1k1 =  Y[k1,] 
           Yk = Y[k,]; Yk[ Yk==0 ] = 1;
           tmp   = Wjk*Wj1k1 * ((j==j1)*(k==k1)/Yk-1/m/Y.) * diff(N..)/m/Y. 
           S1[jk,j1k1] = sum( tmp ) 
       }
   }

   #----- Find S2
   S2  <- array( dim=c(m,K,K) )
   for ( j in 1:m ) {
      for (k in 1:K) {
         for (k1 in 1:K) {
             Wjk  =  Y[k,]  * (Y./n)^rho
             Wjk1 =  Y[k1,] * (Y./n)^rho
             ## Wjk  =  Y[k, ]  
             ## Wjk1 =  Y[k1,] 
             Yk = Y[k,]; Yk[ Yk==0 ] = 1;
             tmp  = Wjk*Wjk1 * ((k==k1)/Yk-1/Y.) * diff(Nj.[j,])/Y. 
             S2[j,k,k1] = sum( tmp )
         }
      }
   }

   #----- Find S3
   S3  <- array( dim=c(m,m,K) )
   for ( j in 1:m ) {
      for (j1 in 1:m) {
         for (k in 1:K) {
             Wjk  =  Y[k,] * (Y./n)^rho
             Wj1k =  Y[k,] * (Y./n)^rho
             ## Wjk  =  Y[k,] 
             ## Wj1k =  Y[k,] 
             Yk = Y[k,]; Yk[ Yk==0 ] = 1;
             tmp  = ((j==j1)-1/m) * Wjk*Wj1k /Yk * diff(N.k[k,])/m/Yk
             S3[j,j1,k] = sum(tmp)
         }
      }
   }

   #----- Find T1
   Z  <- as.vector(t(Z1))[-mK]
   S  <- S1[-mK,-mK]
   T1 <- sum( Z * ginv(S) %*% Z )
   ## T1 <- sum(Z*solve(S,Z))
  
   #----- Find T2
   T2 <- 0;
   for (j in 1:m) {
       Z <- Z2[j,-K]
       S <- S2[j,,]
       S <- S[-K,-K]
       ## T2 <- T2 + sum(Z*solve(S,Z))
       T2 <- T2 + sum( Z * ginv(S) %*% Z )
   }

   #----- Find T3
   T3 <- 0;
   for (k in 1:K) {
       Z <- Z3[-m,k]
       S <- S3[,,k]
       S <- S[-m,-m]
       ## T3 <- T3 + sum(Z*solve(S,Z))
       T3 <- T3 + sum( Z * ginv(S) %*% Z )
   }

   list( T1=T1, T2=T2, T3=T3, df=c(m*K-1, m*(K-1), (m-1)*K) )

}
##=================================================================
