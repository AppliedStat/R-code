#=====================================================================
#  Complete masking
#  ** EXPONENTIAL dist.
#  EM algorithm for censored and masked exponential sample
#     with J competing risks and K groups 
#  Written by C. Park
#---------------------------------------------------------------------
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
#  Date       : June 1, 2003
#
#  Usage      :
#             # data coding 
#               X1 = c(1, 2, 2, 4, 3)    
#               X2 = c(1, 2, 2, 8, 2)
#               d1 = c(1, 2, 1, 3, 0)
#               d2 = c(0, 2, 1, 3,-1)
#               X  = list(GR1=X1,GR2=X2)
#               d  = list(d1,d2)
#            expo.cr.EM(X,d,start=1)
#---------------------------------------------------------------------
expo.cr.EM <-
function(X,d, start=1, maxits=100, eps=0.0001)  {
   if (!is.list(X))  X = list(X)
   if (!is.list(d))  d = list(d)
   X.names = names(X)
   K   = length(X); 
   J   = max(unlist(d))
   if ( length(start) < J*K ) start = rep(start, l=J*K)
   lam = matrix (start, nrow=J, byrow=TRUE )
   if (is.null(X.names)) X.names = paste("G",1:length(X),sep="") 
   d.names = paste("c", 1:J, sep="")
   dimnames(lam) = list(d.names, X.names)
   newlam = lam
#--------------------------------------
   iter <- 0
   converged <- FALSE
   while ((iter<maxits)&(!converged)){
      for ( k in 1:K ) {
          nk.miss = sum( d[[k]] < 0 )
          nk.s    = sum( d[[k]] > 0 )
          nk      = sum( d[[k]] ==0 ) + nk.miss + nk.s
          sumx    = sum( X[[k]] )
          ## tmp     =  d[[k]] 
          tmp     = unique( d[[k]] )
          jj      = tmp[ tmp>0 ]
          lam.s   = sum( lam[jj,k] )
          for ( j in jj ) {
             nk.j = sum( d[[k]] == j )
             xk   = X[[k]]
             NUM  = nk - nk.miss + lam[j,k]/lam.s*nk.miss
             DEN  = sumx + 
                   (nk - nk.miss -nk.j)/lam[j,k]
             newlam[j,k] = NUM / DEN
          }
      }
      iter = iter + 1
      lam = newlam
   }
   return( newlam )
}
#==========================================================================
#            expo.cr.EM(X,d,start=1)
#            expo.cr(X,d)
# The same as expo.cr.EM, different seq is used.
expo.cr.EM2 <-
function(X,d, start=1, maxits=100, eps=0.0001)  {
   if (!is.list(X))  X = list(X)
   if (!is.list(d))  d = list(d)
   X.names = names(X)
   K   = length(X); 
   J   = max(unlist(d))
   if ( length(start) < J*K ) start = rep(start, l=J*K)
   lam = matrix (start, nrow=J, byrow=TRUE )
   if (is.null(X.names)) X.names = paste("G",1:length(X),sep="") 
   d.names = paste("c", 1:J, sep="")
   dimnames(lam) = list(d.names, X.names)
   newlam = lam
#--------------------------------------
   iter <- 0
   converged <- FALSE
   while ((iter<maxits)&(!converged)){
      for ( k in 1:K ) {
          nk.miss = sum( d[[k]] < 0 )
          nk.s    = sum( d[[k]] > 0 )
          nk      = sum( d[[k]] ==0 ) + nk.miss + nk.s
          sumx    = sum( X[[k]] )
          ## tmp     =  d[[k]] 
          tmp     = unique( d[[k]] )
          jj      = tmp[ tmp>0 ]
          lam.s   = sum( lam[jj,k] )
          for ( j in jj ) {
             nk.j = sum( d[[k]] == j )
             xk   = X[[k]]
             DEN  = sumx + (nk-nk.j)/lam[j,k] - nk.miss/lam.s
             newlam[j,k] = nk / DEN
          }
      }
      converged =  (abs(newlam-lam)<eps*abs(lam)) 
      iter = iter + 1
      lam = newlam
   }
   return( newlam )
}
#==========================================================================
#            expo.cr.EM(X1,d,start=1)
#            expo.cr(X1,d)
#==========================================================================
# With Single group, 
expo.cr.EM0 <-  
function(X,d, start=1, maxits=100, eps=0.0001)  {
   J   = max(unlist(d))
   if ( length(start) < J ) lam = rep(start,l=J)
   newlam = numeric(J)
#--------------------------------------
   iter <- 0
   converged <- FALSE

   ## idx = d
   idx     = unique( d )
   jj      = idx[ idx>0 ]
   nk.miss = sum( d < 0 )
   nk.s    = sum( d > 0 )
   nk      = sum( d ==0 ) + nk.miss + nk.s
   sumx    = sum( X )
   while ((iter<maxits)&(!converged)){
      lam.s   = sum( lam[jj] )
      for ( j in jj ) {
         nk.j = sum( d == j )
         DEN  = sumx + (nk-nk.j)/lam[j] - nk.miss/lam.s
         newlam[j] = nk / DEN
      }
      iter = iter + 1
      lam = newlam
   }
   return( newlam )
}
#=====================================================================
#  GENERAL MASKING 
#  ** EXPONENTIAL dist.
#  EM algorithm for censored and GENERALLY masked exponential sample
#     with J competing risks and a SINGLE group
#  Written by C. Park
#---------------------------------------------------------------------
#
#    <input> :
#              X =  lifetimes  
#              M =  subset of 1, 2, ..., J
#                  -1 : complete masking (equivalent to c(1:J))
#                   0 : censored
#
#              X = vector
#              M = list of vectors
#   <output> :
#
#---------------------------------------------------------------------
#
#  Programmer : Park, Chanseok
#  Date       : January 1, 2004
#
#  Usage      :
#             # data coding
#               X = c(2,5,4,1, 2, 2, 4, 3)
#               M = list(1,2,3, c(1,2), 1, 0, c(2,3), c(1,2,3) )
#                                     0 = censoring 
#            expo.cm.EM1(X,M,start=1, maxits=100, eps=1e-6)
#            expo.cm.EM2(X,M,start=1, maxits=100, eps=1e-6) 
#---------------------------------------------------------------------
expo.cm.EM1 <-  
function(X, M, start=1, maxits=200, eps=0.0001)  {
   J   = max(unlist(M))
   if ( length(start) < J ) lam = rep(start,l=J)
   newlam = numeric(J)
#--------------------------------------
   iter <- 0
   converged <- FALSE

   ## idx = d
   idx   = unique( unlist(M) )
   jj    = idx[ idx>0 ]
   nk    = length(X) 
   sumx  = sum( X )
   ## for ( i in 1:nk )  if ( M[[i]]==-1 )  M[[i]] = c(1:J)
   for ( i in 1:nk )  if ( all(M[[i]]<0) )  M[[i]] = c(1:J)
   while ((iter<maxits)&(!converged)){
      for ( j in jj ) {
          a = rep(0,nk)
          for ( i in 1:nk ) { 
              if (any(M[[i]]==j))  a[i] = lam[j] / sum(lam[M[[i]]])
          }
          newlam[j] = nk / ( sumx + sum(1-a)/lam[j] )
          ## newlam[j] = sum(a) /  sumx 
      }
      ## converged =  ( sum((newlam-lam)^2) < eps^2*sum(lam^2) ) 
      converged =  sum( abs(newlam-lam)) < eps 
      iter = iter + 1
      lam = newlam
   }
   list( para=newlam, iter=iter, conv=converged)
   ## return( newlam )
}

##
#----------------------------------------------------
## Using Albert & Baxter EM 
expo.cm.EM2 <-
function(X, M, start=1, maxits=200, eps=0.0001)  {
   J   = max(unlist(M))
   if ( length(start) < J ) lam = rep(start,l=J)
   newlam = numeric(J)
#--------------------------------------
   iter <- 0
   converged <- FALSE

   ## idx = d
   idx   = unique( unlist(M) )
   jj    = idx[ idx>0 ]
   nk    = length(X)
   sumx  = sum( X )
   ## for ( i in 1:nk )  if ( M[[i]]==-1 )  M[[i]] = c(1:J)
   for ( i in 1:nk )  if ( all(M[[i]]<0) )  M[[i]] = c(1:J)
   while ((iter<maxits)&(!converged)){
      for ( j in jj ) {
          a = rep(0,nk)
          for ( i in 1:nk ) {
              if (any(M[[i]]==j))  a[i] = lam[j] / sum(lam[M[[i]]])
          }
          ## newlam[j] = nk / ( sumx + sum(1-a)/lam[j] )
          newlam[j] = sum(a) /  sumx
      }
      ## converged =  ( sum((newlam-lam)^2) < eps^2*sum(lam^2) )
      converged =  sum( abs(newlam-lam)) < eps 
      iter = iter + 1
      lam = newlam
   }
   list( para=newlam, iter=iter, conv=converged)
   ## return( newlam )
}
#----------------------------------------------------


#=====================================================================
#
#  COMPLETE MASKING
#
#
#  ** NORMAL distribution **
#  EM algorithm for censored at the rignt and masked normal sample
#     with J competing risks and single group
#  For multiple groups, use this repeatedly for each group
#  Written by C. Park
#---------------------------------------------------------------------
#
#    <input> :
#              X =  lifetime (group, 1,2...K)
#              d =  missing(-1), censor (0) and type of cause (1,2...m)
#   <output> :
#
#---------------------------------------------------------------------
#
#  Programmer : Park, Chanseok
#  Date       : June 1, 2003
#
#  Usage      :
#             # data coding 
#               X1 = c(1, 1, 1, 2, 2, 2)    
#               X1 = c(-2, -2, -2, 4, 4, 4)    
#               d1 = c(1, 1, 1, 2, 2, 2)
#            source(" ... URL ...")
#            norm.cr.EM(X1,d1, maxits=10)
#            norm.cr.EM(-X1/100,d1, maxits=20)
#            norm.cr.EM(rev(X1),d1, maxits=50)
#---------------------------------------------------------------------
#==========================================================================
#            norm.cr.EM(X,d,start=1)
#            norm.cr(X,d)
#==========================================================================
# With Single group
norm.cr.EM <-  
function(X,d, mu0=0, sd0=1, maxits=100, eps=0.0001)  {
   nk  = length(X) 
   J   = max(unlist(d))
   if ( length(mu0) < J ) mu0 = rep(mu0, l=J)
   if ( length(sd0) < J ) sd0 = rep(sd0, l=J)
   m1 <- m2 <- array( dim=c(nk,J) )
   newmu = rep(NA,l=J)
   newsd = rep(NA,l=J)
   mu = mu0; sd = sd0
#--------------------------------------
   iter = 0
   converged = FALSE

   ## idx = d
   idx = unique( d )
   jj  = idx[ idx>0 ]

   while ((iter<maxits)&&(!converged)){
      a = array(0, dim=c(nk,J) )
      for ( i in 1:nk ) { 
          w = dnorm((X[i]-mu[jj])/sd[jj]) / 
                ( 1-pnorm((X[i]-mu[jj])/sd[jj]) )
          m1[i,jj]= mu[jj] + sd[jj] * w
          m2[i,jj]= mu[jj]^2 + sd[jj]^2 + sd[jj]*(mu[jj]+X[i])*w
          EU      = (w/sd[jj]) / sum(w/sd[jj]) 
          Ij = (d[i]==jj)
          ## I1 = (d[i]==rep(-1,length(jj)) )
          I1 = (d[i]==-1)
          a[i,jj] = Ij + I1*EU
      }
      for ( j in jj ) {
         newmu[j] = mean( a[,j]*X + (1-a[,j])*m1[,j] )
         tmp      = mean( a[,j]*(X-newmu[j])^2 +
                         (1-a[,j])*(m2[,j]-2*newmu[j]*m1[,j]+newmu[j]^2) )
         newsd[j] = sqrt(tmp)
      }
      iter = iter + 1
      ## add checking convergence 
      converged = ( (sum(abs(newmu[jj]-mu[jj]))<=eps) &&
                    (sum(abs(newsd[jj]-sd[jj]))<=eps) )
      mu = newmu; sd = newsd ; 
   }
   return( list(mu=newmu, sd=newsd, iter=iter, conv=converged) )
}
#----------------------------------------------------

#=====================================================================
#
#  GENERAL MASKING
#
#  ** NORMAL distribution ** (general MASKING)
#  EM algorithm for censored at the rignt and generally- masked
#           normal sample
#     with J competing risks and single group
#  For multiple groups, use this repeatedly for each group
#  Written by C. Park
#---------------------------------------------------------------------
#
#    <input> :
#              X =  lifetime (group, 1,2...K)
#              M =  censor (0) and type of cause (1,2...m)
#   <output> :
#
#---------------------------------------------------------------------
#
#  Programmer : Park, Chanseok
#  Date       : June 1, 2003
#
#  Usage      :
#             # data coding 
#               X1 = c(1, 1, 1, 2, 2, 2)    
#               d1 = c(1, 1, 1, 2, 0, -1)
#               M1 = list(1, 1, 1, 2, 0, 1:2)
#            source(" ... URL ...")
#            norm.cr.EM (X1,d1, maxits=10)
#            norm.cm.EM (X1,M1, maxits=10)
#==========================================================================
# With Single group
norm.cm.EM <-  
function(X, M, mu0=0, sd0=1, D, maxits=100, eps=0.0001)  {
   if (is.vector(M)) M <- as.list(M)
   nk = length(X) 

   if ( !missing(D) ) {
         if ( all(D==0) || any(D<0) || any(D>1) ) stop("D should be in (0,1]\n")
         if ( length(D) != nk ) D=rep(D,nk)
   }
   J  = max(unlist(M))
   if ( length(mu0) < J ) mu0 = rep(mu0, l=J)
   if ( length(sd0) < J ) sd0 = rep(sd0, l=J)
   m1 <- m2 <- array( dim=c(nk,J) )
   newmu = rep(NA,l=J)
   newsd = rep(NA,l=J)
   mu = mu0; sd = sd0
#--------------------------------------
   iter = 0
   converged = FALSE

   ## idx = d
   idx = unique( unlist(M) )
   jj  = idx[ idx>0 ]

   for ( i in 1:nk )  if ( all(M[[i]]<0) )  M[[i]] = 1:J
   while ((iter<maxits)&&(!converged)){
      for ( i in 1:nk ) { 
        w = dnorm((X[i]-mu[jj])/sd[jj]) / ( 1-pnorm((X[i]-mu[jj])/sd[jj]) )
        m1[i,jj]= mu[jj] + sd[jj] * w
        m2[i,jj]= mu[jj]^2 + sd[jj]^2 + sd[jj]*(mu[jj]+X[i])*w
      }
      w1    = rep(NA,l=J)
      for ( j in jj ) {
          if ( missing(D) ) {
              U = rep(0,nk)
              for ( i in 1:nk ) {
                  if (any(M[[i]]==j)) {
                     w1[jj] = dnorm((X[i]-mu[jj])/sd[jj]) /
                                ( 1-pnorm((X[i]-mu[jj])/sd[jj]) )
                     U[i] = (w1[j]/sd[j]) / sum(w1[M[[i]]] /sd[M[[i]]] )
                  }
              }
           } else {U=D}
          newmu[j] = mean( U*X + (1-U)*m1[,j] )
          ## tmp      = mean( U*(X-newmu[j])^2 +
          ##                (1-U)*(m2[,j]-2*newmu[j]*m1[,j]+newmu[j]^2) )
          tmp      = mean( U*X^2 + (1-U)*m2[,j] ) - newmu[j]^2
          newsd[j] = sqrt(tmp)
      }
      iter = iter + 1
      ## add checking convergence 
      converged = ( (sum(abs(newmu[jj]-mu[jj]))<=eps) &&
                    (sum(abs(newsd[jj]-sd[jj]))<=eps) )
      mu = newmu; sd = newsd 
   }
   return( list(mu=newmu, sd=newsd, iter=iter, conv=converged) )
}
#-----------------------------------------------------------------
