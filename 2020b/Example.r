#=====================================================================
# EM: Load-Sharing 
# Written by C. Park
#---------------------------------------------------------------------
#   <input> : Y = lifetimes between the j-th and (j+1)th failure
#---------------------------------------------------------------------
# Programmer : Park, Chanseok
# Date       : August 12, 2020
#
# Usage      : Lindley.LS.EM(Y, start=1)
#=====================================================================
Lindley.LS.EM <-   # LS = Load-Sharing
function(Y, start=1, maxits=1000L, eps=1.0E-15) {
   J = ncol(Y)
   converged = logical(J)
   para = numeric(J)
   THETA = start
   iter = numeric(J)
   if( length(start) < J ) THETA = rep_len(start,J)
   for ( j in 0L:(J-1L) ) {
       thetaj = THETA[j+1L]
       V = (J-j-1L) / (J-j)
       y = Y[,j+1L]
       while ( (iter[j+1L]<maxits)&&(!converged[j+1L]) ) {
          w = mean(y+V*(thetaj+thetaj*y+2)/thetaj/(thetaj+thetaj*y+1))
          newtheta = (1-w+sqrt( (w-1)^2+8*w ))/(2*w)
          converged[j+1L] = abs(newtheta-thetaj) < eps
          iter[j+1L] = iter[j+1L] + 1L
          thetaj = newtheta
       }
       para[j+1L] = newtheta
   }
   list(para=para, iter=iter, conv=converged)
}
#=====================================================================

#=====================================================================
# This data set is from:
#  Title: Load-sharing system model and its application to thereal data set
#  Authors: Bhupendra Singh, Puneet Kumar Gupta 
#    http://doi.org/10.1016/j.matcom.2012.02.010
# See also: http://www.stat.sc.edu/~pena/TechReports/KvamPena2003.pdf
##--------------------------------------------------------------------
Y0 = c(21.02, 24.25, 6.55, 15.35, 39.08, 16.2, 34.59, 19.1, 28.22,
32, 11.25, 17.39, 28.47, 23.42, 42.06, 28.51, 34.56, 40.33, 27.56,
9.54, 27.09, 40.36, 41.44, 32.23, 7.53, 28.34, 26.32, 30.47)

Y1 =  c(30.22, 45.54, 19.47, 16.37, 30.32, 4.16, 46.44, 38.4, 37.43,
45.52, 19.09, 25.43, 31.15, 31.28, 23.21, 33.59, 32.53, 15.35,
46.21, 36.21, 11.11, 33.21, 36.28, 8.17, 37.31, 35.58, 28.02,
40.4)

Y2 =  c(43.43, 17.19, 23.28, 25.4, 43.53, 39.52, 16.33, 20.17, 25.41,
39.11, 11.59, 22.51, 2.41, 40.03, 45.36, 16.2, 40.44, 28.33,
28.05, 28.12, 23.33, 17.04, 19.13, 41.27, 13.43, 41.48, 29.33,
42.13)
Y = cbind(Y0, Y1, Y2)
#=====================================================================



Lindley.LS.EM(Y)


