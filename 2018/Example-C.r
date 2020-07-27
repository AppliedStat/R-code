#======================================================================
# Example-C: Censored Rayleigh data
#======================================================================

#----------------------------------------------------------------------
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2018/Rs2.R")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
para <- 5; U <- runif(20)
x <- sort( round( sqrt( -2*para^2*log(U) ),  3) )

y <- c( 1.950, 2.295, 4.282, 4.339, 4.411, 4.460, 4.699, 5.319, 5.440, 5.777,
        7.485, 7.620, 8.181, 8.443,10.627)
R <- rep(10.627,5) 

#----------------------------------------------------------------------
## para <- 6.134117             ## MLE
para <- 1                    ## MCEM
#----------------------

set.seed(1)
MCEM= rc.rayleigh.MCEM(y,R, start=1, eps=0, K=1000, maxits=10 )
 QEM= rc.rayleigh.QEM (y,R, start=1, eps=0, K=1000, maxits=10 )

out1 = round( cbind(MCEM$out, QEM$out), 4)
colnames(out1) = c("beta:MCEM", "beta:QEM")
#----------------------------------------------------------------------

set.seed(1)
MCEM= rc.rayleigh.MCEM(y,R, start=10, eps=0, K=1000, maxits=10 )
 QEM= rc.rayleigh.QEM (y,R, start=10, eps=0, K=1000, maxits=10 )

out2 = round( cbind(MCEM$out, QEM$out), 4)
colnames(out2) = c("beta:MCEM", "beta:QEM")

OUT = rbind( c(1,1,10,10), cbind(out1, out2) )

rownames(OUT) = 0:10


cat("\n\n ********************* \n\n\n")
OUT

