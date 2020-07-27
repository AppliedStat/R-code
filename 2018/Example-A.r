#======================================================================
# Example-A: Censored normal data
#======================================================================

#----------------------------------------------------------------------
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2018/Rs2.R")
#----------------------------------------------------------------------

 y <- c(1.613, 1.644, 1.663, 1.732, 1.740, 1.763, 1.778)
 R <- c(1.778, 1.778, 1.778)
# (MLE:mu=1.742, 0.0792)
# (MLE:mu=1.74222999, 0.07913894)

set.seed(1)
start0 = c(0,1)
EM   = rc.norm.EM  (y,R,start=start0, maxits=10, eps=0 )
MCEM = rc.norm.MCEM(y,R,start=start0, maxits=10, eps=0, K=1000 )
QEM  = rc.norm.QEM (y,R,start=start0, maxits=10, eps=0, K=1000 )


# NOTE: The MCEM results can be slightly different from the draft because
#           MCEM uses random numbers. 

out=round(cbind(EM$out[,1],MCEM$out[,1],QEM$out[,1], EM$out[,2],MCEM$out[,2],QEM$out[,2]),4)
colnames(out) = c("mu:EM", "mu:MCEM", "mu:QEM", "sigma:EM", "sigma:MCEM", "sigma:QEM")


OUT = rbind( c(0,0,0,1,1,1), out )
rownames(OUT) = 0:10

cat("\n\n ********************* \n\n\n")
OUT


