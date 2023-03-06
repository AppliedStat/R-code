#============================================================================
source ("https://raw.githubusercontent.com/AppliedStat/R-code/master/2024b/fnmw41.R")
# If the above does not work, then just cut and past the R codes from 
#     https://github.com/AppliedStat/R-code/blob/master/2024b/fnmw41.R

##  Data from Table C1: No coating 
##  Harwell, M. (1995).  Microbond Tests for Ribbon Fibers, 
##  M.S. Thesis, Dept of Chemical Engineering, Clemson University. 
##  D = Debonding,  B = Fiber break
## 
MODEC1 = c(
 "B", "D", "B", "D", "D", "D", "D", "D", "B", "D",
 "D", "D", "D", "D", "D", "D", "D", "D", "D", "D",
 "D", "D", "D", "D", "D", "D", "D", "D", "D", "D",
 "D", "D", "D", "D", "D", "D", "D", "D" )
MicrobondC1 = c(
0.198, 0.212, 0.33 , 0.321, 0.371, 0.216, 0.285, 0.259, 0.356, 0.338,
0.268, 0.219, 0.211, 0.206, 0.253, 0.264, 0.266, 0.247, 0.234, 0.285,
0.32 , 0.275, 0.298, 0.334, 0.295, 0.281, 0.222, 0.199, 0.283, 0.217,
0.282, 0.246, 0.181, 0.183, 0.283, 0.244, 0.224, 0.286 )
modeC1 <- numeric( length(MODEC1) )
modeC1[ MODEC1 == "D" ] <- 1
modeC1[ MODEC1 == "B" ] <- 2

#============================================================================
data1 = MicrobondC1[modeC1==1] ; data2 = MicrobondC1[modeC1==2] ; data3 = MicrobondC1[modeC1==3] ; 
# plot(  log(sort(data1)),  log(-log(1-ppoints(data1))), pch=2, xlab="logx", ylab="loglog"  )
# plot(  log(sort(data2)),  log(-log(1-ppoints(data2))), pch=2, xlab="logx", ylab="loglog"  )
# plot(  log(sort(data3)),  log(-log(1-ppoints(data3))), pch=2, xlab="logx", ylab="loglog"  )
##==========================================================================
idx  = modeC1
data.all = MicrobondC1 
idx.sort = order( MicrobondC1 )
idx1 = ( idx[idx.sort]==1)
idx2 = ( idx[idx.sort]==2)
data.sort = data.all[idx.sort]
Fall     = ppoints( data.all )
##-----------------------------------------------------------------------
 X = sort(data.all) ; d = idx[ order(data.all) ]
##-----------------------------------------------------------------------

##-----------------------------------------------------------------------
## BS Model 
##-----------------------------------------------------------------------
 SBS <- function(q, alpha, beta)  {   # Survival function of BS
    y1 = sqrt(q/beta[1])
    y2 = sqrt(q/beta[2])
    pnorm(y1-1/y1,sd=alpha[1],lower.tail=FALSE) * pnorm(y2-1/y2,sd=alpha[2],lower.tail=FALSE)
 }
 para.BS = BS.cm.QEM(X,d, maxits=500, eps=1.0E-5, K=1000)
 FBS = 1-SBS(X, alpha=para.BS$alpha, beta=para.BS$beta)
 MSE = mean( (FBS-Fall)^2 )
##=========================================================================


##=========================================================================
cat("\n\n ====================================\n")
cat(" MSE =", MSE)

cat("\n\n ====================================\n")
cat(" BS parameter\n")
print(para.BS)


